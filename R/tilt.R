#----------------------------------------------#
# Exponential-tilt (tilted cell moment) engine
# for feglm -- all families
#----------------------------------------------#

#
# The method: within the IRLS/Fisher-scoring iterations of feglm, the weighted
# least-squares problem solved at each iteration,
#
#     (A' W A) beta = A' W z,
#
# with A the N x q model matrix, W = diag(w_i) the IRLS weights
# w_i = weights_i * (d mu/d eta)_i^2 / V(mu_i) and z the working response,
# decomposes exactly when a subset of the columns of A is constant within
# the J cells of the cross-classification of the categorical covariates.
# Writing A = [X, D] with X the p individual-level columns and D the dummy
# expansion of the k cell-level columns (D_i = W_cell[j(i), ]), the normal
# equations have blocks
#
#     X'WX = Q                    (p x p Gram matrix, one O(N p^2) pass)
#     X'WD = M1' W_cell           (M1_j = sum_{i in j} w_i x_i, tilted first moments)
#     D'WD = W_cell' diag(M0) W_cell   (M0_j = sum_{i in j} w_i, tilted mass)
#     X'Wz, D'Wz                  (from cell sums of w_i z_i)
#
# The w_i-weighting is the exponential tilt of the within-cell covariate
# measure: for Poisson w_i = mu_i and the cell moments are the moments of the
# Esscher transform of the within-cell exposure measure; for other families
# the tilt is the Fisher/IRLS weight (variance tilt for logit, identically 1
# for Gamma-log, etc.). Hence one O(N(p^2 + p + 1)) streaming pass over the
# microdata plus O(J(k^2 + kp)) cell algebra replaces the O(N(p+k)^2)
# weighted cross-product of the full design -- while producing the *same*
# normal equations, so the IRLS iterates (and the final estimates, standard
# errors, and likelihood) are identical to the standard feglm algorithm up
# to floating-point summation order. This works for every feglm family,
# canonical or not, because the decomposition only touches the WLS step that
# all families share.
#

# Builds the cell structure from the model matrix.
# Cell-level columns are low-cardinality columns (intercept, factor dummies
# and their interactions, binary or otherwise discrete regressors): any such
# set is constant within the cells defined by its own joint row patterns, so
# treating it at the cell level is exact by construction. The patterns are
# indexed by cpp_to_index, which hashes the raw values with exact equality
# confirmation -- no arithmetic on the keys, hence no collision risk. The
# low-cardinality screen runs on a row sample only: it is a selection
# device, not a correctness requirement, because a rich column that slips
# through simply multiplies the cell count until the max_cells cap pushes
# it back to the individual level (still exact).
tilt_setup = function(X, max_cells = 1e5){

  n = nrow(X)
  q = ncol(X)

  if(is.null(n) || n == 0 || q == 0) return(NULL)

  max_cells = min(max_cells, ceiling(n / 2))

  # screening on a row sample
  rows_sample = if(n > 5000) as.integer(seq(1, n, length.out = 5000)) else seq_len(n)
  cand = integer(0)
  for(j in seq_len(q)){
    v = X[rows_sample, j]
    if(!anyNA(v) && length(unique(v)) <= 2) cand = c(cand, j)
  }

  if(length(cand) < 2) return(NULL)

  # joint exact indexing of the candidate columns: chunks of 8 columns are
  # indexed by cpp_to_index (fast lookup-table path), and the resulting
  # dense integer ids are packed into exact double keys -- products of
  # ranges are kept below 2^50, so the packing arithmetic is exact integer
  # arithmetic and collision-free -- which are indexed jointly at the end
  chunk_idx = list()
  chunk_J = integer(0)
  for(k0 in seq(1, length(cand), by = 8)){
    chunk = cand[seq(k0, min(k0 + 7, length(cand)))]
    info = cpp_to_index(lapply(chunk, function(j) X[, j]))
    chunk_idx[[length(chunk_idx) + 1]] = info$index
    chunk_J[length(chunk_J) + 1] = length(info$first_obs)
  }

  keys = list()
  keyv = NULL
  mult = 1
  for(t in seq_along(chunk_idx)){
    Jc = as.numeric(chunk_J[t])
    # flush before the packing arithmetic could become inexact
    if(!is.null(keyv) && mult * Jc > 2^53){
      keys[[length(keys) + 1]] = keyv
      keyv = NULL
    }
    if(is.null(keyv)){
      keyv = as.numeric(chunk_idx[[t]])
      mult = Jc
    } else {
      keyv = keyv + mult * (chunk_idx[[t]] - 1)
      mult = mult * Jc
    }
  }
  keys[[length(keys) + 1]] = keyv

  info = cpp_to_index(keys)

  if(length(info$first_obs) <= max_cells){
    id = info$index
    J = length(info$first_obs)
    first_rows = info$first_obs
    ccols = cand
  } else {
    # greedy fallback: add columns one at a time under the cell cap, so
    # that only the cap-breaking columns stay at the individual level
    id = rep(1L, n)
    J = 1
    first_rows = NULL
    ccols = integer(0)
    for(j in cand){
      info = cpp_to_index(list(id, X[, j]))
      if(length(info$first_obs) > max_cells) next
      id = info$index
      J = length(info$first_obs)
      first_rows = info$first_obs
      ccols = c(ccols, j)
    }
  }

  # tilting pays off only if there are cell-level columns beyond the
  # intercept and the cells actually compress the data
  if(length(ccols) < 2 || J >= n / 2) return(NULL)

  icols = setdiff(seq_len(q), ccols)

  # cell ids are dense in first-occurrence order, so the rows of W_cell are
  # aligned with the ids 1:J; any selected column is constant within cells
  # by construction (the cells refine the selected columns' own patterns)
  W_cell = X[first_rows, ccols, drop = FALSE]

  list(cells = id, J = J, W = W_cell,
       ccols = ccols, icols = icols,
       Xi = X[, icols, drop = FALSE])
}

# One weighted least-squares solve via tilted cell moments. Returns an object
# with the same fields as the feols(fromGLM = TRUE) result that feglm.fit
# consumes, so it can be swapped in place of the feols call inside the IRLS
# loop.
tilt_wls = function(z, X, w, tilt, collin.tol, nthreads){

  cells = tilt$cells
  W_cell = tilt$W
  Xi = tilt$Xi
  icols = tilt$icols
  ccols = tilt$ccols

  p = length(icols)
  kc = length(ccols)
  q = p + kc

  wz = w * z

  # single grouped pass: tilted mass M0, residual sums Rz, first moments M1
  if(p > 0){
    S = rowsum(cbind(w, wz, Xi * w), cells)
  } else {
    S = rowsum(cbind(w, wz), cells)
  }
  M0 = S[, 1]
  Rz = S[, 2]

  xwx = matrix(0, q, q)
  xwy = numeric(q)

  xwx[ccols, ccols] = cpp_crossprod(W_cell, M0, nthreads)
  xwy[ccols] = crossprod(W_cell, Rz)

  if(p > 0){
    M1 = S[, -(1:2), drop = FALSE]
    Hgd = crossprod(M1, W_cell)

    xwx[icols, icols] = cpp_crossprod(Xi, w, nthreads)
    xwx[icols, ccols] = Hgd
    xwx[ccols, icols] = t(Hgd)
    xwy[icols] = crossprod(Xi, wz)
  }

  # non-finite values (e.g. exploding working response): let feglm's
  # divergence handling take over
  if(anyNA(xwx) || anyNA(xwy)){
    return(list(coefficients = NA_real_, multicol = FALSE))
  }

  info_inv = cpp_cholesky(xwx, collin.tol, nthreads)

  if(!is.null(info_inv$all_removed)){
    stopi("All variables are collinear with each other. ",
          "Without doubt, your model is misspecified.")
  }

  is_excluded = info_inv$id_excl
  multicol = any(is_excluded)

  if(multicol){
    beta = as.vector(info_inv$XtX_inv %*% xwy[!is_excluded])
    names(beta) = colnames(X)[!is_excluded]
  } else {
    beta = as.vector(info_inv$XtX_inv %*% xwy)
    names(beta) = colnames(X)
  }

  # fitted values from the cell structure: O(N p + J k) instead of O(N (p+k))
  beta_full = numeric(q)
  beta_full[!is_excluded] = beta
  fitted.values = as.vector(W_cell %*% beta_full[ccols])[cells]
  if(p > 0){
    fitted.values = fitted.values + cpp_xbeta(Xi, beta_full[icols], nthreads)
  }

  list(coefficients = beta, fitted.values = fitted.values,
       residuals = z - fitted.values,
       multicol = multicol, is_excluded = is_excluded,
       collin.min_norm = info_inv$min_norm,
       means = NULL, X_demean = X, xwx = xwx)
}
