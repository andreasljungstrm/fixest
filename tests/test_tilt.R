#----------------------------------------------------------------------------#
# Test of the exponential-tilt (tilted cell moment) algorithm of feglm       #
#                                                                            #
# Two parts:                                                                 #
#  1) EXACTNESS: for every feglm family (canonical and non-canonical links,  #
#     with and without weights/offsets, under collinearity, and in           #
#     adversarial designs), feglm(..., tilt = TRUE) must reproduce the       #
#     original algorithm's coefficients, standard errors, vcov, deviance     #
#     and log-likelihood to floating-point accuracy. The tilt is an exact    #
#     reorganization of the IRLS arithmetic, not an approximation, so the    #
#     tolerance below is orders of magnitude tighter than the glm            #
#     convergence tolerance itself.                                          #
#  2) GAINS: what is gained is per-iteration work. The tilted decomposition  #
#     computes each weighted least-squares step in O(N p^2) + O(J k^2)       #
#     instead of O(N (p+k)^2), with p the continuous-involving columns,      #
#     k the dummy columns and J the number of cells, and it reuses the       #
#     assembled information matrix so the final O(N q^2) Hessian             #
#     cross-product is avoided as well. The script fits the same models      #
#     with tilt = FALSE / TRUE and prints total time, per-iteration time,    #
#     and the maximum coefficient/standard-error discrepancies.              #
#----------------------------------------------------------------------------#

library(fixest)
setFixest_notes(FALSE)

tol_exact = 1e-8

max_diff = function(a, b) max(abs(a - b))

compare_fits = function(fml, family, data, label, tol = tol_exact, ...){
  orig = feglm(fml, data, family = family, ...)
  tilt = feglm(fml, data, family = family, tilt = TRUE, ...)

  d_coef = max_diff(coef(orig), coef(tilt))
  d_se   = max_diff(se(orig), se(tilt))
  d_dev  = abs(deviance(orig) - deviance(tilt)) / (1 + abs(deviance(orig)))
  ll_o = logLik(orig)
  d_ll = if(is.na(ll_o)) 0 else abs(ll_o - logLik(tilt)) / (1 + abs(ll_o))

  cat(sprintf("  %-34s dcoef = %.1e  dse = %.1e  iters = %d/%d\n",
              label, d_coef, d_se, orig$iterations, tilt$iterations))

  stopifnot(d_coef < tol, d_se < tol, d_dev < tol, d_ll < tol,
            orig$iterations == tilt$iterations)
  invisible(list(orig = orig, tilt = tilt))
}

#
# Part 1: exactness across all feglm families ####
#

cat("PART 1: tilt = TRUE reproduces the original algorithm exactly\n\n")

set.seed(1)
n = 10000
base = data.frame(x1 = rnorm(n), x2 = rnorm(n),
                  f1 = factor(sample(letters[1:12], n, TRUE)),
                  f2 = factor(sample(1:6, n, TRUE)),
                  w  = runif(n, 0.5, 2),
                  off = runif(n))
eta = 0.3*base$x1 - 0.2*base$x2 + 0.05*as.integer(base$f1) + 0.05*as.integer(base$f2)
base$y_pois = rpois(n, exp(eta - 0.3))
base$y_bin  = rbinom(n, 1, plogis(eta - 0.6))
base$y_gam  = rgamma(n, shape = 2, rate = 2 / exp(eta))
base$y_gaus = rnorm(n, eta)

# canonical links
compare_fits(y_pois ~ x1 + x2 + f1 + f2, "poisson", base, "poisson (log)")
compare_fits(y_bin ~ x1 + x2 + f1 + f2, "logit", base, "binomial (logit)")
compare_fits(y_gaus ~ x1 + x2 + f1 + f2, gaussian(), base, "gaussian (identity)")

# non-canonical links: the tilt is the Fisher/IRLS weight
compare_fits(y_bin ~ x1 + f1 + f2, "probit", base, "binomial (probit)")
compare_fits(y_bin ~ x1 + f1, binomial(link = "cloglog"), base, "binomial (cloglog)")
compare_fits(y_gam ~ x1 + f1 + f2, Gamma(link = "log"), base, "Gamma (log)")
compare_fits(y_gam ~ x1 + f2, inverse.gaussian(link = "log"), base, "inv. gaussian (log)")

# quasi families (dispersion estimated)
compare_fits(y_pois ~ x1 + f1 + f2, quasipoisson(), base, "quasipoisson")
compare_fits(y_bin ~ x1 + f1, quasibinomial(), base, "quasibinomial")

# weights and offset
compare_fits(y_pois ~ x1 + f1 + f2, "poisson", base, "poisson, weights + offset",
             weights = ~w, offset = ~off)

# interactions: factor-factor cells and continuous-factor individual columns
compare_fits(y_pois ~ x1*f2 + f1, "poisson", base, "poisson, x1*f2 interaction")
compare_fits(y_pois ~ x1 + f1*f2, "poisson", base, "poisson, f1*f2 cells")

# purely categorical model: p = 0, everything runs on the cells
compare_fits(y_pois ~ f1*f2, "poisson", base, "poisson, cell-level only (p = 0)")

# fixest's i() syntax
compare_fits(y_pois ~ x1 + i(f1) + i(f2), "poisson", base, "poisson, i() dummies")

# collinearity: the tilted solve must reproduce the exclusion behavior
base$x1_dup = base$x1
res = compare_fits(y_pois ~ x1 + x1_dup + f1, "poisson", base, "poisson, collinear column")
stopifnot(identical(res$orig$collin.var, "x1_dup"),
          identical(res$tilt$collin.var, "x1_dup"))

# adversarial design: a column that is 0/1 on the screening sample but
# continuous elsewhere -- must be pushed back to the individual level by the
# cell cap, leaving the results exact
adv = base
adv$xs = rnorm(n)
adv$xs[seq(1, n, length.out = 5000)] = rbinom(5000, 1, 0.5)
compare_fits(y_pois ~ x1 + xs + f1, "poisson", adv, "poisson, screen-evading column")

# cluster-robust vcov computed on the tilted fit
orig = feglm(y_pois ~ x1 + f1 + f2, base, "poisson")
tilt = feglm(y_pois ~ x1 + f1 + f2, base, "poisson", tilt = TRUE)
stopifnot(max_diff(vcov(orig, cluster = ~f1), vcov(tilt, cluster = ~f1)) < tol_exact)
cat("  cluster-robust vcov identical\n")

# no cell structure at all: silent fallback to the standard algorithm
plain = feglm(y_gaus ~ x1 + x2, base, gaussian(), tilt = TRUE)
stopifnot(max_diff(coef(plain), coef(feglm(y_gaus ~ x1 + x2, base, gaussian()))) < tol_exact)
cat("  fallback without cell structure OK\n")

# fixed-effects absorption must be refused with a clear message
err = tryCatch(feglm(y_pois ~ x1 | f1, base, "poisson", tilt = TRUE),
               error = function(e) conditionMessage(e))
stopifnot(is.character(err), grepl("not compatible with fixed-effects", err))
cat("  fixed-effects guard OK\n")

cat("\nPart 1 passed: tilt reproduces the original algorithm for all families.\n\n")

#
# Part 2: what is gained ####
#

# The gain is in the estimation stage. For each design we fit the identical
# model twice and report: total fit time, time per IRLS iteration, and the
# discrepancy of the results (always at machine precision). The original
# algorithm's per-iteration weighted cross-product costs O(N (p+k)^2)
# (with sparsity tricks), plus a final dense O(N q^2) Hessian; the tilted
# algorithm's iteration costs O(N p^2) for the microdata pass plus
# O(J k^2) cell algebra, and the final Hessian is reused from the last
# iteration's cell assembly.

cat("PART 2: what is gained (identical results, cheaper iterations)\n\n")

bench_one = function(fml, family, data, label, reps = 3){
  # environments are built once so that the timing isolates the estimation
  # stage that the tilt reorganizes; model-matrix construction is identical
  # in both cases
  t_orig = t_tilt = Inf
  for(r in 1:reps){
    env_o = feglm(fml, data, family = family, only.env = TRUE)
    t_orig = min(t_orig, system.time(orig <- feglm.fit(env = env_o))[3])
    env_t = feglm(fml, data, family = family, tilt = TRUE, only.env = TRUE)
    t_tilt = min(t_tilt, system.time(tilt <- feglm.fit(env = env_t))[3])
  }

  it = orig$iterations
  stopifnot(max_diff(coef(orig), coef(tilt)) < 1e-6,
            max_diff(se(orig), se(tilt)) < 1e-6)

  cat(sprintf("  %-38s\n", label))
  cat(sprintf("    original: %6.2fs total, %6.3fs/iter   (%d iterations)\n",
              t_orig, t_orig / it, it))
  cat(sprintf("    tilt:     %6.2fs total, %6.3fs/iter   speedup x%.1f\n",
              t_tilt, t_tilt / it, t_orig / t_tilt))
  cat(sprintf("    max |coef diff| = %.1e, max |se diff| = %.1e\n\n",
              max_diff(coef(orig), coef(tilt)), max_diff(se(orig), se(tilt))))

  c(orig = t_orig, tilt = t_tilt)
}

set.seed(99)
nb = 4e5
db = data.frame(x1 = rnorm(nb), x2 = rnorm(nb),
                f1 = factor(sample(1:100, nb, TRUE)),
                f2 = factor(sample(1:20, nb, TRUE)),
                f3 = factor(sample(1:25, nb, TRUE)))
eta_b = 0.3*db$x1 - 0.2*db$x2 + 0.02*as.integer(db$f1) - 0.01*as.integer(db$f2)
db$y_pois = rpois(nb, exp(eta_b))
db$y_bin = rbinom(nb, 1, plogis(eta_b - 0.5))

cat(sprintf("N = %s observations\n\n", formatC(nb, big.mark = ",", format = "d")))

bench_one(y_pois ~ x1 + x2 + f1 + f2, "poisson", db,
          "Poisson, additive factors (k = 120, J = 2,000, p = 2)")

bench_one(y_bin ~ x1 + x2 + f1 + f2, "logit", db,
          "Logit, additive factors (k = 120, J = 2,000, p = 2)")

bench_one(y_pois ~ x1 + x2 + f3*f2, "poisson", db,
          "Poisson, interacted factors (k = 500, J = 500, p = 2)")

bench_one(y_pois ~ x1*f2 + f1, "poisson", db,
          "Poisson, continuous-factor interaction (k = 120, p = 21)")

cat("All tilt tests passed.\n")

#
# Part 3 (optional): heavyweight logit stress test ####
#

# Run with FIXEST_TILT_HEAVY=1. Entity-level binary attributes: N
# observations belong to E entities and carry k binary attributes constant
# within entity, with 65% prevalence -- so fixest's dummy-sparsity shortcut
# does not apply and the original algorithm runs its dense O(N q^2)
# weighted cross-product at every IRLS iteration (>10 minutes at these
# sizes), while the tilted algorithm detects the E cells and runs each
# iteration at O(N p^2) + O(E k^2) (under a minute), with identical output.

if(isTRUE(Sys.getenv("FIXEST_TILT_HEAVY") == "1")){

  cat("\nPART 3: heavyweight logit stress test\n\n")

  set.seed(2026)
  nh = 1e5; E = 2500; kh = 1700
  gh = sample.int(E, nh, TRUE)
  A = matrix(rbinom(E*kh, 1, 0.65), E, kh)
  Xh = cbind(x1 = rnorm(nh), x2 = rnorm(nh), A[gh, ])
  colnames(Xh) = c("x1", "x2", paste0("a", 1:kh))
  beta_h = c(0.4, -0.3, rnorm(kh, 0, 1.2/sqrt(kh)))
  yh = rbinom(nh, 1, plogis(drop(Xh %*% beta_h) - 1.2))
  rm(A); invisible(gc())

  cat(sprintf("  N = %s, E = %s cells, k = %d entity-level columns, q = %d\n",
              formatC(nh, big.mark = ","), formatC(E, big.mark = ","),
              kh, ncol(Xh) + 1))

  t_o = system.time(orig <- feglm.fit(yh, Xh, family = "logit"))
  cat(sprintf("  original: %7.1f s  (%d iterations, %.1f s/iter)\n",
              t_o[3], orig$iterations, t_o[3]/orig$iterations))

  t_t = system.time(tilt <- feglm.fit(yh, Xh, family = "logit", tilt = TRUE))
  cat(sprintf("  tilt:     %7.1f s  (%d iterations, %.1f s/iter)   speedup x%.1f\n",
              t_t[3], tilt$iterations, t_t[3]/tilt$iterations, t_o[3]/t_t[3]))
  cat(sprintf("  max |coef diff| = %.1e, max |se diff| = %.1e\n",
              max_diff(coef(orig), coef(tilt)), max_diff(se(orig), se(tilt))))

  stopifnot(max_diff(coef(orig), coef(tilt)) < 1e-8,
            max_diff(se(orig), se(tilt)) < 1e-8,
            orig$iterations == tilt$iterations)
  cat("\nHeavy stress test passed: identical fits.\n")
}
