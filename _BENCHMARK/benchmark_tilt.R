#----------------------------------------------------------------------------#
# Parametric benchmark of feglm's exponential-tilt algorithm                 #
#                                                                            #
# Sweeps each parameter that drives computation time along a grid, holding  #
# the others at a base configuration, and records total / per-iteration     #
# times of the original and tilted algorithms plus the (machine-precision)  #
# discrepancy of the results:                                               #
#                                                                            #
#   N      : number of observations                                          #
#   k      : number of cell-level columns (entity attributes or dummies)     #
#   p      : number of individual-level continuous columns                   #
#   J (=E) : number of cells (entities)                                      #
#   prev   : prevalence of the binary attribute columns -- crosses the       #
#            threshold of fixest's dummy-sparsity shortcut (rows > 50%       #
#            zeros), i.e. the regime switch between the sparse and dense     #
#            paths of the original algorithm                                 #
#   family : gaussian / poisson / logit / probit / Gamma-log                 #
#   design : "attributes" (dense entity-level binaries, prevalence 0.65)     #
#            vs "dummies" (classic one-hot factor dummies, sparse)           #
#                                                                            #
# Base configuration: attributes design, logit, N = 1e5, k = 200, E = 1e3,  #
# p = 2, prevalence 0.65. Results are written to                            #
# tilt_benchmark_results.csv and plotted in tilt_benchmark.png.             #
#----------------------------------------------------------------------------#

library(fixest)
setFixest_notes(FALSE)

out_csv = "tilt_benchmark_results.csv"
out_png = "tilt_benchmark.png"

gen_data = function(N, E, k, p, prevalence, design, family, seed){
  set.seed(seed)
  g = sample.int(E, N, TRUE)
  if(design == "attributes"){
    A = matrix(rbinom(E * k, 1, prevalence), E, k)
    C = A[g, , drop = FALSE]
  } else {
    # one-hot dummies of the entity factor itself: k = E - 1 sparse columns
    C = matrix(0, N, k)
    sel = g > 1
    C[cbind(which(sel), g[sel] - 1L)] = 1
  }
  colnames(C) = paste0("a", seq_len(k))
  if(p > 0){
    Xp = matrix(rnorm(N * p), N, p, dimnames = list(NULL, paste0("x", seq_len(p))))
    X = cbind(Xp, C)
    beta = c(rnorm(p, 0, 0.3), rnorm(k, 0, 1.2 / sqrt(k)))
  } else {
    X = C
    beta = rnorm(k, 0, 1.2 / sqrt(k))
  }
  eta = drop(X %*% beta)
  y = switch(family,
    gaussian = eta + rnorm(N),
    poisson  = rpois(N, exp(0.5 * eta - 0.5)),
    logit    = rbinom(N, 1, plogis(eta - 1.2)),
    probit   = rbinom(N, 1, pnorm(0.7 * eta - 1)),
    gamma    = rgamma(N, shape = 2, rate = 2 / exp(0.5 * eta)))
  list(y = y, X = X)
}

fam_obj = function(family){
  switch(family,
    gaussian = gaussian(), poisson = "poisson", logit = "logit",
    probit = "probit", gamma = Gamma(link = "log"))
}

time_fit = function(y, X, family, tilt){
  t = system.time(fit <- feglm.fit(y, X, family = fam_obj(family), tilt = tilt))[3]
  if(t < 20){
    t = min(t, system.time(fit <- feglm.fit(y, X, family = fam_obj(family), tilt = tilt))[3])
  }
  list(t = as.numeric(t), fit = fit)
}

results = list()
run_point = function(sweep, design = "attributes", family = "logit",
                     N = 1e5, E = 1000, k = 200, p = 2, prevalence = 0.65){
  if(design == "dummies") k = E - 1
  seed = 31000 + length(results)
  d = gen_data(N, E, k, p, prevalence, design, family, seed)

  o = time_fit(d$y, d$X, family, tilt = FALSE)
  b = time_fit(d$y, d$X, family, tilt = TRUE)

  dcoef = max(abs(coef(o$fit) - coef(b$fit)))
  dse = max(abs(se(o$fit) - se(b$fit)))
  stopifnot(dcoef < 1e-6, dse < 1e-4)

  row = data.frame(sweep = sweep, design = design, family = family,
                   N = N, E = E, k = k, p = p, q = ncol(d$X) + 1,
                   prevalence = prevalence,
                   iters_orig = o$fit$iterations, iters_tilt = b$fit$iterations,
                   t_orig = o$t, t_tilt = b$t,
                   per_iter_orig = o$t / o$fit$iterations,
                   per_iter_tilt = b$t / b$fit$iterations,
                   speedup = o$t / b$t, dcoef = dcoef, dse = dse)
  results[[length(results) + 1]] <<- row
  write.csv(do.call(rbind, results), out_csv, row.names = FALSE)
  cat(sprintf("[%s] %s %s N=%.0e E=%d k=%d p=%d prev=%.2f | orig %6.1fs (%d it) | tilt %6.1fs (%d it) | x%.1f | dcoef %.0e\n",
      format(Sys.time(), "%H:%M:%S"), sweep, family, N, E, k, p, prevalence,
      o$t, o$fit$iterations, b$t, b$fit$iterations, o$t / b$t, dcoef))
  invisible(gc())
}

cat("Base configuration: attributes design, logit, N=1e5, E=1000, k=200, p=2, prev=0.65\n\n")

# 1) N sweep
for(N in c(2.5e4, 5e4, 1e5, 2e5, 4e5, 8e5)) run_point("N", N = N)

# 2) k sweep (E scaled to stay above k)
for(k in c(25, 50, 100, 200, 400, 800)) run_point("k", k = k, E = max(1000, 2 * k))

# 3) p sweep (including the pure cell-level model p = 0)
for(p in c(0, 1, 2, 4, 8, 16, 32, 64)) run_point("p", p = p)

# 4) J sweep (number of cells; E stays above k to keep the design full rank)
for(E in c(500, 1000, 2000, 8000, 32000)) run_point("J", E = E)

# 5) prevalence sweep (crosses the sparse-path threshold of the original)
for(pr in c(0.03, 0.06, 0.12, 0.25, 0.35, 0.5, 0.65)) run_point("prevalence", prevalence = pr)

# 6) family sweep
for(f in c("gaussian", "poisson", "logit", "probit", "gamma")) run_point("family", family = f)

# 7) classic sparse factor dummies (k = E - 1 one-hot columns)
for(E in c(26, 51, 101, 201, 401, 801)) run_point("dummies", design = "dummies", E = E)

res = do.call(rbind, results)
write.csv(res, out_csv, row.names = FALSE)
cat("\nResults written to", out_csv, "\n")

#
# Plot: the full picture ####
#

png(out_png, width = 1500, height = 950, res = 105)
par(mfrow = c(2, 3), mar = c(4.2, 4.2, 2.5, 1), mgp = c(2.6, 0.8, 0))

panel = function(sub, xvar, xlab, log = "xy", xvals = NULL){
  s = res[res$sweep == sub, ]
  x = if(is.null(xvals)) s[[xvar]] else xvals
  ylim = range(c(s$t_orig, s$t_tilt))
  plot(x, s$t_orig, type = "b", pch = 19, col = "firebrick", log = log,
       xlab = xlab, ylab = "fit time (s)", ylim = ylim,
       main = paste0("sweep: ", sub))
  lines(x, s$t_tilt, type = "b", pch = 17, col = "royalblue")
  text(x, sqrt(s$t_orig * s$t_tilt), labels = sprintf("x%.1f", s$speedup),
       cex = 0.75, col = "gray30")
  legend("topleft", c("original", "tilt"), col = c("firebrick", "royalblue"),
         pch = c(19, 17), lty = 1, bty = "n", cex = 0.9)
}

panel("N", "N", "N (observations)")
panel("k", "k", "k (cell-level columns)")
panel("p", "p", "p (individual-level columns)", log = "y")
panel("J", "E", "J (number of cells)")
panel("prevalence", "prevalence", "attribute prevalence", log = "y")

s = res[res$sweep == "family", ]
bp = barplot(rbind(s$t_orig, s$t_tilt), beside = TRUE, names.arg = s$family,
             col = c("firebrick", "royalblue"), ylab = "fit time (s)",
             main = "sweep: family", legend.text = c("original", "tilt"),
             args.legend = list(bty = "n"))
text(colMeans(bp), pmax(s$t_orig, s$t_tilt) * 0.5,
     labels = sprintf("x%.1f", s$speedup), cex = 0.8)

dev.off()
cat("Plot written to", out_png, "\n")

#
# Console summary ####
#

cat("\n===== SUMMARY =====\n")
for(sw in unique(res$sweep)){
  s = res[res$sweep == sw, ]
  cat(sprintf("\n%s sweep: speedup %.1f - %.1f (median %.1f)\n",
      sw, min(s$speedup), max(s$speedup), median(s$speedup)))
}
cat(sprintf("\nmax coefficient discrepancy over all %d runs: %.1e\n",
    nrow(res), max(res$dcoef)))
