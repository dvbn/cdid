# Development simulation support for cdid (not package API).
#
# Seed contract. with_seed(seed, expr) evaluates `expr` after
#   set.seed(seed, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion",
#            sample.kind = "Rejection")
# and then restores the caller's RNG kinds and .Random.seed exactly (or
# removes .Random.seed if it did not exist). Simulation results therefore
# depend only on `seed`, and calling a helper never alters the caller's
# random-number stream. The MC harness reserves the next nonoverlapping
# L'Ecuyer stream for estimation/inference, separately from data generation.

with_seed <- function(seed, expr) {
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed != round(seed) || abs(seed) > .Machine$integer.max) {
    stop("`seed` must be a single whole number in integer range.",
         call. = FALSE)
  }
  genv <- globalenv()
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = genv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = genv, inherits = FALSE)
  on.exit({
    RNGkind(old_kind[1L], old_kind[2L], old_kind[3L])
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = genv)
    } else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
      rm(".Random.seed", envir = genv)
    }
  })
  set.seed(seed, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion",
           sample.kind = "Rejection")
  expr
}
