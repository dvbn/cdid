# Test-only loading of immutable original algebra with ONE matched-start
# adaptation. This is not unadapted legacy byte parity or a production oracle.
ipw_reference_sha <- c(
  chained_ipw_link_original_6550dc.R = "6550dc592451cb7e9b7f45573d083fa14cacc1f3a09421a3330d31997183b0b5",
  weighted_pair_ipw_derivative_reference.R = "c96beff841db195f9958b9e60a8bd32ee8356c4f06b90e4ce0e5fc66bb2dcb45"
)

ipw_reference_match_start <- function(original, file) {
  # Exact original calls, including their distinct controls and algebra.
  expected <- switch(file,
    chained_ipw_link_original_6550dc.R = quote(stats::glm.fit(
      x = Z_use, y = D_use, weights = fit_weights,
      family = stats::binomial(link = "logit"), control = ctrl, intercept = FALSE)),
    weighted_pair_ipw_derivative_reference.R = quote(glm.fit(
      Z[use,,drop=FALSE], D[use], weights=case[use]/rr[use],
      family=binomial(), control=glm.control(epsilon=1e-12,maxit=100))),
    stop("unsupported matched-start reference: ", file, call. = FALSE)
  )
  matched_call <- expected
  design <- if (file == "chained_ipw_link_original_6550dc.R") expected$x else expected[[2L]]
  matched_call$start <- substitute(rep(0, ncol(design)), list(design = design))

  replace_one <- function(code, before, after) {
    count <- 0L
    walk <- function(node) {
      if (!is.call(node) && !is.expression(node)) return(node)
      if (is.call(node) && (identical(node[[1L]], quote(glm.fit)) ||
                            identical(node[[1L]], quote(stats::glm.fit)))) {
        count <<- count + 1L
        if (!identical(node, before)) stop("wrong glm.fit call in reference", call. = FALSE)
        return(after)
      }
      for (i in seq_along(node)) {
        # Single-bracket replacement retains literal NULL call arguments.
        if (!identical(node[[i]], quote(expr = ))) node[i] <- list(walk(node[[i]]))
      }
      node
    }
    result <- walk(code)
    if (count != 1L) stop("reference must contain exactly one glm.fit call", call. = FALSE)
    result
  }
  matched <- replace_one(original, expected, matched_call)
  reversed <- replace_one(matched, matched_call, expected)
  if (!identical(reversed, original)) stop("matched-start reversal differs from parsed original", call. = FALSE)
  matched
}

ipw_reference_source <- function(file, envir) {
  name <- basename(file)
  if (!name %in% names(ipw_reference_sha)) stop("unsupported matched-start reference: ", name, call. = FALSE)
  if (!file.exists(file)) stop("missing packaged reference: ", name, call. = FALSE)
  if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required to bind packaged references.", call. = FALSE)
  if (!identical(digest::digest(file = file, algo = "sha256"), ipw_reference_sha[[name]])) {
    stop("packaged reference bytes differ: ", name, call. = FALSE)
  }
  eval(ipw_reference_match_start(parse(file, keep.source = FALSE), name), envir = envir)
  invisible(envir)
}
