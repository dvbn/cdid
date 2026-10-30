# Forwarding to did without re-evaluating or re-creating arguments.
#
# `forward_call()` evaluates a call to `fun` in the wrapper frame `env`. That
# call passes each formal the caller supplied by symbol (`yname = yname`) and
# the wrapper's dots as `...`. Consequences:
#   * each argument expression is evaluated at most once, lazily, in its
#     original environment (wrappers, do.call() and dots behave as for a
#     direct did call);
#   * formals the caller did not supply are not passed, so did applies its
#     own defaults and missing() logic;
#   * `override` replaces the forwarded expression for supplied formals only.
# did records the delegated call (e.g. `did::att_gt(yname = yname, ...)`) in
# the objects it returns; nothing else in the result is touched.

arg_missing <- function(name, env) {
  eval(call("missing", as.name(name)), env)
}

forward_call <- function(fun, formal_names, env, override = list(),
                         dots = TRUE) {
  supplied <- formal_names[!vapply(formal_names, arg_missing, logical(1),
                                   env = env)]
  args <- lapply(supplied, as.name)
  names(args) <- supplied
  for (nm in intersect(names(override), supplied)) {
    args[[nm]] <- override[[nm]]
  }
  if (dots) {
    args <- c(args, list(as.name("...")))
  }
  eval(as.call(c(list(fun), args)), env)
}

cdid_abort <- function(message, class) {
  stop(structure(
    class = c(class, "cdid_error", "error", "condition"),
    list(message = message, call = NULL)
  ))
}
