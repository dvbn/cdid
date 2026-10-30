# Formals of did 2.5.1 aggte(), in order.
upstream_aggte_formals <- c(
  "MP", "type", "balance_e", "min_e", "max_e", "na.rm", "bstrap", "biters",
  "cband", "alp", "clustervars"
)

aggte <- function(MP,
                  type = "group",
                  balance_e = NULL,
                  min_e = -Inf,
                  max_e = Inf,
                  na.rm = FALSE,
                  bstrap = NULL,
                  biters = NULL,
                  cband = NULL,
                  alp = NULL,
                  clustervars = NULL) {
  # A chained result (cdid_MP) is aggregated from its retained state without
  # refitting, with did's selector, inheritance and flag conventions
  # (R/chained_did_request.R). did::aggte() forces `MP` first in any case, so
  # inspecting its class here changes nothing for ordinary did results, which
  # are still aggregated by did itself, unchanged and lazily forwarded.
  if (!arg_missing("MP", environment()) && inherits(MP, "cdid_MP")) {
    return(chained_aggte_request(MP, type = type, balance_e = balance_e, min_e = min_e, max_e = max_e,
                                 na.rm = na.rm, bstrap = bstrap, biters = biters, cband = cband, alp = alp,
                                 clustervars = clustervars, call = sys.call()))
  }
  forward_call(quote(did::aggte), upstream_aggte_formals, environment(),
               dots = FALSE)
}
