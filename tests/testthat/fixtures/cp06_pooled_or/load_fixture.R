# Real-owner fixture loader. This never calls an OR fit/replay function and
# never replaces or injects bindings in the installed package namespace.
por_fixture_ns <- function(name) {
  ns <- asNamespace("cdid")
  stopifnot(isNamespace(ns), environmentIsLocked(ns), bindingIsLocked(name, ns))
  fn <- get(name, envir = ns, inherits = FALSE)
  stopifnot(is.function(fn), identical(environment(fn), ns))
  fn
}

por_fixture_csv <- function(directory, name) {
  read.csv(file.path(directory, name), stringsAsFactors = FALSE,
           check.names = FALSE, na.strings = c("NA", "NaN", "nan"))
}

load_pooled_or_reference <- function(directory, case = "NYT_F", order = NULL) {
  raw <- por_fixture_csv(directory, "raw.csv")
  raw$cohort <- as.numeric(raw$cohort)
  N <- nrow(raw)
  if (is.null(order)) order <- seq_len(N)
  stopifnot(is.numeric(order), length(order) == N,
            identical(sort(as.integer(order)), seq_len(N)))
  cases <- por_fixture_csv(directory, "cases.csv")
  meta <- cases[cases$case == case, , drop = FALSE]
  stopifnot(nrow(meta) == 1L)
  calendar <- por_fixture_csv(directory, "calendar.csv")
  pairs <- por_fixture_csv(directory, "pairs.csv")
  all_blocks <- por_fixture_csv(directory, "blocks.csv")
  blocks <- all_blocks[all_blocks$case == case, , drop = FALSE]
  all_links <- por_fixture_csv(directory, "links.csv")
  links <- all_links[all_links$case == case, , drop = FALSE]
  stopifnot(nrow(blocks) == meta$unique_blocks, nrow(links) == meta$links,
            !anyDuplicated(blocks$scientific_block_key), !anyDuplicated(links$link_id))
  control <- por_fixture_ns("ipw_owner_control")(
    max_chunk_rows = 23L, max_retained_entries = 1e7,
    max_matrix_entries = 2e7, max_work_ops = 1e10, max_output_entries = 1e7)
  Y <- as.matrix(raw[paste0("Y", seq_len(nrow(calendar)))])
  long <- data.frame(id = rep(raw$id, each = nrow(calendar)),
                     time = rep(calendar$time, N), y = as.vector(t(Y)),
                     g = rep(raw$cohort, each = nrow(calendar)))
  zero <- which(rowSums(!is.na(Y)) == 0L)
  # Keep one actual all-NA panel identity. Other zero-observer identities enter
  # only through the representative linked frame, using the real owner policy.
  if (length(zero) > 1L) long <- long[!long$id %in% raw$id[zero[-1L]], , drop = FALSE]
  panel <- por_fixture_ns("construct_ipw_panel_owner")(
    long, "id", "time", "y", "g", control)
  slope_columns <- c("x", "x2")[seq_len(meta$shared_slopes)]
  X <- cbind(x = raw$x, x2 = raw$x^2)[, seq_len(meta$shared_slopes), drop = FALSE]
  if (meta$intercept == 1L) X <- cbind(`(Intercept)` = 1, X)
  X[raw$A == 0L, ] <- NA_real_
  frame <- por_fixture_ns("construct_ipw_frame_owner")(
    panel, raw$id[order], X[order, , drop = FALSE], weights = raw$w[order],
    control = control, frame_g = raw$cohort[order],
    allow_frame_only = TRUE, design_on_A_only = TRUE)
  selection <- por_fixture_csv(directory, paste0("selection_", meta$reader_universe, ".csv"))
  selection <- selection[match(raw$id, selection$id), , drop = FALSE]
  Q <- as.matrix(selection[c("selection_x_neg", "selection_x_zero", "selection_x_pos")])
  U <- selection$U == 1L
  rho <- selection$rho
  reader <- por_fixture_ns("construct_cp04_frame_reader")(
    frame, raw$id[order], Q[order, , drop = FALSE], raw$cohort[order], raw$w[order],
    control, U = U[order], rho = rho[order], universe = meta$reader_universe,
    justification = "known transport under the frozen independent fixed atom law")
  views <- lapply(seq_len(nrow(links)), function(i) {
    link <- links[i, , drop = FALSE]
    block <- blocks[match(link$block_id, blocks$block_id), , drop = FALSE]
    K <- as.numeric(strsplit(block$controls, "|", fixed = TRUE)[[1L]])
    stopifnot(all(K > link$t), !link$target %in% K)
    por_fixture_ns("construct_ipw_link_view")(
      frame, link$s, link$t, g = link$target, K = K, control = control,
      link_id = link$link_id, functional_key = link$functional_key)
  })
  names(views) <- links$link_id
  rows <- por_fixture_csv(directory, paste0(case, "_rows.csv"))
  link_rows <- por_fixture_csv(directory, paste0(case, "_link_rows.csv"))
  all_points <- por_fixture_csv(directory, "points.csv")
  points <- all_points[all_points$case == case, , drop = FALSE]
  all_coordinates <- por_fixture_csv(directory, "coordinates.csv")
  coordinates <- all_coordinates[all_coordinates$case == case, , drop = FALSE]
  ids <- reader$state$id
  rows <- rows[match(ids, rows$id), , drop = FALSE]
  link_rows <- link_rows[match(ids, link_rows$id), , drop = FALSE]
  stopifnot(identical(as.numeric(rows$id), as.numeric(ids)),
            identical(as.numeric(link_rows$id), as.numeric(ids)))
  point_group <- function(group) {
    x <- points[points$group == group, , drop = FALSE]
    setNames(x$value[order(x$group_index)], x$name[order(x$group_index)])
  }
  list(panel = panel, frame = frame, reader = reader, views = views,
       raw = raw[match(ids, raw$id), , drop = FALSE], control = control,
       Q = Q[match(ids, raw$id), , drop = FALSE], rho = rho[match(ids, raw$id)],
       U = U[match(ids, raw$id)], X = X[match(ids, raw$id), , drop = FALSE],
       catalog = list(case = meta, calendar = calendar, pairs = pairs,
                      blocks = blocks, links = links),
       outcome_recipe = list(intercept = meta$intercept == 1L,
                             shared_columns = slope_columns,
                             pair_transitions = setNames(pairs$modeled_transitions, pairs$pair_id)),
       expectations = list(points = points, coordinates = coordinates,
                           rows = as.matrix(rows[, -1L, drop = FALSE]),
                           phi = as.matrix(link_rows[, -1L, drop = FALSE]),
                           upstream = point_group("upstream"), gamma = point_group("gamma"),
                           mu = point_group("mu"), tau = point_group("tau"),
                           eta = point_group("eta"), lambda = point_group("lambda")))
}
