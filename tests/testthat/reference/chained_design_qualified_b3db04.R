# Private CP02 augmented incidence algebra. Functional keys are caller-declared
# population-operator identities, NEVER certified from fitted rates or metadata.
# Default storage is O(vertices + links): two endpoints and one parent forest.

fmt_design_num <- function(x) sprintf("%.17g", x) # Display only, never identity.

design_vector <- function(x, n, type = "numeric", finite = TRUE) {
  is.atomic(x) && is.null(dim(x)) && length(x) == n && !is.object(x) &&
    (if (type == "numeric") is.numeric(x) else if (type == "logical") is.logical(x) else is.character(x)) &&
    !anyNA(x) && (type != "numeric" || !finite || all(is.finite(x)))
}

design_columns <- function(data, required, optional = character()) {
  if (!is.data.frame(data) || any(vapply(required, function(x) sum(names(data) == x, na.rm = TRUE) != 1L, logical(1L))) ||
      any(vapply(intersect(optional, names(data)), function(x) sum(names(data) == x, na.rm = TRUE) != 1L, logical(1L)))) {
    cdid_abort("Missing or ambiguous selected design column.", "cdid_design_input_error")
  }
}

design_bound <- function(x) {
  if (!design_vector(x, 1L) || x < 0 || x != floor(x) || x > 2^53 - 1) {
    cdid_abort("Supply an explicit finite integer entry bound.", "cdid_design_bound_error")
  }
  as.double(x)
}

# Truthful metadata validation only. No key generation, equality inference,
# inferred control group, normalized-rate fingerprint, or boolean certification.
validate_functional_metadata <- function(links, check_consistency = TRUE) {
  n <- nrow(links)
  qual_cols <- c("s_qualified", "t_qualified", "s_coordinate_role", "t_coordinate_role",
                 "max_qualified_horizon", "uses_post_nuisance")
  has_qual_any <- any(qual_cols %in% names(links))
  # Legacy inputs retain unknown columns as data; select pool metadata only
  # when explicit qualification metadata selects the new branch.
  optional <- c("control_group", "horizon", "profile_id", "composition_corrected",
                "post_shareable", if (has_qual_any) "control_pool", qual_cols)
  design_columns(links, "functional_key", optional)
  fk <- links$functional_key
  if (!design_vector(fk, n, "character") || any(!nzchar(fk))) {
    cdid_abort("Declare nonempty opaque population-functional keys explicitly.", "cdid_design_functional_error")
  }
  if (has_qual_any) {
    if (!all(c("s_qualified", "t_qualified") %in% names(links)) ||
        !design_vector(links$s_qualified, n, "logical") ||
        !design_vector(links$t_qualified, n, "logical")) {
      cdid_abort("Qualified endpoint declarations require non-missing scalar logical `s_qualified` and `t_qualified` columns.",
        c("cdid_design_qualification_error", "cdid_design_post_share_error", "cdid_design_functional_error", "cdid_design_input_error"))
    }
    if (any(c("s_coordinate_role", "t_coordinate_role") %in% names(links))) {
      if (!all(c("s_coordinate_role", "t_coordinate_role") %in% names(links)) ||
          !design_vector(links$s_coordinate_role, n, "character") || any(!nzchar(links$s_coordinate_role)) ||
          !design_vector(links$t_coordinate_role, n, "character") || any(!nzchar(links$t_coordinate_role))) {
        cdid_abort("Coordinate role declarations require non-empty scalar character `s_coordinate_role` and `t_coordinate_role` columns.",
          c("cdid_design_qualification_error", "cdid_design_functional_error", "cdid_design_input_error"))
      }
    }
    if ("max_qualified_horizon" %in% names(links)) {
      mh <- links$max_qualified_horizon
      if (!is.atomic(mh) || !is.null(dim(mh)) || is.object(mh) || length(mh) != n ||
          !(is.numeric(mh) || (is.logical(mh) && all(is.na(mh)))) ||
          any(is.nan(mh) | mh == -Inf, na.rm = TRUE)) {
        cdid_abort("Invalid `max_qualified_horizon` metadata.",
          c("cdid_design_horizon_error", "cdid_design_qualification_error", "cdid_design_functional_error", "cdid_design_input_error"))
      }
    }
    if ("uses_post_nuisance" %in% names(links) && !design_vector(links$uses_post_nuisance, n, "logical")) {
      cdid_abort("Invalid `uses_post_nuisance` logical metadata.",
        c("cdid_design_qualification_error", "cdid_design_functional_error", "cdid_design_input_error"))
    }
  }
  if (has_qual_any && "control_pool" %in% names(links)) {
    cp <- links$control_pool
    if (!design_vector(cp, n, "character") || any(!nzchar(cp) | cp == "{}")) {
      cdid_abort("Invalid `control_pool` metadata.",
        c("cdid_design_qualification_error", "cdid_design_functional_error", "cdid_design_input_error"))
    }
  }
  for (nm in intersect(c("control_group", "profile_id"), names(links))) {
    v <- links[[nm]]
    if (!design_vector(v, n, "character") || any(!nzchar(v))) cdid_abort("Invalid functional metadata.", "cdid_design_functional_error")
  }
  if ("control_group" %in% names(links) && any(!links$control_group %in% c("nevertreated", "notyettreated"))) {
    cdid_abort("Invalid control-group metadata.", "cdid_design_functional_error")
  }
  for (nm in intersect(c("composition_corrected", "post_shareable"), names(links))) {
    if (!design_vector(links[[nm]], n, "logical")) cdid_abort("Invalid logical metadata.", "cdid_design_functional_error")
  }
  if ("horizon" %in% names(links) &&
      (!design_vector(links$horizon, n, finite = FALSE) || any(links$horizon == -Inf))) {
    cdid_abort("Invalid horizon metadata.", "cdid_design_horizon_error")
  }
  rates <- attr(links, "profile_rates", exact = TRUE)
  if (!is.null(rates)) {
    if (is.matrix(rates) && is.numeric(rates) && nrow(rates) == n && ncol(rates) > 0L) {
      good <- all(is.finite(rates) & rates > 0)
    } else if (is.list(rates) && length(rates) == n) {
      good <- all(vapply(rates, function(v) design_vector(v, length(v)) && length(v) > 0L && all(v > 0), logical(1L)))
    } else good <- FALSE
    if (!good) cdid_abort("Invalid rate metadata; rates do not certify operators.", "cdid_design_functional_error")
  }
  # Reject contradictory descriptive metadata; agreement does not prove equality.
  if (check_consistency && "g" %in% names(links) && n > 1L) {
    ord <- order(links$g, fk, method = "radix")
    a <- ord[-n]; b <- ord[-1L]
    same <- links$g[a] == links$g[b] & fk[a] == fk[b]
    for (nm in intersect(c("control_group", "horizon", "profile_id", "composition_corrected",
                           if (has_qual_any) "control_pool"), names(links))) {
      if (any(same & links[[nm]][a] != links[[nm]][b])) {
        cdid_abort("A declared key has contradictory descriptive metadata.",
          c("cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error"))
      }
    }
    if ("max_qualified_horizon" %in% names(links)) {
      mh <- as.double(links$max_qualified_horizon)
      same_mh <- (is.na(mh[a]) & is.na(mh[b])) | (!is.na(mh[a]) & !is.na(mh[b]) & mh[a] == mh[b])
      if (any(same & !same_mh)) {
        cdid_abort("A declared key has contradictory `max_qualified_horizon` metadata.",
          c("cdid_design_horizon_error", "cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error"))
      }
    }
    # Optional per-key rate descriptions must be internally consistent. This
    # is a conservative metadata contract, not a test of operator equality or
    # inequality. Agreement creates no key; a correction flag cannot bypass
    # conflicting metadata. Raw nuisance rates that are not an operator's
    # profile should not be supplied as its per-key profile description.
    if (!is.null(rates)) {
      for (j in which(same)) {
        ra <- if (is.matrix(rates)) rates[a[j], ] else rates[[a[j]]]
        rb <- if (is.matrix(rates)) rates[b[j], ] else rates[[b[j]]]
        if (!identical(as.double(ra), as.double(rb))) cdid_abort(
          "A declared key has inconsistent rate-profile metadata; no operator equality is certified.",
          c("cdid_design_profile_merge_error", "cdid_design_functional_error"))
      }
    }
  }
  invisible(TRUE)
}

validate_design_links <- function(links) {
  if (!is.data.frame(links)) cdid_abort("`links` must be a data.frame.", "cdid_design_input_error")
  if (!("ref" %in% names(links)) && "reference" %in% names(links)) {
    design_columns(links, "reference"); links$ref <- links$reference
  }
  if (!("functional_key" %in% names(links))) cdid_abort(
    "Explicit functional keys are required; rates, labels and booleans cannot certify them.",
    c("cdid_design_profile_merge_error", "cdid_design_functional_error"))
  design_columns(links, c("g", "s", "t", "ref", "functional_key"))
  n <- nrow(links)
  if (n == 0L) cdid_abort("At least one declared link is required.", "cdid_design_input_error")
  for (nm in c("g", "s", "t", "ref")) {
    if (!design_vector(links[[nm]], n)) cdid_abort("Design dates/cohorts must be finite scalar numeric columns.", "cdid_design_input_error")
  }
  validate_functional_metadata(links, check_consistency = FALSE)
  if (any(links$s >= links$t)) cdid_abort("Each link requires s < t.", "cdid_design_endpoint_error")
  ord <- order(links$g, method = "radix")
  if (any(links$ref >= links$g) || (n > 1L && any(
      links$g[ord[-1L]] == links$g[ord[-n]] & links$ref[ord[-1L]] != links$ref[ord[-n]]))) {
    cdid_abort("Each cohort requires a unique reference preceding treatment.", "cdid_design_reference_error")
  }
  has_qual <- "s_qualified" %in% names(links)
  post <- links$t >= links$g | links$s >= links$g
  pre <- (links$s < links$g & links$s != links$ref) | (links$t < links$g & links$t != links$ref)
  if (!has_qual) {
    if ("post_shareable" %in% names(links) && any(post & !links$post_shareable)) cdid_abort(
      "Post sharing requires the declared full CS5/no-anticipation/overlap population qualification.", "cdid_design_post_share_error")
    if ("horizon" %in% names(links) && any(links$t >= links$g & links$horizon < pmax(links$g, links$t))) {
      cdid_abort("Post control horizon must cover max(g,t).", "cdid_design_horizon_error")
    }
    if (any(pre & grepl("^(nyt:)?(h|horizon)=[0-9.+eE-]+$", links$functional_key))) cdid_abort(
      "A bare horizon label does not declare a complete population comparison operator.",
      c("cdid_design_profile_merge_error", "cdid_design_functional_error"))
  } else {
    sq <- links$s_qualified; tq <- links$t_qualified
    if (any((links$s == links$ref & !sq) | (links$t == links$ref & !tq))) {
      cdid_abort("Reference endpoints must be declared qualified.",
        c("cdid_design_qualification_error", "cdid_design_reference_error", "cdid_design_functional_error", "cdid_design_input_error"))
    }
    if (any((links$s < links$g & links$s != links$ref & sq) | (links$t < links$g & links$t != links$ref & tq))) {
      cdid_abort("Early non-reference coordinates cannot be shared or qualified under PT-post.",
        c("cdid_design_qualification_error", "cdid_design_post_share_error", "cdid_design_functional_error", "cdid_design_input_error"))
    }
    if ("s_coordinate_role" %in% names(links)) {
      check_role <- function(u, g_v, ref_v, q_v, r_v) {
        ok_ref <- u == ref_v & r_v == "reference"
        ok_pre <- u < g_v & u != ref_v & !q_v & r_v %in% c("pre_nuisance", "pre")
        ok_att <- u >= g_v & q_v & r_v %in% c("post_att", "post")
        ok_pnu <- u >= g_v & !q_v & r_v == "post_nuisance"
        all(ok_ref | ok_pre | ok_att | ok_pnu)
      }
      if (!check_role(links$s, links$g, links$ref, sq, links$s_coordinate_role) ||
          !check_role(links$t, links$g, links$ref, tq, links$t_coordinate_role)) {
        cdid_abort("Endpoint coordinate role contradicts its declared qualification or date.",
          c("cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error", "cdid_design_input_error"))
      }
    }
    if ("uses_post_nuisance" %in% names(links)) {
      expected_upn <- (links$s >= links$g & !sq) | (links$t >= links$g & !tq)
      if (any(links$uses_post_nuisance != expected_upn)) {
        cdid_abort("`uses_post_nuisance` contradicts endpoint qualification.",
          c("cdid_design_qualification_error", "cdid_design_functional_error", "cdid_design_input_error"))
      }
    }
    if ("post_shareable" %in% names(links)) {
      expected_ps <- (links$s < links$g | sq) & (links$t < links$g | tq)
      if (any(links$post_shareable != expected_ps)) {
        cdid_abort("`post_shareable` must be TRUE iff every post endpoint on the link is qualified.",
          c("cdid_design_post_share_error", "cdid_design_qualification_error", "cdid_design_functional_error", "cdid_design_input_error"))
      }
    }
    if ("horizon" %in% names(links) && any(links$t >= links$g & links$horizon < pmax(links$g, links$t))) {
      cdid_abort("Post control horizon must cover max(g,t).", "cdid_design_horizon_error")
    }
    if ("max_qualified_horizon" %in% names(links)) {
      mh <- as.double(links$max_qualified_horizon)
      if (any(!is.na(mh) & mh < links$g) ||
          any(links$s >= links$g & sq & (is.na(mh) | mh < links$s)) ||
          any(links$t >= links$g & tq & (is.na(mh) | mh < links$t)) ||
          any(links$s >= links$g & !sq & !is.na(mh) & mh >= links$s) ||
          any(links$t >= links$g & !tq & !is.na(mh) & mh >= links$t)) {
        cdid_abort("`max_qualified_horizon` contradicts endpoint qualification or cohort treatment date.",
          c("cdid_design_horizon_error", "cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error", "cdid_design_input_error"))
      }
    }
    nuis <- pre | (links$s >= links$g & !sq) | (links$t >= links$g & !tq)
    if (any(nuis & grepl("^(nyt:)?(h|horizon)=[0-9.+eE-]+$", links$functional_key))) cdid_abort(
      "A bare horizon label does not declare a complete population comparison operator.",
      c("cdid_design_profile_merge_error", "cdid_design_functional_error"))
  }
  if (has_qual && "control_pool" %in% names(links)) {
    cp <- links$control_pool
    mh_vec <- if ("max_qualified_horizon" %in% names(links)) as.double(links$max_qualified_horizon) else NULL
    for (i in seq_len(n)) {
      if (grepl("^\\{.*\\}$", cp[i])) {
        inner <- substring(cp[i], 2L, nchar(cp[i]) - 1L)
        parts <- strsplit(inner, ",", fixed = TRUE)[[1L]]
        vals <- suppressWarnings(as.double(parts))
        if (!length(parts) || grepl(",[[:space:]]*$", inner) || any(!nzchar(parts)) || anyNA(vals) || any(is.nan(vals) | vals == -Inf) ||
            anyDuplicated(vals) || any(vals <= links$g[i])) {
          cdid_abort("`control_pool` must contain unique control cohorts strictly after target g.",
            c("cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error", "cdid_design_input_error"))
        }
        if ("control_group" %in% names(links) &&
            all(is.infinite(vals)) != (links$control_group[i] == "nevertreated")) {
          cdid_abort("`control_pool` contradicts `control_group` metadata.",
            c("cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error"))
        }
        if (!is.null(mh_vec) && !is.na(mh_vec[i]) && any(is.finite(vals) & vals <= mh_vec[i])) {
          cdid_abort("`control_pool` contradicts `max_qualified_horizon`.",
            c("cdid_design_horizon_error", "cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error"))
        }
        if (has_qual) {
          if ((links$s[i] >= links$g[i] && links$s_qualified[i] != all(vals > links$s[i])) ||
              (links$t[i] >= links$g[i] && links$t_qualified[i] != all(vals > links$t[i]))) {
            cdid_abort("Endpoint qualification contradicts declared `control_pool`.",
              c("cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error", "cdid_design_input_error"))
          }
        }
      }
    }
  }
  validate_functional_metadata(links)
  # Injective integer tuple identities; original doubles are matched exactly.
  gi <- match(links$g, sort(unique(links$g)))
  dates <- sort(unique(c(links$s, links$t)))
  fi <- match(links$functional_key, sort(unique(links$functional_key), method = "radix"))
  if (has_qual) {
    ep_gi <- c(gi, gi); ep_fi <- c(fi, fi)
    ep_di <- match(c(links$s, links$t), dates)
    ep_qual <- c(links$s_qualified, links$t_qualified)
    ord_ep <- order(ep_gi, ep_fi, ep_di, method = "radix")
    a_ep <- ord_ep[-length(ord_ep)]; b_ep <- ord_ep[-1L]
    same_ep <- ep_gi[a_ep] == ep_gi[b_ep] & ep_fi[a_ep] == ep_fi[b_ep] & ep_di[a_ep] == ep_di[b_ep]
    if (any(same_ep & ep_qual[a_ep] != ep_qual[b_ep])) {
      cdid_abort("Contradictory endpoint qualification for repeated (g, functional_key, date).",
        c("cdid_design_qualification_error", "cdid_design_post_share_error", "cdid_design_profile_merge_error", "cdid_design_functional_error", "cdid_design_input_error"))
    }
    if ("s_coordinate_role" %in% names(links)) {
      ep_role <- c(links$s_coordinate_role, links$t_coordinate_role)
      if (any(same_ep & ep_role[a_ep] != ep_role[b_ep])) {
        cdid_abort("Contradictory endpoint coordinate role for repeated (g, functional_key, date).",
          c("cdid_design_qualification_error", "cdid_design_profile_merge_error", "cdid_design_functional_error", "cdid_design_input_error"))
      }
    }
  }
  keys <- paste(gi, match(links$s, dates), match(links$t, dates), fi, sep = ":")
  if (anyDuplicated(keys)) cdid_abort("Duplicate declared links.", "cdid_design_duplicate_link_error")
  links
}

normalize_target_post_times <- function(target_post_times, cohorts, refs) {
  empty <- data.frame(g = numeric(), time = numeric())
  if (is.null(target_post_times)) return(empty)
  if (design_vector(target_post_times, length(target_post_times))) {
    grid <- expand.grid(time = sort(unique(target_post_times)), g = cohorts,
      KEEP.OUT.ATTRS = FALSE)[, c("g", "time")]
  } else if (is.data.frame(target_post_times)) {
    tc <- if ("time" %in% names(target_post_times)) "time" else "t"
    design_columns(target_post_times, c("g", tc))
    gv <- target_post_times$g; tv <- target_post_times[[tc]]
    if (!design_vector(gv, nrow(target_post_times)) || !design_vector(tv, nrow(target_post_times))) {
      cdid_abort("Targets require finite scalar numeric columns.", "cdid_design_target_error")
    }
    grid <- data.frame(g = as.double(gv), time = as.double(tv))
  } else cdid_abort("Invalid requested post targets.", "cdid_design_target_error")
  ci <- match(grid$g, cohorts)
  if (anyNA(ci) || any(grid$time < grid$g | grid$time <= refs[ci])) cdid_abort("Invalid cohort or post target date.", "cdid_design_target_error")
  # Exact tuple deduplication, avoiding data-frame numeric string conversion.
  key <- paste(ci, match(grid$time, sort(unique(grid$time))), sep = ":")
  grid <- grid[!duplicated(key), , drop = FALSE]
  grid <- grid[order(grid$g, grid$time, method = "radix"), , drop = FALSE]
  rownames(grid) <- NULL
  grid
}

# Temporary CSR adjacency is linear size and discarded after BFS. Sorted
# neighbours and canonical link ranks give deterministic parents under row order.
incidence_forest <- function(nv, from, to, anchors, edge_order = seq_along(from)) {
  m <- length(from)
  if (2 * as.double(m) > .Machine$integer.max) cdid_abort("Graph exceeds supported adjacency index range.", "cdid_design_bound_error")
  u <- c(from, to); v <- c(to, from)
  e <- rep(seq_len(m), 2L); sign <- c(rep(1L, m), rep(-1L, m))
  ord <- order(u, v, rep(edge_order, 2L), method = "radix")
  degree <- tabulate(u, nbins = nv)
  v <- v[ord]; e <- e[ord]; sign <- sign[ord]
  end <- cumsum(as.double(degree)); start <- end - degree + 1
  comp <- integer(nv); parent <- parent_edge <- integer(nv); parent_sign <- integer(nv)
  queue <- integer(nv); roots <- integer(nv); nc <- 0L
  for (root in c(anchors, seq_len(nv))) {
    if (comp[root] != 0L) next
    nc <- nc + 1L; roots[nc] <- root; comp[root] <- nc
    head <- tail <- 1L; queue[1L] <- root
    while (head <= tail) {
      a <- queue[head]; head <- head + 1L
      if (degree[a] == 0L) next
      for (k in seq.int(start[a], end[a])) {
        b <- v[k]
        if (comp[b] == 0L) {
          comp[b] <- nc; parent[b] <- a; parent_edge[b] <- e[k]; parent_sign[b] <- sign[k]
          tail <- tail + 1L; queue[tail] <- b
        }
      }
    }
  }
  list(component_id = comp, n_components = nc, roots = roots[seq_len(nc)],
       parent_vertex = parent, parent_link = parent_edge, parent_sign = parent_sign, degree = degree)
}

build_chained_design <- function(links, target_post_times = NULL) {
  links <- validate_design_links(links)
  has_qual <- "s_qualified" %in% names(links)
  g <- as.double(links$g); s <- as.double(links$s); t <- as.double(links$t)
  ref <- as.double(links$ref); fk <- links$functional_key; m <- length(g)
  cohorts <- sort(unique(g)); gi <- match(g, cohorts); refs <- ref[match(cohorts, g)]
  targets <- normalize_target_post_times(target_post_times, cohorts, refs)
  fn <- sort(unique(fk), method = "radix"); fi <- match(fk, fn)
  dates <- sort(unique(c(refs, s, t, targets$time)))
  sq <- if (has_qual) links$s_qualified else (s == ref | s >= g)
  tq <- if (has_qual) links$t_qualified else (t == ref | t >= g)
  vertex_key <- function(ci, time, f, qual) {
    # -1 is reference; 0 is shared qualified post; positive indices are declared pre/post nuisance keys.
    kind <- ifelse(time == refs[ci], -1L, ifelse(qual, 0L, f))
    paste(ci, match(time, dates), kind, sep = ":")
  }
  sk <- vertex_key(gi, s, fi, sq); tk <- vertex_key(gi, t, fi, tq)
  ci <- c(seq_along(cohorts), gi, gi, match(targets$g, cohorts))
  time <- c(refs, s, t, targets$time)
  f <- c(rep(0L, length(cohorts)), fi, fi, rep(0L, nrow(targets)))
  qual <- c(rep(TRUE, length(cohorts)), sq, tq, rep(TRUE, nrow(targets)))
  key <- vertex_key(ci, time, f, qual)
  kind <- ifelse(time == refs[ci], 1L, ifelse(time < cohorts[ci], 2L, ifelse(!qual, 3L, 4L)))
  fv <- ifelse(kind %in% c(2L, 3L), f, 0L)
  keep <- !duplicated(key)
  ord <- which(keep)[order(ci[keep], kind[keep], time[keep], fv[keep], method = "radix")]
  vertices <- data.frame(g = cohorts[ci[ord]], time = time[ord],
    functional_key = ifelse(kind[ord] %in% c(2L, 3L), fn[pmax(fv[ord], 1L)], NA_character_),
    role = c("reference", "pre", "post_nuisance", "post")[kind[ord]], vertex_index = seq_along(ord),
    stringsAsFactors = FALSE)
  vk <- key[ord]; nv <- nrow(vertices)
  pv <- which(vertices$role != "reference"); p <- length(pv)
  vertices$param_index <- NA_integer_; vertices$param_index[pv] <- seq_len(p)
  # Unique display labels are not used to construct/match graph vertices.
  label <- ifelse(vertices$role == "reference", paste0("ref(g=", fmt_design_num(vertices$g), ",t=", fmt_design_num(vertices$time), ")"),
    ifelse(vertices$role == "post", paste0("ATT(g=", fmt_design_num(vertices$g), ",t=", fmt_design_num(vertices$time), ")"),
      ifelse(vertices$role == "post_nuisance",
        paste0("beta_post_nuisance(g=", fmt_design_num(vertices$g), ",t=", fmt_design_num(vertices$time), ",fn=", vertices$functional_key, ")"),
        paste0("beta(g=", fmt_design_num(vertices$g), ",t=", fmt_design_num(vertices$time), ",fn=", vertices$functional_key, ")"))))
  vertices$vertex_id <- make.unique(label)
  from <- match(sk, vk); to <- match(tk, vk)
  canonical <- order(gi, s, t, fi, method = "radix"); edge_order <- order(canonical)
  anchors <- which(vertices$role == "reference")
  forest <- incidence_forest(nv, from, to, anchors, edge_order)
  anchored <- forest$component_id[anchors]
  vertices$degree <- forest$degree; vertices$is_isolated <- forest$degree == 0L
  vertices$component_id <- forest$component_id
  vertices$connected_to_ref <- vertices$component_id %in% anchored
  vertices$anchored_value <- ifelse(vertices$role == "reference", 0, NA_real_)
  parameters <- vertices[pv, , drop = FALSE]; rownames(parameters) <- NULL
  parameters$param_id <- parameters$vertex_id
  parameters$level_estimable <- parameters$connected_to_ref
  # One additional linear-size time graph exposes false collapsed connectivity.
  all_ci <- match(vertices$g, cohorts)
  time_key <- paste(all_ci, match(vertices$time, dates), sep = ":")
  ut <- unique(time_key); ti <- match(time_key, ut)
  naive <- incidence_forest(length(ut), ti[from], ti[to], ti[anchors], edge_order)
  parameters$time_graph_connected_to_ref <- naive$component_id[ti[pv]] %in% naive$component_id[ti[anchors]]
  links$link_index <- seq_len(m); links$link_id <- paste0("link_", edge_order)
  links$start_vertex_index <- from; links$end_vertex_index <- to
  links$start_param_index <- vertices$param_index[from]; links$end_param_index <- vertices$param_index[to]
  links$link_type <- ifelse(t < g, "pre_only", ifelse(s == ref | s >= g, "post_only", "straddling"))
  links$component_id <- forest$component_id[from]; rownames(links) <- NULL
  post <- parameters[parameters$role == "post", , drop = FALSE]
  post_targets <- post[, c("g", "time", "param_id", "param_index", "degree", "is_isolated", "component_id", "connected_to_ref", "level_estimable", "time_graph_connected_to_ref"), drop = FALSE]
  post_targets$is_supported <- post_targets$degree > 0L
  post_targets$status <- ifelse(post_targets$level_estimable, "estimable",
    ifelse(post_targets$is_isolated, "isolated_unsupported", "disconnected_component"))
  rownames(post_targets) <- NULL
  unanchored <- setdiff(seq_len(forest$n_components), anchored)
  rank <- as.integer(nv - forest$n_components)
  # Retain inexpensive per-cohort diagnostics without per-cohort dense graphs.
  nh <- length(cohorts)
  pci <- match(parameters$g, cohorts)
  vcount <- tabulate(all_ci, nh); ccount <- tabulate(all_ci[forest$roots], nh)
  pcount <- tabulate(pci, nh); mcount <- tabulate(gi, nh)
  pre_flag <- parameters$role == "pre"
  post_nuis_flag <- parameters$role == "post_nuisance"
  post_flag <- parameters$role == "post"
  by_cohort <- data.frame(g = cohorts, ref = refs, n_links = mcount,
    n_vertices = vcount, n_parameters = pcount,
    n_pre_parameters = tabulate(pci[pre_flag], nh),
    n_post_parameters = tabulate(pci[post_flag], nh), n_components = ccount,
    rank = vcount - ccount, nullity = pcount - (vcount - ccount),
    cycle_rank = mcount - (vcount - ccount),
    n_estimable_post = tabulate(pci[post_flag & parameters$level_estimable], nh),
    n_isolated_post = tabulate(pci[post_flag & parameters$is_isolated], nh),
    n_disconnected_post = tabulate(pci[post_flag & !parameters$level_estimable & !parameters$is_isolated], nh))
  if (has_qual) {
    by_cohort$n_post_nuisance_parameters <- tabulate(pci[post_nuis_flag], nh)
  }
  rank_diag <- list(n_links = m, n_vertices = nv, n_parameters = p,
    n_cohorts = length(cohorts), n_components = forest$n_components, rank = rank,
    nullity = length(unanchored), cycle_rank = m - rank, has_full_column_rank = length(unanchored) == 0L,
    unanchored_components = unanchored, n_post_targets = nrow(post_targets),
    n_estimable_post_targets = sum(post_targets$level_estimable),
    n_isolated_post_targets = sum(post_targets$is_isolated),
    n_disconnected_post_targets = sum(!post_targets$level_estimable & !post_targets$is_isolated),
    by_cohort = by_cohort)
  if (has_qual) {
    rank_diag$n_pre_parameters <- sum(pre_flag)
    rank_diag$n_post_nuisance_parameters <- sum(post_nuis_flag)
    rank_diag$n_post_parameters <- sum(post_flag)
  }
  incidence <- list(start = ifelse(is.na(links$start_param_index), 0L, links$start_param_index),
    end = ifelse(is.na(links$end_param_index), 0L, links$end_param_index),
    n_links = m, n_parameters = p, link_ids = links$link_id, param_ids = parameters$param_id)
  structure(list(incidence = incidence, forest = forest, parameters = parameters,
    vertices = vertices, links = links, post_targets = post_targets,
    rank_diagnostics = rank_diag), class = "cdid_chained_design")
}

construct_chained_design <- build_chained_design

incidence_multiply <- function(design, x, transpose = FALSE) {
  inc <- design$incidence
  n <- if (transpose) inc$n_links else inc$n_parameters
  vector_input <- is.null(dim(x))
  if (!is.numeric(x) || any(!is.finite(x)) ||
      (vector_input && (length(x) != n || is.object(x))) ||
      (!vector_input && (!is.matrix(x) || nrow(x) != n))) cdid_abort("Invalid incidence product input.", "cdid_design_input_error")
  if (vector_input) x <- matrix(x, ncol = 1L)
  if (!transpose) {
    grounded <- rbind(rep(0, ncol(x)), x)
    out <- grounded[inc$end + 1L, , drop = FALSE] - grounded[inc$start + 1L, , drop = FALSE]
    rownames(out) <- inc$link_ids
  } else {
    out <- matrix(0, inc$n_parameters, ncol(x))
    # rowsum aggregates each endpoint without allocating an m-by-p object.
    for (side in c("start", "end")) {
      j <- inc[[side]]; keep <- j > 0L
      if (any(keep)) {
        sums <- rowsum(x[keep, , drop = FALSE], j[keep], reorder = FALSE)
        dest <- as.integer(rownames(sums))
        out[dest, ] <- out[dest, , drop = FALSE] + (if (side == "start") -1 else 1) * sums
      }
    }
    rownames(out) <- inc$param_ids
  }
  if (any(!is.finite(out))) cdid_abort("Incidence product overflow.", "cdid_design_product_error")
  colnames(out) <- colnames(x)
  if (vector_input) drop(out) else out
}

resolve_design_parameter <- function(design, x) {
  if (is.character(x) && is.null(dim(x)) && length(x) == 1L && !is.na(x)) {
    j <- match(x, design$parameters$param_id)
  } else if (design_vector(x, 1L) && x == floor(x) && x >= 1 && x <= nrow(design$parameters)) {
    j <- as.integer(x)
  } else j <- NA_integer_
  if (is.na(j)) cdid_abort("Unknown or noninteger parameter index.", "cdid_design_input_error")
  j
}

# Sparse signed root path, constructed only when requested. Length <= vertices.
root_path_certificate <- function(design, vertex) {
  f <- design$forest; k <- 0L; cur <- vertex
  while (f$parent_vertex[cur] != 0L) { k <- k + 1L; cur <- f$parent_vertex[cur] }
  edges <- signs <- integer(k); cur <- vertex
  for (j in seq_len(k)) {
    edges[j] <- f$parent_link[cur]; signs[j] <- f$parent_sign[cur]; cur <- f$parent_vertex[cur]
  }
  ord <- order(edges)
  data.frame(link_index = edges[ord], sign = signs[ord])
}

check_level_estimable <- function(design, param) {
  j <- resolve_design_parameter(design, param)
  ok <- design$parameters$level_estimable[j]
  list(estimable = ok, witness_links = if (ok) root_path_certificate(design, design$parameters$vertex_index[j]) else NULL)
}

check_contrast_estimable <- function(design, param_a, param_b) {
  if (!inherits(design, "cdid_chained_design")) cdid_abort("Expected chained design.", "cdid_design_input_error")
  a <- resolve_design_parameter(design, param_a); b <- resolve_design_parameter(design, param_b)
  same <- design$parameters$component_id[a] == design$parameters$component_id[b]
  both <- design$parameters$level_estimable[a] && design$parameters$level_estimable[b]
  ok <- same || both; witness <- NULL
  if (ok) {
    aa <- root_path_certificate(design, design$parameters$vertex_index[a])
    bb <- root_path_certificate(design, design$parameters$vertex_index[b])
    edges <- sort(unique(c(aa$link_index, bb$link_index)))
    sign <- integer(length(edges))
    sign[match(aa$link_index, edges)] <- -aa$sign
    sign[match(bb$link_index, edges)] <- sign[match(bb$link_index, edges)] + bb$sign
    witness <- data.frame(link_index = edges[sign != 0L], sign = sign[sign != 0L])
  }
  list(param_a = design$parameters$param_id[a], param_b = design$parameters$param_id[b],
    estimable = ok, same_component = same, both_level_estimable = both, witness_links = witness)
}

nullspace_certificate <- function(design, component) {
  if (!design_vector(component, 1L) || component != floor(component) ||
      !component %in% design$rank_diagnostics$unanchored_components) cdid_abort("Expected unanchored component index.", "cdid_design_input_error")
  list(component_id = component, param_index = which(design$parameters$component_id == component), value = 1L)
}

# Debug/test reference only: check the total dense entries BEFORE any allocation.
materialize_chained_design <- function(design, max_entries) {
  bound <- design_bound(max_entries)
  inc <- design$incidence; m <- as.double(inc$n_links); p <- as.double(inc$n_parameters)
  nv <- as.double(nrow(design$vertices)); q <- as.double(design$rank_diagnostics$nullity)
  terms <- c(m * p, m * p, m * nv, p * q)
  if (any(!is.finite(terms)) || any(terms > bound) || sum(terms) > bound) cdid_abort("Dense reference exceeds explicit entry bound.", "cdid_design_bound_error")
  W <- matrix(0L, m, p, dimnames = list(inc$link_ids, inc$param_ids))
  k <- which(inc$start > 0L); W[cbind(k, inc$start[k])] <- -1L
  k <- which(inc$end > 0L); W[cbind(k, inc$end[k])] <- 1L
  roots <- matrix(0L, m, nv)
  for (v in seq_len(nv)) {
    cert <- root_path_certificate(design, v)
    roots[cert$link_index, v] <- cert$sign
  }
  paths <- matrix(0L, m, p, dimnames = dimnames(W))
  k <- which(design$parameters$level_estimable)
  paths[, k] <- roots[, design$parameters$vertex_index[k], drop = FALSE]
  B <- matrix(0L, p, q, dimnames = list(inc$param_ids, if (q > 0) paste0("comp_", design$rank_diagnostics$unanchored_components) else character()))
  for (j in seq_len(q)) B[nullspace_certificate(design, design$rank_diagnostics$unanchored_components[j])$param_index, j] <- 1L
  out <- design; out$W <- W
  out$rank_diagnostics$root_paths <- roots; out$rank_diagnostics$path_witnesses <- paths
  out$rank_diagnostics$nullspace_basis <- B
  out
}

build_unaugmented_time_design <- function(links, target_post_times = NULL) {
  v <- validate_design_links(links)
  v <- v[, c("g", "s", "t", "ref"), drop = FALSE]
  v$functional_key <- "unaugmented_shared_time"
  dates <- sort(unique(c(v$s, v$t)))
  key <- paste(match(v$g, sort(unique(v$g))), match(v$s, dates), match(v$t, dates), sep = ":")
  out <- build_chained_design(v[!duplicated(key), , drop = FALSE], target_post_times)
  out$is_valid_augmented <- FALSE
  out
}
