# New D27-only tests. Pure base R entry points; no generator calls on sourcing.
# Independent finite references are packaged in fixtures/cp02-gap/; see its
# provenance.dcf for the pinned sources, exporter and byte hashes.
# Replay is explicitly external: a frozen replay script owns the 16-call budget.

# Locators follow test-cp02-continuous-dgp.R: CDID_CP02_DEV_DIR, then the source
# tree (tests/testthat/../../inst/development), then the installed package
# (R CMD check, installed tests). Fixtures travel with tests/testthat. Missing
# or altered inputs are errors, never skips.
cp02_gap_dev_files <- c("seed.R", "dgp.R", "cp02_baseline_variants_dgp.R", "cp02_gap_dgp.R")
cp02_gap_fixture_files <- c("catalog.csv", "rows.csv", "coordinates.csv", "targets.csv",
                            "references.csv", "truth.csv", "restricted.csv")

cp02_gap_dev_dir <- function(test_dir, installed = system.file("development", package = "cdid"),
                             override = Sys.getenv("CDID_CP02_DEV_DIR", "")) {
  complete <- function(d) nzchar(d) && all(file.exists(file.path(d, cp02_gap_dev_files)))
  if (nzchar(override)) {
    if (!complete(override)) stop("CDID_CP02_DEV_DIR lacks the D27 gap development files.", call. = FALSE)
    return(normalizePath(override))
  }
  for (d in c(file.path(test_dir, "..", "..", "inst", "development"), installed))
    if (complete(d)) return(normalizePath(d))
  stop("inst/development/cp02_gap_dgp.R and its accepted helpers not found.", call. = FALSE)
}

# MD5 of the exported LF bytes (CR removed), checked with base R tools only.
cp02_gap_fixture_dir <- function(test_dir) {
  dir <- file.path(test_dir, "fixtures", "cp02-gap")
  prov <- file.path(dir, "provenance.dcf")
  if (!file.exists(prov)) stop("Packaged D27 gap fixtures not found: ", dir, call. = FALSE)
  p <- read.dcf(prov); p <- p[p[, "Record"] == "fixture", , drop = FALSE]
  if (!identical(sort(p[, "File"]), sort(cp02_gap_fixture_files)) ||
      !all(file.exists(file.path(dir, p[, "File"]))))
    stop("Packaged D27 gap fixtures are incomplete: ", dir, call. = FALSE)
  md5 <- vapply(file.path(dir, p[, "File"]), function(f) {
    b <- readBin(f, "raw", file.size(f)); tmp <- tempfile(); on.exit(unlink(tmp))
    writeBin(b[b != as.raw(13L)], tmp); unname(tools::md5sum(tmp))
  }, "")
  if (!identical(unname(md5), unname(p[, "MD5"])))
    stop("Packaged D27 gap fixture bytes differ from provenance.dcf.", call. = FALSE)
  normalizePath(dir)
}

cp02_gap_test_load <- function(dev_dir) {
  e <- new.env(parent = baseenv())
  for (f in cp02_gap_dev_files) sys.source(file.path(dev_dir, f), envir = e)
  e
}

cp02_gap_read_fixtures <- function(dir) {
  lapply(stats::setNames(file.path(dir, cp02_gap_fixture_files), sub("[.]csv$", "", cp02_gap_fixture_files)),
    utils::read.csv, stringsAsFactors = FALSE, colClasses = "character")
}

cp02_gap_test_assert <- function(ok, label) {
  if (!isTRUE(ok)) stop(paste("D27 check failed:", label), call. = FALSE)
}
cp02_gap_test_equal <- function(a, b, label, tol = 1e-12) {
  if (inherits(a, "table")) a <- unclass(a)
  if (inherits(b, "table")) b <- unclass(b)
  cp02_gap_test_assert(isTRUE(all.equal(unname(a), unname(b), tolerance = tol,
                                      check.attributes = FALSE)), label)
}
cp02_gap_test_error <- function(expr, label) {
  ans <- tryCatch({ force(expr); NULL }, error = identity)
  cp02_gap_test_assert(inherits(ans, "error"), label)
  invisible(ans)
}
cp02_gap_test_rng <- function() {
  list(kind = RNGkind(), seed = if (exists(".Random.seed", globalenv(), inherits = FALSE))
    get(".Random.seed", globalenv(), inherits = FALSE) else NULL)
}

# Independent literal laws and means, from law_map.json (not new settings).
cp02_gap_test_spec <- function(law) {
  cal <- startsWith(law, "CAL_"); no <- startsWith(law, "NO_NT_")
  periods <- if (cal) c(0, 3, 7, 9, 12) else 1:5
  gs <- if (cal) c(6, 10, 12) else c(3, 4, 5)
  codes <- if (no) gs else c(0, gs)
  joint <- if (no) rbind(rep(2/15, 3), rep(1/5, 3)) else
    rbind(c(1/5, 1/10, 1/10, 1/10), c(1/20, 3/20, 3/20, 3/20))
  colnames(joint) <- codes
  masks <- if (cal) c(lapply(seq_len(10), function(i) utils::combn(periods, 2)[, i]),
                         list(7, numeric())) else if (startsWith(law, "MON_")) lapply(2:5, seq_len) else
    if (startsWith(law, "NONADJ_")) list(c(1,3),c(1,4),c(1,5),c(2,4),c(2,5),c(3,5)) else list(1:5)
  probs <- if (cal) c(rep(2/25, 10), 1/10, 1/10) else rep(1/length(masks), length(masks))
  list(cal = periods, gs = gs, codes = codes, joint = joint,
    shift = cal || startsWith(law, "NONADJ_"),
    ant = if (cal) if (law == "CAL_B") 1 else 2 else 0,
    skew = endsWith(law, "_S"), masks = masks, probs = probs)
}

cp02_gap_test_formula <- function(law, g, t, x, z = 0, U = 0, raw = 0) {
  r <- cp02_gap_test_spec(law); j <- match(t, r$cal); gi <- match(g, c(0, r$gs))
  alpha <- c(0, 1/2, -1/2, 1)[gi]
  shift <- if (r$shift) (t == min(r$cal)) *
    (c(0,2,-1,1)[gi] + c(0,7/10,1/5,-2/5)[gi] * x) else rep(0, length(t))
  scale <- if (r$skew) (2+j)/4 else rep(1, length(t))
  y0 <- alpha + j/4 + x*(t - min(r$cal))/2 + z + U + scale*raw + shift
  theta <- c(1,2,3/2)[match(g, r$gs)]
  effect <- ifelse(g != 0 & t >= g, theta + (t-g)/2 + x,
    ifelse(g != 0 & t < g & t >= g-r$ant, 1/4+x/2, 0))
  data.frame(y0 = y0, y_g = y0+effect, effect = effect, eps = scale*raw,
    eps_raw = rep_len(raw, length(t)), eps_scale = scale, shift = shift)
}

cp02_gap_test_operator <- function(law, g, pool, a, b) {
  r <- cp02_gap_test_spec(law)
  if (!length(pool)) return(NA_real_)
  val <- 0
  for (x in 0:1) {
    # Bayes mixture is evaluated from independent joint-cell probabilities.
    probs <- r$joint[x+1, as.character(pool)]
    inc <- cp02_gap_test_formula(law, pool, b, x)$y_g -
      cp02_gap_test_formula(law, pool, a, x)$y_g
    target <- cp02_gap_test_formula(law, g, b, x)$y_g -
      cp02_gap_test_formula(law, g, a, x)$y_g
    val <- val + r$joint[x+1, as.character(g)] / sum(r$joint[, as.character(g)]) *
      (target - sum(probs*inc)/sum(probs))
  }
  val
}

cp02_gap_test_potential <- function(E, law) {
  r <- cp02_gap_test_spec(law)
  grid <- expand.grid(G = r$codes, t = r$cal, x = 0:1)
  # Deterministic nonzero static and period innovations test both paths.
  z <- rep(c(-3/5, 7/5), length.out=nrow(grid))
  U <- seq_len(nrow(grid))/17; raw <- rep(c(-1, 0, 2), length.out=nrow(grid))
  before <- serialize(list(grid,z,U,raw), NULL)
  got <- E$cp02_gap_potential(law, grid$G, grid$t, grid$x, z, U, raw)
  want <- cp02_gap_test_formula(law, grid$G, grid$t, grid$x, z, U, raw)
  cp02_gap_test_equal(got, want, paste(law,"potential formula"))
  cp02_gap_test_assert(identical(before, serialize(list(grid,z,U,raw),NULL)), "input mutation")
  # Common conditional untreated increments from actual references, not date1.
  for (g in r$gs) for (x in 0:1) for (u in r$cal[r$cal > min(r$cal)]) {
    ref <- max(r$cal[r$cal+r$ant < g])
    inc <- E$cp02_gap_potential(law, r$codes, u, x)$y0 -
      E$cp02_gap_potential(law, r$codes, ref, x)$y0
    cp02_gap_test_equal(inc, rep((match(u,r$cal)-match(ref,r$cal))/4 + x*(u-ref)/2,
      length(r$codes)), "common actual-baseline untreated increments")
  }
  invisible(TRUE)
}

# Full clone of the loaded binding with one textual body substitution, so the
# change also reaches internal callers. The substitution must actually apply.
cp02_gap_test_mutant <- function(E, fun, from, to) {
  M <- new.env(parent = parent.env(E))
  for (nm in ls(E, all.names = TRUE)) {
    f <- get(nm, envir = E, inherits = FALSE)
    if (is.function(f) && identical(environment(f), E)) environment(f) <- M
    assign(nm, f, envir = M)
  }
  txt <- paste(deparse(body(M[[fun]]), width.cutoff = 500L), collapse = "\n")
  changed <- sub(from, to, txt, fixed = TRUE)
  cp02_gap_test_assert(!identical(txt, changed), paste("mutation actually applied:", fun))
  f <- M[[fun]]; body(f) <- parse(text = changed)[[1L]]; assign(fun, f, envir = M)
  M
}

# One control's D10/D13 population catalog against the original-source fixture
# and the independent operator formula.
cp02_gap_test_catalog <- function(E, law, ctl, ref, got = E$cp02_gap_comparison_truth(law, ctl)) {
  r <- cp02_gap_test_spec(law); s <- E$cp02_gap_settings(law)
  rr <- ref$catalog[ref$catalog$law==law & ref$catalog$control==ctl,]
  cp02_gap_test_equal(c(nrow(got$rows),nrow(got$coordinates),got$rank),
    as.numeric(rr[,c("rows","parameters","rank")]),"independent catalog size/rank")
  rr <- ref$rows[ref$rows$law==law & ref$rows$control==ctl,]
  rows_key <- function(a) paste(a$group,a$s,a$t,a$pool,sep="|")
  # Pool strings from the independent Python reference have 0 first.
  cp02_gap_test_assert(identical(sort(rows_key(got$rows)),sort(rows_key(rr))),"D10/D13 exact row/pool support")
  k <- match(rows_key(rr),rows_key(got$rows)); blank <- function(v) { v[is.na(v)] <- ""; v }
  cp02_gap_test_assert(identical(blank(got$rows$from[k]),rr$from) && identical(blank(got$rows$to[k]),rr$to),
    "original-source endpoint coordinates")
  v <- got$coordinates; rc <- ref$coordinates[ref$coordinates$law==law & ref$coordinates$control==ctl,]
  cp02_gap_test_assert(nrow(v)==nrow(rc) && setequal(v$coordinate,rc$coordinate),"original-source coordinate set")
  cp02_gap_test_assert(identical(v$available[match(rc$coordinate,v$coordinate)],rc$available=="1"),
    "coordinate availability")
  rr <- ref$references[ref$references$law==law & ref$references$control==ctl,]
  cp02_gap_test_assert(nrow(rr)==length(s$reference),"reference set")
  cp02_gap_test_equal(unname(s$reference[rr$group]),as.numeric(rr$reference),"literal actual references")
  rows <- got$rows
  for(i in seq_len(nrow(rows))) {
    pool <- as.numeric(strsplit(rows$pool[i],",",fixed=TRUE)[[1]])
    bs <- cp02_gap_test_operator(law,rows$group[i],pool,rows$reference[i],rows$s[i])
    bt <- cp02_gap_test_operator(law,rows$group[i],pool,rows$reference[i],rows$t[i])
    cp02_gap_test_equal(c(rows$beta_s[i],rows$beta_t[i],rows$mu[i]),c(bs,bt,bt-bs),"operator-specific endpoint/contrast truths")
    # Original and fixed horizon membership checked separately.
    nt <- r$codes[r$codes==0]
    orig <- c(nt,r$gs[r$gs>rows$group[i] & r$gs>rows$t[i]+r$ant])
    if(ctl=="NT") orig <- nt
    cp02_gap_test_assert(identical(rows$original[i],setequal(pool,orig)),"original row flag")
    hs <- if(ctl=="NT") numeric() else r$cal[r$cal>=rows$group[i] & r$cal>=rows$t[i]]
    hs <- hs[vapply(hs,function(h) setequal(pool,c(nt,r$gs[r$gs>rows$group[i] & r$gs>h+r$ant])),TRUE)]
    cp02_gap_test_equal(rows$horizons[[i]],hs,"D13 horizon provenance")
  }
  cp02_gap_test_equal(as.vector(got$incidence %*% v$beta),rows$mu,"all incidence contrasts")
  # Every coordinate has its own population truth, including private pre.
  for(i in which(grepl("private",v$role))) {
    pool <- as.numeric(strsplit(v$pool[i],",",fixed=TRUE)[[1]])
    cp02_gap_test_equal(v$beta[i],cp02_gap_test_operator(law,v$group[i],pool,
      s$reference[as.character(v$group[i])],v$t[i]),"private endpoint truth")
  }
  invisible(TRUE)
}

cp02_gap_run_deterministic_tests <- function(fixture_dir, E) {
  before_rng <- cp02_gap_test_rng()
  ref <- cp02_gap_read_fixtures(fixture_dir)
  fraction <- function(s) vapply(strsplit(s,"/",fixed=TRUE), function(v)
    if(length(v)==1L) as.numeric(v) else as.numeric(v[1])/as.numeric(v[2]), 0)
  laws <- c("CAL_B","CAL_S","MON_B","MON_S","NONADJ_B","NONADJ_S","NO_NT_B","NO_NT_S")
  cp02_gap_test_assert(identical(E$cp02_gap_catalog()$law,laws),"exact eight-law order")
  for (law in laws) {
    r <- cp02_gap_test_spec(law); s <- E$cp02_gap_settings(law)
    cp02_gap_test_equal(s$periods,r$cal,"literal calendar")
    cp02_gap_test_equal(s$cohort_codes,r$codes,"cohort support")
    cp02_gap_test_equal(s$joint,r$joint,"joint population")
    cp02_gap_test_equal(s$anticipation,r$ant,"anticipation in calendar units")
    cp02_gap_test_equal(s$patterns$prob,r$probs,"mask probability")
    cp02_gap_test_equal(s$patterns$n_obs,lengths(r$masks),"mask size")
    want_mask <- t(vapply(r$masks,function(v) r$cal %in% v,logical(5)))
    cp02_gap_test_equal(s$observation_mask,want_mask,"literal masks")
    # Finite midpoint enumeration proves inverse maps' exact probabilities.
    den <- if (startsWith(law,"NO_NT_")) 15 else 20
    pop <- E$cp02_gap_population((seq_len(den)-1/2)/den,law)
    tab <- table(factor(pop$x,0:1),factor(pop$G,r$codes))/den
    cp02_gap_test_equal(tab,r$joint,"direct finite categorical probabilities")
    denp <- if(startsWith(law,"CAL_")) 50 else length(r$masks)
    idx <- E$cp02_gap_pattern_index((seq_len(denp)-1/2)/denp,law)
    cp02_gap_test_equal(tabulate(idx,nbins=length(r$masks))/denp,r$probs,"finite mask map")
    # Deterministic product measure: every (X,G,mask) cell, never sample shares.
    prod <- outer(as.vector(tab),tabulate(idx,nbins=length(r$masks))/denp)
    cp02_gap_test_equal(prod,outer(as.vector(r$joint),r$probs),"independent joint-mask measure")
    cp02_gap_test_equal(E$cp02_gap_population(c(0,1),law)$G,c(r$codes[1],tail(r$codes,1)),"categorical endpoint map")
    cp02_gap_test_potential(E,law)
    tr <- E$cp02_gap_truth(law)
    cp02_gap_test_equal(tr$cohort$pi,colSums(r$joint),"population cohort shares")
    cp02_gap_test_equal(tr$cohort$mean_x,r$joint[2,]/colSums(r$joint),"population E[X|g]")
    fam <- if(startsWith(law,"CAL_")) "CAL" else "BASE"
    rt <- ref$truth[ref$truth$family==fam,]
    key <- function(g,t) paste(g,t,sep=",")
    cp02_gap_test_assert(identical(key(tr$att_gt$group,tr$att_gt$t),key(rt$group,rt$t)),"exact target set/order")
    cp02_gap_test_equal(tr$att_gt$att,fraction(rt$rational),"independent rational ATT")
    if(fam=="CAL") {
      ra <- ref$truth[ref$truth$family=="CAL_anticipation",]
      cp02_gap_test_equal(tr$anticipation[,c("group","t","att")],data.frame(group=as.numeric(ra$group),
        t=as.numeric(ra$t),att=fraction(ra$rational)),"separately labeled causal anticipation")
      cp02_gap_test_assert(all(tr$att_gt$e == tr$att_gt$t-tr$att_gt$group),"literal event time")
    } else cp02_gap_test_assert(nrow(tr$anticipation)==0,"no new anticipation")
    for(ctl in c("NT","NYT")) {
      cp02_gap_test_catalog(E,law,ctl,ref,tr$comparison[[ctl]])
      rr <- ref$targets[ref$targets$law==law & ref$targets$control==ctl,]
      k <- match(key(rr$group,rr$t),key(tr$att_gt$group,tr$att_gt$t))
      cp02_gap_test_assert(!anyNA(k) && length(k)==nrow(tr$att_gt),"target availability set")
      cp02_gap_test_assert(identical(tr$att_gt[[paste0("available_",ctl)]][k],rr$available=="1"),"target availability")
    }
    # Independently recompute BOTH reporting rules, including early pools.
    for(i in seq_len(nrow(tr$reporting))) {
      a <- tr$reporting[i,]; j <- match(a$t,r$cal); reference <- max(r$cal[r$cal+r$ant<a$group])
      base <- if(a$base_period=="varying" && a$t<a$group) r$cal[j-1] else reference
      cutoff <- max(base,a$t)+r$ant
      K <- if(a$control=="NT") r$codes[r$codes==0] else
        r$codes[r$codes==0 | (r$codes!=a$group & r$codes>cutoff)]
      normalized <- a$base_period=="universal" && a$t==reference
      want <- if(normalized) 0 else cp02_gap_test_operator(law,a$group,K,base,a$t)
      cp02_gap_test_equal(c(a$base,a$cutoff),c(base,cutoff),"reporting dates/cutoff")
      cp02_gap_test_assert(a$pool==paste(K,collapse=",") && a$operator_defined==(length(K)>0),"reporting actual pool")
      cp02_gap_test_equal(a$truth,want,"varying/universal population truth")
      prob <- sum(r$probs[vapply(r$masks,function(mask) all(c(base,a$t) %in% mask),TRUE)])
      cp02_gap_test_assert(a$direct_pair_supported == (normalized || (length(K)>0 && prob>0)),"direct reporting pair support")
    }
    cp02_gap_test_assert(nrow(tr$reporting)==54L,"complete both-baseline reporting grid")
    # Independent population aggregate formulas and map support.
    check_agg <- function(values,weights,cells,restricted=FALSE) {
      for(i in seq_len(nrow(values))) {
        a <- values[i,]; m <- weights[weights$aggregate==a$aggregate & weights$component==a$component,]
        kk <- match(key(m$group,m$t),key(cells$group,cells$t))
        cp02_gap_test_assert(!anyNA(kk) && all(m$weight>0),"aggregate selector support")
        cp02_gap_test_equal(sum(m$weight),1,"aggregate weights sum")
        cp02_gap_test_equal(a$att,sum(m$weight*cells$att[kk]),"aggregate map truth")
        cp02_gap_test_assert(a$available==all(cells$available[kk]),"aggregate availability flag")
        # Analytic D5 weight per selected cell, independent of implementation.
        pi <- colSums(r$joint)[as.character(cells$group)]
        use <- if(a$component=="overall") seq_len(nrow(cells)) else if(a$aggregate=="group")
          which(cells$group==as.numeric(sub("g=","",a$component,fixed=TRUE))) else if(a$aggregate=="dynamic")
          which(cells$e==as.numeric(sub("e=","",a$component,fixed=TRUE))) else
          which(cells$t==as.numeric(sub("t=","",a$component,fixed=TRUE)))
        b <- cells[use,]; p <- pi[use]
        wantw <- if(a$aggregate=="simple") p/sum(p) else if(a$aggregate=="group") {
          gps <- unique(b$group); masses <- colSums(r$joint)[as.character(gps)]
          p / sum(masses) / as.numeric(table(b$group)[as.character(b$group)])
        } else {
          lev <- if(a$aggregate=="dynamic") b$e else b$t
          p / as.numeric(tapply(p,lev,sum)[as.character(lev)]) / length(unique(lev))
        }
        matchw <- match(key(m$group,m$t),key(b$group,b$t))
        cp02_gap_test_assert(nrow(m)==nrow(b) && !anyNA(matchw),"exact aggregate target set")
        cp02_gap_test_equal(m$weight,wantw[matchw],"analytic population aggregate weights")
      }
    }
    check_agg(tr$aggregates,tr$aggregate_weights,tr$att_gt)
    overall <- tr$aggregates[tr$aggregates$component=="overall",]
    want <- if(fam=="CAL") c(31/10,269/90,31/10,259/90) else c(47/20,47/20,73/30,131/60)
    cp02_gap_test_equal(overall$att,want,"independent aggregate rational truths")
    if(startsWith(law,"NO_NT_")) {
      cp02_gap_test_assert(all(!overall$available) && sum(tr$att_gt$available)==3L,"noNT unavailable full/terminal levels")
      res <- tr$restricted
      cp02_gap_test_assert(!any(res$weights$t==5) && all(res$values$available),"explicit restricted aggregate selector")
      cp02_gap_test_assert(all(res$weights$target_set!="full_post"),"restricted labels distinct")
      check_agg(res$values,res$weights,tr$att_gt[tr$att_gt$available,],TRUE)
      ov <- res$values[res$values$component=="overall",]
      k <- match(ifelse(ov$aggregate=="simple","simple",paste0(ov$aggregate,"_overall")),ref$restricted$aggregate)
      cp02_gap_test_assert(!anyNA(k) && nrow(ov)==nrow(ref$restricted),"restricted aggregate set")
      cp02_gap_test_equal(ov$att,fraction(ref$restricted$rational[k]),"pinned restricted aggregate rationals")
    } else cp02_gap_test_assert(is.null(tr$restricted),"no implicit aggregate restriction")
    # Exact Gaussian/exponential moments; no distribution sample testing.
    scale <- if(r$skew) (2+1:5)/4 else rep(1,5)
    cp02_gap_test_equal(tr$noise$var_eps,scale^2,"innovation variance")
    cp02_gap_test_equal(tr$noise$third_eps,if(r$skew) 2*scale^3 else rep(0,5),"innovation third moment")
    cp02_gap_test_equal(tr$noise$fourth_eps,if(r$skew) 9*scale^4 else rep(3,5),"innovation fourth moment")
    s$joint[] <- 0
    cp02_gap_test_equal(E$cp02_gap_settings(law)$joint,r$joint,"settings isolation")
  }
  # Accepted exponential quantile identity at fixed probabilities, no RNG.
  u <- c(0,1/10,1/2,9/10)
  cp02_gap_test_equal(E$cp02_centered_exp_quantile(u),-log(1-u)-1,"accepted exponential quantile")
  for(n in list(0,-1,1.5,NA_real_,Inf,TRUE,"20",numeric(),c(1,2),2^31))
    cp02_gap_test_error(E$cp02_gap_validate_draw_inputs("MON_B",n,702103),"invalid n before draws")
  for(seed in list(NA_real_,Inf,1.5,TRUE,"2",numeric(),c(1,2),2^31))
    cp02_gap_test_error(E$cp02_gap_validate_draw_inputs("MON_B",200,seed),"invalid seed before draws")
  for(law in list("CAL",NA_character_,c("MON_B","MON_S"),1))
    cp02_gap_test_error(E$cp02_gap_settings(law),"invalid law")
  # The accepted sim_spec guard is not weakened: CAL's literal design still fails it.
  cp02_gap_test_error(E$sim_spec(periods=c(0,3,7,9,12),cohorts=c(6,10,12)),"sim_spec guard retained")
  # Meaningful formula mutations of a full clone; each independent check must reject.
  muts <- list(
    list(label="halved binary-X post effect", fun="cp02_gap_potential", from="s$kappa * x[post]",
      to="0.5 * s$kappa * x[post]", law="MON_B", check="potential", expect="MON_B potential formula"),
    list(label="first-date shift on the untreated path only", fun="cp02_gap_potential",
      from="y_g = y0 + effect", to="y_g = y0 - shift + effect", law="NONADJ_B", check="potential",
      expect="NONADJ_B potential formula"),
    list(label="anticipation window beyond its literal dates", fun="cp02_gap_potential",
      from="period >= G - s$anticipation", to="period >= G - 2 * s$anticipation", law="CAL_S",
      check="potential", expect="CAL_S potential formula"),
    list(label="pool weights ignoring X-conditional composition", fun="cp02_gap_operator_truth",
      from="weights <- s$joint[x + 1L, as.character(pool)]", to="weights <- colSums(s$joint)[as.character(pool)]",
      law="NONADJ_B", check="catalog", expect="operator-specific endpoint/contrast truths"))
  errs <- vapply(muts, function(m) {
    M <- cp02_gap_test_mutant(E, m$fun, m$from, m$to)
    err <- cp02_gap_test_error(if (m$check == "potential") cp02_gap_test_potential(M, m$law) else
      cp02_gap_test_catalog(M, m$law, "NYT", ref), paste(m$label, "rejected"))
    cp02_gap_test_assert(grepl(m$expect, conditionMessage(err), fixed = TRUE),
      paste(m$label, "fails its independent check"))
    conditionMessage(err)
  }, "")
  cp02_gap_test_assert(identical(before_rng,cp02_gap_test_rng()),"deterministic checks consume no RNG")
  invisible(list(laws=length(laws), catalog_cases=16L,
    mutations=vapply(muts, `[[`, "", "label"), mutation_errors=unname(errs), generator_calls=0L))
}

# Exercise both locators on throwaway copies of the source-tree, R CMD check
# and installed-test layouts. Explicit `installed`/`override` arguments stand
# in for system.file() and the environment; no package or namespace is faked.
cp02_gap_test_locators <- function(dev_dir, fixture_dir) {
  root <- tempfile("cp02-gap-layout-"); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  put <- function(to, files, from) {
    dir.create(to, recursive = TRUE)
    cp02_gap_test_assert(all(file.copy(file.path(from, files), to)), "layout copy")
    normalizePath(to)
  }
  tests <- function(top) file.path(root, top, "tests", "testthat")
  src_t <- tests("src/cdid"); chk_t <- tests("cdid.Rcheck"); ins_t <- tests("lib/cdid")
  for (t in c(src_t, chk_t, ins_t))
    put(file.path(t, "fixtures", "cp02-gap"), c(cp02_gap_fixture_files, "provenance.dcf"), fixture_dir)
  src_dev <- put(file.path(root, "src/cdid/inst/development"), cp02_gap_dev_files, dev_dir)
  chk_dev <- put(file.path(root, "cdid.Rcheck/cdid/development"), cp02_gap_dev_files, dev_dir)
  ins_dev <- put(file.path(root, "lib/cdid/development"), cp02_gap_dev_files, dev_dir)
  loc <- function(t, installed = "", override = "") cp02_gap_dev_dir(t, installed, override)
  cp02_gap_test_assert(identical(loc(src_t), src_dev), "source-tree development locator")
  cp02_gap_test_assert(identical(loc(chk_t, chk_dev), chk_dev), "R CMD check installed locator")
  cp02_gap_test_assert(identical(loc(ins_t, ins_dev), ins_dev), "installed-tests locator")
  cp02_gap_test_assert(identical(loc(src_t, ins_dev), src_dev), "source tree before installed copy")
  cp02_gap_test_assert(identical(loc(chk_t, chk_dev, ins_dev), ins_dev), "explicit override honoured")
  for (t in c(src_t, chk_t, ins_t))
    cp02_gap_test_assert(identical(cp02_gap_fixture_dir(t),
      normalizePath(file.path(t, "fixtures", "cp02-gap"))), "packaged fixture locator")
  f <- file.path(src_t, "fixtures", "cp02-gap", "truth.csv")
  writeLines(readLines(f), f, sep = "\r\n")
  cp02_gap_test_assert(identical(cp02_gap_fixture_dir(src_t),
    normalizePath(file.path(src_t, "fixtures", "cp02-gap"))), "CRLF checkout accepted")
  # Missing or altered inputs fail; nothing is skipped.
  bad <- put(file.path(root, "incomplete"), setdiff(cp02_gap_dev_files, "cp02_gap_dgp.R"), dev_dir)
  fails <- list(
    list(function() loc(chk_t), "not found"), list(function() loc(chk_t, bad), "not found"),
    list(function() loc(src_t, ins_dev, bad), "CDID_CP02_DEV_DIR"),
    list(function() {
      cat("CAL,6,7,0\n", file = file.path(chk_t, "fixtures", "cp02-gap", "truth.csv"), append = TRUE)
      cp02_gap_fixture_dir(chk_t)
    }, "differ"),
    list(function() {
      unlink(file.path(ins_t, "fixtures", "cp02-gap", "rows.csv")); cp02_gap_fixture_dir(ins_t)
    }, "incomplete"),
    list(function() cp02_gap_fixture_dir(file.path(root, "nowhere")), "not found"))
  for (x in fails) {
    err <- cp02_gap_test_error(x[[1]](), paste("locator refusal:", x[[2]]))
    cp02_gap_test_assert(grepl(x[[2]], conditionMessage(err), fixed = TRUE), paste("locator refusal:", x[[2]]))
  }
  invisible(c("source", "check", "installed"))
}

# Check saved replay outputs, never call a constructor here. With the loaded
# binding E, also require the sample-free truth and declared settings.
cp02_gap_check_saved_draw <- function(d, law, n=200L, seed, E=NULL) {
  r <- cp02_gap_test_spec(law); p <- d$data$panel; o <- d$oracle$panel; u <- d$data$units
  cp02_gap_test_assert(identical(names(d),c("data","oracle","settings")),"draw top-level style")
  cp02_gap_test_assert(nrow(p)==5*n && nrow(o)==5*n && nrow(u)==n,"full-frame rows")
  cp02_gap_test_assert(identical(u$id,seq_len(n)),"all frame units")
  cp02_gap_test_equal(p$period,rep(r$cal,n),"literal draw dates")
  cp02_gap_test_assert(identical(p$id,rep(seq_len(n),each=5L)),"calendar unit order")
  cp02_gap_test_assert(all(u$G %in% r$codes) && all(u$x %in% 0:1),"finite draw support")
  cp02_gap_test_assert(d$settings$law==law && d$settings$n==n && d$settings$seed==seed,"draw inputs")
  labels <- vapply(r$masks,function(v) if(length(v)) paste(v,collapse="-") else "none","")
  want_mask <- t(vapply(u$pattern,function(a) r$cal %in% r$masks[[match(a,labels)]],logical(5)))
  cp02_gap_test_equal(d$data$observation_mask,want_mask,"saved observation mask")
  cp02_gap_test_assert(identical(p$observed,as.vector(t(want_mask))),"observed rows")
  cp02_gap_test_equal(u$n_obs,rowSums(want_mask),"unit observation counts")
  cp02_gap_test_assert(all(u$eligible_ge2==(u$n_obs>=2)),"auxiliary eligibility")
  cp02_gap_test_assert(identical(is.na(p$y),!p$observed),"missing-outcome mask")
  cp02_gap_test_assert(identical(d$data$observed,p[p$observed,,drop=FALSE]),"observed frame selector")
  cp02_gap_test_assert(!any(c("y0","y_g","effect","eps","eps_raw","unit_effect") %in%
    c(names(p),names(u))),"oracle isolation")
  for(v in c("G","x","z","pattern","n_obs"))
    cp02_gap_test_assert(identical(p[[v]],u[[v]][p$id]),paste("unit invariant",v))
  cp02_gap_test_assert(identical(p$id,o$id) && identical(p$period,o$period),"potential path alignment")
  want <- cp02_gap_test_formula(law,p$G,p$period,p$x,p$z,d$oracle$units$unit_effect[p$id],o$eps_raw)
  cp02_gap_test_equal(o[,names(want)],want,"saved draw potential formulas")
  cp02_gap_test_equal(p$y[p$observed],want$y_g[p$observed],"observed outcome path")
  cnt <- table(factor(u$G,r$codes),factor(u$pattern,labels))
  cn <- d$data$counts
  cp02_gap_test_equal(cn$units,cnt[cbind(match(cn$G,r$codes),match(cn$pattern,labels))],"realized unforced counts")
  cp02_gap_test_assert(sum(cn$units)==n && nrow(cn)==length(r$codes)*length(labels),"all count cells")
  if (!is.null(E)) {
    cp02_gap_test_assert(identical(d$oracle$truth,E$cp02_gap_truth(law)),"sample-free oracle truth")
    s <- E$cp02_gap_settings(law); s$n <- as.integer(n); s$seed <- seed; s$rng$r_version <- R.version.string
    cp02_gap_test_assert(identical(d$settings,s),"declared settings")
  }
  invisible(TRUE)
}

# testthat::is_testing() without loading testthat when sourced by plain R.
if (identical(Sys.getenv("TESTTHAT"), "true")) {
  test_that("D27 gap locators resolve source, check and installed layouts; missing inputs fail", {
    td <- test_path()
    expect_identical(cp02_gap_test_locators(cp02_gap_dev_dir(td), cp02_gap_fixture_dir(td)),
                     c("source", "check", "installed"))
  })
  test_that("D27 gap bindings satisfy independent deterministic finite references", {
    td <- test_path()
    res <- cp02_gap_run_deterministic_tests(cp02_gap_fixture_dir(td),
                                            cp02_gap_test_load(cp02_gap_dev_dir(td)))
    expect_identical(res$laws, 8L)
    expect_identical(res$generator_calls, 0L)
    expect_length(res$mutation_errors, 4L)
  })
}
