.bhf_small_article_synthetic <- function(...) {
  arguments <- list(
    profile = "low", rho = 0.5, n_states = 3L, n_strata = 2L,
    psus_per_stratum = c(2L, 3L), observations_per_cell = 1L,
    vhat_state = c(0.0002, 0.0003, 0.0004),
    effect_seed = 7121L, outcome_seed = 7122L
  )
  do.call(generate_article_synthetic, utils::modifyList(arguments, list(...)))
}

test_that("article low and high profiles freeze published simulation sigmas", {
  low <- generate_article_synthetic(profile = "low", rho = 0)
  high <- generate_article_synthetic(profile = "high", rho = 0.5)
  expect_identical(low$truth$config$alpha, -1.5)
  expect_identical(low$truth$sigmas,
                   c(state = 0.3, stratum = 0.4, psu = 0.5))
  expect_identical(high$truth$sigmas,
                   c(state = 0.6, stratum = 0.7, psu = 0.9))
  expect_identical(low$truth$config$rho, 0)
  expect_identical(high$truth$config$rho, 0.5)
})

test_that("generator is reproducible and restores the caller RNG state", {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(91)
  before <- .Random.seed
  first <- .bhf_small_article_synthetic()
  expect_identical(.Random.seed, before)
  second <- .bhf_small_article_synthetic()
  expect_identical(first, second)

  changed_outcome <- .bhf_small_article_synthetic(outcome_seed = 9001L)
  expect_identical(first$truth$effects, changed_outcome$truth$effects)
  expect_identical(first$data$probability, changed_outcome$data$probability)
  expect_false(identical(first$data$outcome, changed_outcome$data$outcome))
})

test_that("states are crossed, PSUs nested, and every stratum has two PSUs", {
  generated <- .bhf_small_article_synthetic()
  data <- generated$data
  truth <- generated$truth
  crossing <- xtabs(~ state_id + stratum_id, data = data)
  expect_true(all(crossing > 0))
  psu_to_stratum <- split(data$stratum_id, data$psu_flat_id)
  expect_true(all(vapply(psu_to_stratum, function(x) {
    length(unique(x)) == 1L
  }, logical(1))))
  expect_true(all(truth$config$psus_per_stratum >= 2L))
  expect_identical(
    sort(unique(truth$psu_structure$psu_flat_id)),
    seq_len(truth$config$dimensions[["J"]])
  )
})

test_that("all article centering constraints hold at machine precision", {
  generated <- .bhf_small_article_synthetic()
  truth <- generated$truth
  expect_equal(
    sum(truth$population_shares * truth$effects$state), 0,
    tolerance = 1e-15
  )
  expect_equal(mean(truth$effects$stratum), 0, tolerance = 1e-15)
  for (h in seq_len(truth$config$dimensions[["H"]])) {
    block <- truth$psu_structure$stratum_id == h
    expect_true(sum(block) >= 2L)
    expect_equal(mean(truth$effects$psu[block]), 0, tolerance = 1e-15)
  }
  expect_equal(truth$centering_residuals$state_weighted, 0,
               tolerance = 1e-15)
  expect_equal(truth$centering_residuals$stratum_mean, 0,
               tolerance = 1e-15)
  expect_equal(truth$centering_residuals$psu_within_stratum,
               rep(0, truth$config$dimensions[["H"]]), tolerance = 1e-15)
})

test_that("informative tilt and global mean-one weight truth are exact", {
  uninformative <- generate_article_synthetic(profile = "low", rho = 0)
  informative <- .bhf_small_article_synthetic()
  expect_identical(uninformative$data$raw_weight,
                   rep(1, nrow(uninformative$data)))
  expected_tilt <- exp(
    -0.5 * informative$truth$effects$stratum[informative$data$stratum_id]
  )
  expect_equal(informative$data$raw_weight, unname(expected_tilt),
               tolerance = 0)
  expect_equal(
    informative$data$w_lik / informative$data$raw_weight,
    rep(nrow(informative$data) / sum(informative$data$raw_weight),
        nrow(informative$data)),
    tolerance = 1e-15
  )
  expect_equal(sum(informative$data$w_lik), nrow(informative$data),
               tolerance = 1e-14)
  expect_identical(informative$truth$weights$normalization,
                   "global_mean_one")
})

test_that("binary outcomes and independent A/A-star/B identities are valid", {
  generated <- .bhf_small_article_synthetic()
  truth <- generated$truth
  estimands <- truth$estimands
  expect_true(all(generated$data$outcome %in% c(0L, 1L)))
  expect_equal(generated$data$probability, stats::plogis(generated$data$eta),
               tolerance = 0)
  expect_equal(estimands$A$summary[["total"]],
               estimands$A$summary[["mean"]] *
                 (1 - estimands$A$summary[["mean"]]),
               tolerance = 1e-15)
  expect_equal(estimands$B$summary[["total"]],
               estimands$B$summary[["mean"]] *
                 (1 - estimands$B$summary[["mean"]]),
               tolerance = 1e-15)
  expect_equal(estimands$A_star$summary[["between"]],
               max(0, estimands$A$summary[["between"]] -
                     estimands$A_star$correction),
               tolerance = 1e-15)
  expect_equal(estimands$B$summary[["within"]],
               estimands$B$summary[["within_binomial"]] +
                 estimands$B$summary[["within_mixture"]],
               tolerance = 1e-15)
})

test_that("primitive truth calculation is invariant to row permutation", {
  generated <- .bhf_small_article_synthetic()
  truth <- generated$truth
  order <- .bhf_synthetic_with_seed(
    7008L, function() sample.int(nrow(generated$data))
  )
  permuted <- .bhf_article_truth_from_primitives(
    data = generated$data[order, , drop = FALSE],
    alpha = truth$config$alpha,
    u_state = truth$effects$state,
    u_stratum = truth$effects$stratum,
    u_psu = truth$effects$psu,
    sigmas = truth$sigmas,
    population_shares = truth$population_shares,
    vhat_state = truth$vhat_state
  )
  expect_equal(permuted, truth$estimands, tolerance = 1e-14)
})

test_that("committed text truth fixture and checksums round-trip", {
  fixture_path <- testthat::test_path(
    "..", "fixtures", "article_synthetic_truth.dput"
  )
  fixture <- dget(fixture_path)
  generated <- do.call(
    generate_article_synthetic, fixture$generator_arguments
  )
  expect_identical(fixture$schema_version,
                   "article-synthetic-fixture-1.0.0")

  # The numeric contract is the substantive one: regenerating from the recorded
  # generator arguments must reproduce the committed truth. Compare it without
  # the checksum block, which is compared separately below.
  fixture_truth <- fixture$truth
  generated_truth <- generated$truth
  fixture_truth$checksums <- NULL
  generated_truth$checksums <- NULL
  expect_equal(fixture_truth, generated_truth, tolerance = 1e-14)

  expect_identical(generated$truth$checksums$algorithm,
                   "adler32-canonical-dput")
  expect_identical(names(generated$truth$checksums),
                   names(fixture$truth$checksums))

  # Integer outcomes serialise identically everywhere, so their checksum is a
  # portable contract.
  expect_identical(generated$truth$checksums$outcomes, "7da00c25")

  # The data and truth checksums are a reference-platform fingerprint. dput()
  # writes full double precision, so last-bit arithmetic differences between
  # platforms change the serialised text without changing any value beyond the
  # 1e-14 tolerance asserted above. Pin the minted literals only where the
  # regenerated data reproduces them bit for bit.
  minted <- identical(
    generated$truth$checksums$data, fixture$truth$checksums$data
  )
  if (minted) {
    expect_identical(generated$truth$checksums$data, "c4274f95")
    expect_identical(generated$truth$checksums$truth_without_checksums,
                     "558c3644")
  } else {
    expect_match(generated$truth$checksums$data, "^[0-9a-f]{8}$")
    expect_match(generated$truth$checksums$truth_without_checksums,
                 "^[0-9a-f]{8}$")
  }

  # Self-consistency holds on every platform: the stored checksum must equal a
  # fresh computation over the same object.
  without_checksums <- generated$truth
  without_checksums$checksums <- NULL
  expect_identical(
    .bhf_synthetic_checksum(without_checksums),
    generated$truth$checksums$truth_without_checksums
  )
  expect_false(any(grepl("\\x00", readLines(fixture_path, warn = FALSE),
                         fixed = TRUE)))
})

test_that("regeneration script is read-only by default", {
  script <- testthat::test_path(
    "..", "..", "data-raw", "generate_article_synthetic.R"
  )
  if (!file.exists(script)) {
    succeed("data-raw is intentionally excluded from the source archive")
    return(invisible(NULL))
  }
  fixture <- testthat::test_path(
    "..", "fixtures", "article_synthetic_truth.dput"
  )
  before_hash <- unname(tools::md5sum(fixture))
  before_time <- file.info(fixture)$mtime
  old <- Sys.getenv("BHFVAR_REGENERATE_SYNTHETIC_FIXTURE", unset = NA_character_)
  on.exit({
    if (is.na(old)) Sys.unsetenv("BHFVAR_REGENERATE_SYNTHETIC_FIXTURE") else
      Sys.setenv(BHFVAR_REGENERATE_SYNTHETIC_FIXTURE = old)
  }, add = TRUE)
  Sys.setenv(BHFVAR_REGENERATE_SYNTHETIC_FIXTURE = "false")
  environment <- new.env(parent = globalenv())
  expect_no_error(capture.output(sys.source(script, envir = environment)))
  expect_identical(unname(tools::md5sum(fixture)), before_hash)
  expect_identical(file.info(fixture)$mtime, before_time)
  expect_true(exists("article_synthetic_truth_fixture", envir = environment,
                     inherits = FALSE))
})

test_that("generator rejects non-article and non-nested configurations", {
  expect_error(generate_article_synthetic(profile = "lo"),
               "exactly 'low' or 'high'")
  expect_error(generate_article_synthetic(rho = 0.25),
               "exactly 0 or 0.5")
  expect_error(generate_article_synthetic(alpha = -1),
               "exactly -1.5")
  expect_error(generate_article_synthetic(psus_per_stratum = 1L),
               "at least two")
  expect_error(generate_article_synthetic(n_states = 1L),
               "n_states")
  expect_error(generate_article_synthetic(vhat_state = -1),
               "vhat_state")
})
