# Regenerate only the small committed article-aligned truth fixture.
# Default execution is read-only; writing requires an explicit opt-in.

root_candidates <- normalizePath(c(".", "..", "../.."), mustWork = FALSE)
root_valid <- vapply(root_candidates, function(path) {
  file.exists(file.path(path, "DESCRIPTION")) &&
    file.exists(file.path(path, "R", "synthetic_dgp.R"))
}, logical(1))
if (!any(root_valid)) stop("Cannot locate the bhfvar package root.")
root <- root_candidates[which(root_valid)[[1L]]]
sys.source(file.path(root, "R", "synthetic_dgp.R"), envir = environment())

fixture_arguments <- list(
  profile = "low", rho = 0.5, n_states = 3L, n_strata = 2L,
  psus_per_stratum = c(2L, 3L), observations_per_cell = 1L,
  vhat_state = c(0.0002, 0.0003, 0.0004),
  effect_seed = 7121L, outcome_seed = 7122L
)
generated <- do.call(generate_article_synthetic, fixture_arguments)
article_synthetic_truth_fixture <- list(
  schema_version = "article-synthetic-fixture-1.0.0",
  generator_arguments = fixture_arguments,
  truth = generated$truth
)
cat("article synthetic fixture checksums:\n")
print(article_synthetic_truth_fixture$truth$checksums)

write_enabled <- tolower(Sys.getenv(
  "BHFVAR_REGENERATE_SYNTHETIC_FIXTURE", "false"
)) %in% c("1", "true", "yes")
if (write_enabled) {
  dput(
    article_synthetic_truth_fixture,
    file = file.path(root, "tests", "fixtures", "article_synthetic_truth.dput")
  )
  cat("wrote the small text fixture; no binary or large result was created.\n")
} else {
  cat("read-only mode; set BHFVAR_REGENERATE_SYNTHETIC_FIXTURE=true to write.\n")
}
