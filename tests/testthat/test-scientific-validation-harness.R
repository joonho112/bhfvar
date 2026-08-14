test_that("Gate G1 tiny and medium plans are frozen exactly", {
  plan <- bhf_scientific_validation_plan()
  expect_identical(plan$tiny$mcmc,
                   list(chains=4L, iter=1000L, warmup=500L,
                        kept_per_chain=500L, adapt_delta=.95,
                        max_treedepth=12L))
  expect_identical(plan$tiny$thresholds$ebfmi_min, .20)
  expect_identical(plan$tiny$thresholds$rhat_max, 1.05)
  expect_identical(plan$tiny$thresholds$ess_bulk_min, 100)
  expect_identical(nrow(plan$medium$jobs), 12L)
  expect_setequal(plan$medium$jobs$profile, c("low","high"))
  expect_setequal(plan$medium$jobs$rho, c(0,.5))
  expect_identical(plan$medium$mcmc$iter, 2000L)
  expect_identical(plan$medium$thresholds$core_rhat_max, 1.01)
  expect_identical(plan$medium$thresholds$ess_bulk_min, 400)
  expect_identical(plan$medium$thresholds$overall_coverage_min, .80)
  expect_true(plan$constraints$oracle_immutable)
  expect_false(plan$constraints$restricted_application_reproduction_claim)
})

test_that("validation materialization uses DGP and current data contract only", {
  plan <- bhf_scientific_validation_plan()
  material <- suppressWarnings(bhf_materialize_validation_job(
    plan$tiny,
    list(n_states=4L,n_strata=3L,psus_per_stratum=2L,
         observations_per_cell=1L)
  ))
  expect_s3_class(material$synthetic, "bhf_article_synthetic")
  expect_s3_class(material$prepared, "bhf_data")
  expect_identical(material$prepared$schema_version, "0.5.0")
  expect_equal(sum(material$prepared$stan_data$w_lik),
               material$prepared$stan_data$N, tolerance=1e-12)
  expect_identical(material$prepared$input_info$deattenuation, "supplied")
  expect_identical(material$prepared$input_info$sigma_state_prior,
                   "half_t3_2.5")
})

test_that("sensitivity grid varies only sigma_state prior", {
  grid <- bhf_scientific_sensitivity_grid()
  expect_identical(grid$variant,
                   c("half_t3_2.5","half_normal_1","half_cauchy_2.5","half_t3_5"))
  expect_identical(grid$sigma_state_prior_code, 1:4)
  expect_true(all(grid$varied_component == "sigma_state"))
  expect_true(all(grid$other_priors_fixed))
})

test_that("diagnostic adjudication uses frozen thresholds", {
  core <- data.frame(variable=c("alpha","sigma_state"),rhat=c(1,1.01),
                     ess_bulk=c(200,150),ess_tail=c(180,140))
  raw <- data.frame(variable="z_state[1]",rhat=1.02,ess_bulk=100,ess_tail=100)
  good <- list(tier="tiny",divergences=0L,treedepth_hits=0L,
               ebfmi=rep(.5,4),core=core,raw_effect=raw)
  expect_true(bhf_adjudicate_validation_diagnostics(good)$pass)
  bad <- good; bad$ebfmi[[2]] <- .1
  verdict <- bhf_adjudicate_validation_diagnostics(bad)
  expect_false(verdict$pass)
  expect_false(unname(verdict$checks[["ebfmi"]]))
})

test_that("artifact status is resumable only for exact complete config", {
  path <- tempfile(fileext=".rds")
  expect_identical(bhf_validation_artifact_status(path,"abc")$status,"missing")
  saveRDS(list(status="complete",config_hash="abc",value=1,
               model_sha256="model-a"),path)
  expect_true(bhf_validation_artifact_status(path,"abc")$reusable)
  expect_true(bhf_validation_artifact_status(
    path,"abc",list(model_sha256="model-a"))$reusable
  )
  stale <- bhf_validation_artifact_status(
    path,"abc",list(model_sha256="model-b")
  )
  expect_false(stale$reusable)
  expect_identical(stale$mismatch,"model_sha256")
  expect_false(bhf_validation_artifact_status(path,"other")$reusable)
  expect_error(
    bhf_validation_artifact_status(path,"abc",list("model-a")),
    "fully named"
  )
  saveRDS(list(status="running",config_hash="abc"),path)
  expect_false(bhf_validation_artifact_status(path,"abc")$reusable)
})

test_that("validation hash fallback hashes serialized content", {
  fallback <- .bhf_validation_hash
  environment(fallback) <- new.env(parent=baseenv())
  first <- fallback(list(a=1L,b="x"))
  second <- fallback(list(a=1L,b="x"))
  changed <- fallback(list(a=2L,b="x"))
  expect_match(first,"^[0-9a-f]{32}$")
  expect_identical(first,second)
  expect_false(identical(first,changed))
})
