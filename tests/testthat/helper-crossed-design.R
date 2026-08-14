# Deterministic microdata used to test the crossed/nested data contract.
#
# This fixture is intentionally assembled without package preparation helpers.
# Its hand-specified mappings remain independent of the current implementation,
# including the legacy assumption that raw PSU labels are globally unique.
make_tiny_crossed_design_fixture <- function() {
  data <- data.frame(
    row_id = sprintf("row%02d", seq_len(12)),
    outcome = c(0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 1L, 1L, 1L, 0L),
    state = rep(c("A", "B", "C"), times = 4),
    stratum = rep(c("H1", "H1", "H2", "H2"), each = 3),
    psu = rep(c("P1", "P2", "P1", "P2"), each = 3),
    weight = as.double(seq_len(12)),
    stringsAsFactors = FALSE
  )

  state_levels <- c("A", "B", "C")
  stratum_levels <- c("H1", "H2")
  composite_psu_levels <- c("H1::P1", "H1::P2", "H2::P1", "H2::P2")
  composite_psu_key <- paste(data$stratum, data$psu, sep = "::")

  state_id <- as.integer(match(data$state, state_levels))
  stratum_id <- as.integer(match(data$stratum, stratum_levels))
  psu_flat_id <- as.integer(match(composite_psu_key, composite_psu_levels))
  psu_within_stratum_id <- rep(rep(1:2, each = 3), times = 2)

  raw_weights <- setNames(data$weight, data$row_id)
  mean_one_weights <- setNames(data$weight / mean(data$weight), data$row_id)

  known <- list(
    population_shares = c(A = 0.50, B = 0.30, C = 0.20),
    vhat = c(A = 0.0025, B = 0.0040, C = 0.0064)
  )

  truth <- list(
    dimensions = c(N = 12L, S = 3L, H = 2L, J = 4L),
    state_levels = state_levels,
    stratum_levels = stratum_levels,
    composite_psu_levels = composite_psu_levels,
    state_id = setNames(state_id, data$row_id),
    stratum_id = setNames(stratum_id, data$row_id),
    composite_psu_key = setNames(composite_psu_key, data$row_id),
    psu_flat_id = setNames(psu_flat_id, data$row_id),
    psu_within_stratum_id = setNames(
      as.integer(psu_within_stratum_id),
      data$row_id
    ),
    composite_psu_map = data.frame(
      composite_key = composite_psu_levels,
      stratum = rep(stratum_levels, each = 2),
      psu = rep(c("P1", "P2"), times = 2),
      flat_id = 1:4,
      within_stratum_id = rep(1:2, times = 2),
      stringsAsFactors = FALSE
    ),
    J_h = c(H1 = 2L, H2 = 2L),
    psu_start = c(H1 = 1L, H2 = 3L),
    raw_weights = raw_weights,
    mean_one_weights = mean_one_weights
  )

  list(data = data, known = known, truth = truth)
}
