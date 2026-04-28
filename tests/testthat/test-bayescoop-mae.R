make_toy_mae <- function(n = 12, p1 = 20, p2 = 15, seed = 1, y_col = "Y") {
  set.seed(seed)
  sample_ids <- paste0("S", seq_len(n))
  x1 <- matrix(
    rnorm(p1 * n),
    nrow = p1,
    dimnames = list(paste0("g", seq_len(p1)), sample_ids)
  )
  x2 <- matrix(
    rnorm(p2 * n),
    nrow = p2,
    dimnames = list(paste0("p", seq_len(p2)), sample_ids)
  )

  col_data <- S4Vectors::DataFrame(row.names = sample_ids)
  col_data[[y_col]] <- rnorm(n)

  MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(rna = x1, prot = x2),
    colData = col_data
  )
}

test_that("MAP fit and predict work on MAE input", {
  mae <- make_toy_mae(seed = 11)

  fit <- suppressWarnings(suppressMessages(BayesCOOP(
    data_train = mae,
    bb = FALSE,
    control = list(rho = c(0, 0.5)),
    maxit = 5,
    warning = FALSE,
    verbose = FALSE
  )))

  expect_s3_class(fit, "BayesCOOP")
  expect_type(fit$rho_MAP, "double")
  expect_true(length(fit$beta_MAP) > 0)

  pred <- suppressWarnings(suppressMessages(predict(
    object = fit,
    newdata = mae,
    bb = FALSE,
    warning = FALSE,
    verbose = FALSE
  )))

  n <- nrow(SummarizedExperiment::colData(mae))
  expect_equal(length(pred$y_pred), n)
  expect_equal(length(pred$y_valid), n)
})

test_that("BB fit and predict return posterior quantities", {
  mae <- make_toy_mae(seed = 22)

  fit <- suppressWarnings(suppressMessages(BayesCOOP(
    data_train = mae,
    bb = TRUE,
    bbiters = 6,
    bbburn = 2,
    maxit = 3,
    warning = FALSE,
    verbose = FALSE
  )))

  expect_equal(nrow(fit$beta_samples), 4)
  expect_equal(length(fit$rho_samples), 4)

  pred <- suppressWarnings(suppressMessages(predict(
    object = fit,
    newdata = mae,
    bb = TRUE,
    warning = FALSE,
    verbose = FALSE
  )))

  n <- nrow(SummarizedExperiment::colData(mae))
  expect_equal(length(pred$y_pred), n)
  expect_equal(length(pred$y_valid), n)
})

test_that("validation errors are raised for unsupported MAE inputs", {
  mae_one <- make_toy_mae(seed = 33)
  mae_one <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = list(rna = MultiAssayExperiment::experiments(mae_one)[[1]]),
    colData = SummarizedExperiment::colData(mae_one)
  )

  expect_error(
    suppressWarnings(suppressMessages(BayesCOOP(
      data_train = mae_one,
      bb = FALSE,
      control = list(rho = c(0, 0.5)),
      maxit = 3,
      warning = FALSE,
      verbose = FALSE
    ))),
    "at least two assays"
  )

  mae_no_y <- make_toy_mae(seed = 44, y_col = "outcome")
  expect_error(
    suppressWarnings(suppressMessages(BayesCOOP(
      data_train = mae_no_y,
      bb = FALSE,
      control = list(rho = c(0, 0.5)),
      maxit = 3,
      warning = FALSE,
      verbose = FALSE
    ))),
    "Response column `Y` is missing"
  )
})
