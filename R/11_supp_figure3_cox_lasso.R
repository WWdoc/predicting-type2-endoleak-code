## Analysis: Supplementary Figure 3 Cox-LASSO sensitivity analysis
## Date: 2026-09-08
## Random seed: 42

set.seed(42)
source("config.R")
source("R/_common.R")
assert_packages(c("dplyr", "ggplot2", "glmnet", "patchwork", "scales", "survival", "tidyr"))

io <- resolve_io()
data <- read_analysis_data(io$input_file)
time_col <- "t2el_followup_years"
event_col <- "t2el"
exposure <- "ostial_area_mm2"
candidate_sets <- lapply(cfg$cox_models[-1], function(x) unique(c(exposure, x)))
require_columns(data, unique(c(time_col, event_col, unlist(candidate_sets))), "Cox-LASSO analysis")
data[[event_col]] <- as_binary(data[[event_col]])
dir.create(io$output_dir, recursive = TRUE, showWarnings = FALSE)

model_titles <- c(
  `Model 2` = "Clinically important baseline variables",
  `Model 3` = "Non-branch anatomical variables",
  `Model 4` = "Baseline plus non-branch anatomical variables",
  `Model 5` = "Comprehensive preoperative variables"
)

run_lasso <- function(model_name, variables) {
  columns <- unique(c(time_col, event_col, variables))
  d <- data[stats::complete.cases(data[, columns, drop = FALSE]), columns, drop = FALSE]
  d <- d[d[[time_col]] > 0, , drop = FALSE]
  x <- stats::model.matrix(~ . - 1, data = d[, variables, drop = FALSE])
  y <- survival::Surv(d[[time_col]], d[[event_col]])
  if (ncol(x) >= sum(d[[event_col]] == 1)) {
    warning(model_name, ": candidate coefficients are numerous relative to events; this is a penalized sensitivity analysis.")
  }
  fold_id <- stratified_folds(d[[event_col]], cfg$cv_folds)
  fit <- glmnet::glmnet(x, y, family = "cox", alpha = 1, standardize = TRUE)
  cv <- glmnet::cv.glmnet(
    x, y, family = "cox", alpha = 1, type.measure = "C",
    foldid = fold_id, standardize = TRUE
  )
  selected_min <- rownames(as.matrix(stats::coef(cv, s = "lambda.min")))[as.matrix(stats::coef(cv, s = "lambda.min"))[, 1] != 0]
  selected_1se <- rownames(as.matrix(stats::coef(cv, s = "lambda.1se")))[as.matrix(stats::coef(cv, s = "lambda.1se"))[, 1] != 0]
  index_min <- which.min(abs(cv$lambda - cv$lambda.min))
  index_1se <- which.min(abs(cv$lambda - cv$lambda.1se))

  beta <- as.data.frame(as.matrix(fit$beta))
  beta$term <- rownames(beta)
  coefficient_path <- beta |>
    tidyr::pivot_longer(-term, names_to = "step", values_to = "coefficient") |>
    dplyr::mutate(
      step = as.integer(sub("V", "", step)), lambda = fit$lambda[step],
      log_lambda = log(lambda), highlight = term == exposure
    )
  path_plot <- ggplot2::ggplot(coefficient_path, ggplot2::aes(log_lambda, coefficient, group = term)) +
    ggplot2::geom_line(ggplot2::aes(color = highlight, linewidth = highlight), alpha = 0.85) +
    ggplot2::geom_vline(xintercept = log(c(cv$lambda.min, cv$lambda.1se)), linetype = "dashed",
                        color = c("#2C7BB6", "#C44E52")) +
    ggplot2::scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "#9B59B6")) +
    ggplot2::scale_linewidth_manual(values = c(`FALSE` = 0.35, `TRUE` = 1.0)) +
    ggplot2::labs(title = "Coefficient path", x = "log(lambda)", y = "Coefficient") +
    theme_rsna(9) + ggplot2::theme(legend.position = "none")

  cv_data <- data.frame(
    log_lambda = log(cv$lambda), c_index = cv$cvm,
    lower = cv$cvm - cv$cvsd, upper = cv$cvm + cv$cvsd
  )
  cv_label <- sprintf(
    "C-index at lambda.min = %.3f (%d variables)\nC-index at lambda.1se = %.3f (%d variables)",
    cv$cvm[index_min], length(selected_min), cv$cvm[index_1se], length(selected_1se)
  )
  cv_plot <- ggplot2::ggplot(cv_data, ggplot2::aes(log_lambda, c_index)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper), color = "grey70", width = 0) +
    ggplot2::geom_line(color = "#2C7BB6") + ggplot2::geom_point(color = "#2C7BB6", size = 1.5) +
    ggplot2::geom_vline(xintercept = log(cv$lambda.min), linetype = "dashed", color = "#2C7BB6") +
    ggplot2::geom_vline(xintercept = log(cv$lambda.1se), linetype = "dashed", color = "#C44E52") +
    ggplot2::annotate("label", x = -Inf, y = -Inf, label = cv_label, hjust = -0.05, vjust = -0.3, size = 2.5) +
    ggplot2::labs(title = "Cross-validated C-index", x = "log(lambda)", y = "C-index") +
    theme_rsna(9)

  selection_count <- setNames(numeric(ncol(x)), colnames(x))
  valid_bootstrap <- 0L
  for (b in seq_len(cfg$bootstrap_repetitions)) {
    index <- sample(seq_len(nrow(d)), nrow(d), replace = TRUE)
    event_b <- d[[event_col]][index]
    if (length(unique(event_b)) < 2 || min(table(event_b)) < 3) next
    fold_b <- stratified_folds(event_b, cfg$cv_folds)
    cv_b <- tryCatch(
      glmnet::cv.glmnet(
        x[index, , drop = FALSE], survival::Surv(d[[time_col]][index], event_b),
        family = "cox", alpha = 1, type.measure = "C", foldid = fold_b, standardize = TRUE
      ), error = function(e) NULL
    )
    if (is.null(cv_b)) next
    coefficients <- as.matrix(stats::coef(cv_b, s = "lambda.1se"))
    selected <- rownames(coefficients)[coefficients[, 1] != 0]
    selection_count[names(selection_count) %in% selected] <- selection_count[names(selection_count) %in% selected] + 1
    valid_bootstrap <- valid_bootstrap + 1L
  }
  if (valid_bootstrap == 0) stop("No valid bootstrap resamples for ", model_name)
  selection <- data.frame(
    term = names(selection_count), frequency = as.numeric(selection_count / valid_bootstrap)
  ) |>
    dplyr::arrange(frequency) |>
    dplyr::mutate(term = factor(term, levels = term), highlight = as.character(term) == exposure)
  selection_plot <- ggplot2::ggplot(selection, ggplot2::aes(term, frequency, fill = highlight)) +
    ggplot2::geom_col(width = 0.7) + ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c(`FALSE` = "grey72", `TRUE` = "#9B59B6")) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    ggplot2::labs(title = "Bootstrap selection frequency", x = NULL, y = "Selection frequency") +
    theme_rsna(9) + ggplot2::theme(legend.position = "none")

  combined <- path_plot + cv_plot + selection_plot +
    patchwork::plot_annotation(title = paste(model_name, model_titles[[model_name]], sep = ": "))
  stem <- paste0("supp_figure3_", gsub(" ", "_", tolower(model_name)))
  save_publication_plot(combined, stem, io$output_dir, 15, 5.5)
  utils::write.csv(coefficient_path, file.path(io$output_dir, paste0(stem, "_coefficient_path.csv")), row.names = FALSE)
  utils::write.csv(cv_data, file.path(io$output_dir, paste0(stem, "_cv_curve.csv")), row.names = FALSE)
  utils::write.csv(selection, file.path(io$output_dir, paste0(stem, "_selection_frequency.csv")), row.names = FALSE)
  data.frame(
    Model = model_name, N = nrow(d), Events = sum(d[[event_col]]), Candidate_coefficients = ncol(x),
    Lambda_min = cv$lambda.min, Lambda_1se = cv$lambda.1se,
    C_index_lambda_min = cv$cvm[index_min], C_index_lambda_1se = cv$cvm[index_1se],
    Selected_lambda_min = paste(selected_min, collapse = "; "),
    Selected_lambda_1se = paste(selected_1se, collapse = "; "),
    Valid_bootstrap = valid_bootstrap, stringsAsFactors = FALSE
  )
}

summary <- dplyr::bind_rows(Map(run_lasso, names(candidate_sets), candidate_sets))
utils::write.csv(summary, file.path(io$output_dir, "supp_figure3_cox_lasso_summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
