## Shared utilities for the public analysis scripts.

set.seed(42)
options(stringsAsFactors = FALSE, scipen = 999)

assert_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing R packages: ", paste(missing, collapse = ", "),
      ". Install them before running this script."
    )
  }
}

resolve_io <- function(default_input = cfg$input_file) {
  args <- commandArgs(trailingOnly = TRUE)
  input_file <- if (length(args) >= 1L && nzchar(args[[1]])) args[[1]] else default_input
  output_dir <- if (length(args) >= 2L && nzchar(args[[2]])) args[[2]] else cfg$output_dir
  list(input_file = input_file, output_dir = output_dir)
}

read_analysis_data <- function(path) {
  if (!file.exists(path)) stop("Input file does not exist: ", path)
  extension <- tolower(tools::file_ext(path))
  if (extension == "csv") {
    assert_packages("readr")
    return(as.data.frame(readr::read_csv(path, show_col_types = FALSE)))
  }
  if (extension %in% c("xlsx", "xls")) {
    assert_packages("readxl")
    return(as.data.frame(readxl::read_excel(path, .name_repair = "unique")))
  }
  stop("Input must be CSV, XLSX, or XLS.")
}

require_columns <- function(data, columns, context = "analysis") {
  columns <- unique(columns[nzchar(columns)])
  missing <- setdiff(columns, names(data))
  if (length(missing)) {
    stop(
      "Missing columns for ", context, ": ", paste(missing, collapse = ", "),
      ". See DATA_DICTIONARY.md or edit config.R."
    )
  }
  invisible(TRUE)
}

as_binary <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) {
    out <- as.numeric(x)
  } else {
    z <- tolower(trimws(as.character(x)))
    out <- rep(NA_real_, length(z))
    out[z %in% c("0", "no", "n", "false", "absent", "negative")] <- 0
    out[z %in% c("1", "yes", "y", "true", "present", "positive", "event")] <- 1
  }
  invalid <- !is.na(out) & !out %in% c(0, 1)
  if (any(invalid)) stop("A binary variable contains values other than 0/1.")
  out
}

prepare_binary_factors <- function(data, columns) {
  require_columns(data, columns, "binary-variable preparation")
  for (column in columns) {
    data[[column]] <- factor(
      as_binary(data[[column]]), levels = c(0, 1), labels = c("No", "Yes")
    )
  }
  data
}

prepare_t2el_group <- function(data) {
  require_columns(data, "t2el", "T2EL grouping")
  data$t2el <- factor(
    as_binary(data$t2el), levels = c(0, 1), labels = c("Non-T2EL", "T2EL")
  )
  data
}

backtick_names <- function(x) paste0("`", gsub("`", "", x, fixed = TRUE), "`")

make_survival_formula <- function(time, event, predictors = character(0)) {
  rhs <- if (length(predictors)) paste(backtick_names(predictors), collapse = " + ") else "1"
  stats::as.formula(
    paste0("survival::Surv(", backtick_names(time), ", ", backtick_names(event), ") ~ ", rhs)
  )
}

formula_list <- function(columns, value) {
  lapply(columns, function(column) {
    rlang::new_formula(rlang::sym(column), value, env = parent.frame())
  })
}

label_formulas <- function(labels) {
  Map(
    function(column, label) rlang::new_formula(rlang::sym(column), label, env = parent.frame()),
    names(labels), unname(labels)
  )
}

format_p_rsna <- function(p) {
  ifelse(
    is.na(p), "NA",
    ifelse(
      p < 0.001, "<.001",
      ifelse(
        p > 0.99, ">.99",
        ifelse(p < 0.01, sub("^0", "", sprintf("%.3f", p)), sub("^0", "", sprintf("%.2f", p)))
      )
    )
  )
}

format_estimate_ci <- function(estimate, lower, upper, digits = 2L) {
  pattern <- paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)")
  sprintf(pattern, estimate, lower, upper)
}

export_gtsummary <- function(table, stem, output_dir, title) {
  assert_packages(c("gtsummary", "flextable", "knitr"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  table_df <- as.data.frame(gtsummary::as_tibble(table, col_labels = FALSE))
  utils::write.csv(
    table_df, file.path(output_dir, paste0(stem, ".csv")),
    row.names = FALSE, fileEncoding = "UTF-8"
  )
  markdown <- c(paste0("# ", title), "", knitr::kable(table_df, format = "pipe"))
  writeLines(markdown, file.path(output_dir, paste0(stem, ".md")), useBytes = TRUE)
  ft <- gtsummary::as_flex_table(table)
  do.call(
    flextable::save_as_docx,
    c(setNames(list(ft), title), list(path = file.path(output_dir, paste0(stem, ".docx"))))
  )
  invisible(table_df)
}

save_publication_plot <- function(plot, stem, output_dir, width = 7, height = 5) {
  assert_packages("ggplot2")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    file.path(output_dir, paste0(stem, ".pdf")), plot = plot,
    width = width, height = height, units = "in", device = "pdf"
  )
  ggplot2::ggsave(
    file.path(output_dir, paste0(stem, ".png")), plot = plot,
    width = width, height = height, units = "in", dpi = 300, bg = "white"
  )
  invisible(plot)
}

theme_rsna <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Arial") +
    ggplot2::theme(
      axis.text = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      legend.title = ggplot2::element_blank()
    )
}

wilson_ci <- function(successes, total, conf_level = 0.95) {
  if (total == 0) return(c(lower = NA_real_, upper = NA_real_))
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  p <- successes / total
  denominator <- 1 + z^2 / total
  center <- (p + z^2 / (2 * total)) / denominator
  half <- z * sqrt((p * (1 - p) + z^2 / (4 * total)) / total) / denominator
  c(lower = max(0, center - half), upper = min(1, center + half))
}

stratified_folds <- function(event, k = 10L) {
  event <- as_binary(event)
  if (sum(event == 1, na.rm = TRUE) < 3 || sum(event == 0, na.rm = TRUE) < 3) {
    stop("At least three events and three non-events are required for cross-validation.")
  }
  k_use <- min(k, sum(event == 1), sum(event == 0))
  fold_id <- integer(length(event))
  for (level in c(0, 1)) {
    index <- which(event == level)
    fold_id[index] <- sample(rep(seq_len(k_use), length.out = length(index)))
  }
  fold_id
}

agreement_statistics <- function(x, y) {
  assert_packages("irr")
  keep <- stats::complete.cases(x, y)
  x <- as.numeric(x[keep])
  y <- as.numeric(y[keep])
  n <- length(x)
  if (n < 3) stop("At least three complete measurement pairs are required.")
  difference <- x - y
  average <- (x + y) / 2
  bias <- mean(difference)
  sd_difference <- stats::sd(difference)
  loa_lower <- bias - 1.96 * sd_difference
  loa_upper <- bias + 1.96 * sd_difference
  bias_se <- sd_difference / sqrt(n)
  t_critical <- stats::qt(0.975, df = n - 1)
  bias_ci <- bias + c(-1, 1) * t_critical * bias_se
  loa_se <- sd_difference * sqrt(1 / n + (1.96^2) / (2 * (n - 1)))
  loa_lower_ci <- loa_lower + c(-1, 1) * t_critical * loa_se
  loa_upper_ci <- loa_upper + c(-1, 1) * t_critical * loa_se
  icc <- irr::icc(
    data.frame(measurement_1 = x, measurement_2 = y),
    model = "twoway", type = "agreement", unit = "single"
  )
  bias_test <- stats::t.test(difference, mu = 0)
  list(
    n = n, x = x, y = y, average = average, difference = difference,
    bias = bias, bias_ci = bias_ci, loa_lower = loa_lower, loa_upper = loa_upper,
    loa_lower_ci = loa_lower_ci, loa_upper_ci = loa_upper_ci,
    bias_p = bias_test$p.value, icc = unname(icc$value),
    icc_lower = icc$lbound, icc_upper = icc$ubound
  )
}

agreement_comparisons <- function(data) {
  require_columns(data, c("W1", "W2", "Z1", "Z2"), "observer-variability analysis")
  definitions <- list(
    `Intra-rater: Observer W (W1 vs W2)` = c("W1", "W2"),
    `Intra-rater: Observer Z (Z1 vs Z2)` = c("Z1", "Z2"),
    `Inter-rater: Time 1 (W1 vs Z1)` = c("W1", "Z1"),
    `Inter-rater: Time 2 (W2 vs Z2)` = c("W2", "Z2")
  )
  lapply(definitions, function(columns) agreement_statistics(data[[columns[1]]], data[[columns[2]]]))
}

agreement_results_frame <- function(results) {
  output <- do.call(rbind, Map(function(label, z) {
    data.frame(
      Comparison = label,
      N = z$n,
      Measurement_1_mean_SD = sprintf("%.2f ± %.2f", mean(z$x), stats::sd(z$x)),
      Measurement_2_mean_SD = sprintf("%.2f ± %.2f", mean(z$y), stats::sd(z$y)),
      Bias_95CI = format_estimate_ci(z$bias, z$bias_ci[1], z$bias_ci[2]),
      LoA = sprintf("%.2f to %.2f", z$loa_lower, z$loa_upper),
      LoA_lower_95CI = sprintf("%.2f to %.2f", z$loa_lower_ci[1], z$loa_lower_ci[2]),
      LoA_upper_95CI = sprintf("%.2f to %.2f", z$loa_upper_ci[1], z$loa_upper_ci[2]),
      ICC_95CI = format_estimate_ci(z$icc, z$icc_lower, z$icc_upper),
      P_value_bias = z$bias_p,
      stringsAsFactors = FALSE
    )
  }, names(results), results))
  output$P_value_bias_Holm <- stats::p.adjust(output$P_value_bias, method = "holm")
  output$P_value_bias_display <- format_p_rsna(output$P_value_bias)
  output$P_value_bias_Holm_display <- format_p_rsna(output$P_value_bias_Holm)
  output
}
