## Analysis: Figures 3 and 4 Kaplan-Meier analyses
## Date: 2026-09-08
## Random seed: 42

set.seed(42)
source("config.R")
source("R/_common.R")
assert_packages(c("broom", "dplyr", "ggplot2", "patchwork", "survival", "survminer", "tidyr"))

io <- resolve_io()
data <- read_analysis_data(io$input_file)
specifications <- list(
  `figure3A_overall_survival_by_t2el` = list(
    time = "followup_years", event = "all_cause_death", group = "t2el",
    group_labels = c("Non-T2EL", "T2EL"), y = "Overall survival (%)"
  ),
  `figure3B_freedom_reintervention_by_t2el` = list(
    time = "reintervention_followup_years", event = "any_reintervention", group = "t2el",
    group_labels = c("Non-T2EL", "T2EL"), y = "Freedom from reintervention (%)"
  ),
  `figure4B_overall_survival_by_sac_status` = list(
    time = "followup_years", event = "all_cause_death", group = "sac_enlargement",
    group_labels = c("Regression/Stability", "Enlargement"), y = "Overall survival (%)"
  ),
  `figure4C_freedom_reintervention_by_sac_status` = list(
    time = "reintervention_followup_years", event = "any_reintervention", group = "sac_enlargement",
    group_labels = c("Regression/Stability", "Enlargement"), y = "Freedom from reintervention (%)"
  ),
  `figure4D_freedom_t2el_reintervention_by_sac_status` = list(
    time = "t2el_reintervention_followup_years", event = "t2el_related_reintervention",
    group = "sac_enlargement", group_labels = c("Regression/Stability", "Enlargement"),
    y = "Freedom from T2EL-related reintervention (%)"
  )
)
required <- unique(unlist(lapply(specifications, function(x) c(x$time, x$event, x$group))))
require_columns(data, required, "Figures 3 and 4")
dir.create(io$output_dir, recursive = TRUE, showWarnings = FALSE)

draw_survival <- function(stem, specification) {
  d <- data[, c(specification$time, specification$event, specification$group), drop = FALSE]
  names(d) <- c("time", "event", "group")
  d$event <- as_binary(d$event)
  d$group <- factor(as_binary(d$group), levels = c(0, 1), labels = specification$group_labels)
  d <- d[stats::complete.cases(d) & d$time >= 0, , drop = FALSE]
  if (length(unique(d$group)) != 2) stop("Two non-empty groups are required for ", stem)
  fit <- survival::survfit(survival::Surv(time, event) ~ group, data = d)
  cox <- survival::coxph(survival::Surv(time, event) ~ group, data = d)
  cox_row <- broom::tidy(cox, exponentiate = TRUE, conf.int = TRUE)[1, ]
  hr_label <- sprintf(
    "HR %.2f (95%% CI: %.2f–%.2f)\nP %s",
    cox_row$estimate, cox_row$conf.low, cox_row$conf.high,
    ifelse(cox_row$p.value < 0.001, "< .001", paste0("= ", sub("^0", "", sprintf("%.3f", cox_row$p.value))))
  )
  plot <- survminer::ggsurvplot(
    fit, data = d, risk.table = TRUE, conf.int = TRUE, censor = TRUE,
    break.time.by = 3, xlim = range(cfg$risk_times), xlab = "Time (years)",
    ylab = specification$y, surv.scale = "percent",
    palette = c("#1F6D8C", "#A84255"), legend.title = "",
    legend.labs = specification$group_labels,
    ggtheme = theme_rsna(10), tables.theme = theme_rsna(9),
    risk.table.height = 0.25
  )
  plot$plot <- plot$plot + ggplot2::annotate(
    "text", x = 0.5, y = 0., label = hr_label,
    hjust = 0, vjust = 0, family = "Arial", size = 3.2
  )
  arranged <- survminer::arrange_ggsurvplots(list(plot), print = FALSE)
  ggplot2::ggsave(file.path(io$output_dir, paste0(stem, ".pdf")), arranged, width = 7, height = 5.2)
  ggplot2::ggsave(file.path(io$output_dir, paste0(stem, ".png")), arranged, width = 7, height = 5.2, dpi = 300, bg = "white")
  data.frame(
    Figure = stem, N = nrow(d), Events = sum(d$event),
    HR = cox_row$estimate, CI_lower = cox_row$conf.low,
    CI_upper = cox_row$conf.high, P_value = cox_row$p.value
  )
}

survival_results <- dplyr::bind_rows(Map(draw_survival, names(specifications), specifications))
utils::write.csv(survival_results, file.path(io$output_dir, "figures3_4_survival_statistics.csv"), row.names = FALSE)

## Figure 4A: sac status by T2EL group.
bar_data <- data.frame(
  t2el = factor(as_binary(data$t2el), levels = c(0, 1), labels = c("Non-T2EL", "T2EL")),
  sac = factor(as_binary(data$sac_enlargement), levels = c(0, 1), labels = c("Regression/Stability", "Enlargement"))
) |>
  tidyr::drop_na()
counts <- bar_data |>
  dplyr::count(t2el, sac) |>
  dplyr::group_by(t2el) |>
  dplyr::mutate(percent = 100 * n / sum(n), label = sprintf("%d (%.1f%%)", n, percent)) |>
  dplyr::ungroup()
test <- stats::chisq.test(table(bar_data$t2el, bar_data$sac), correct = FALSE)
figure4a <- ggplot2::ggplot(counts, ggplot2::aes(t2el, n, fill = sac)) +
  ggplot2::geom_col(width = 0.65) +
  ggplot2::geom_text(ggplot2::aes(label = label), position = ggplot2::position_stack(vjust = 0.5), size = 3) +
  ggplot2::scale_fill_manual(values = c("Regression/Stability" = "#1F6D8C", "Enlargement" = "#A84255")) +
  ggplot2::labs(x = NULL, y = "Number of patients", fill = "Sac status") +
  ggplot2::annotate("text", x = 1.5, y = max(tapply(counts$n, counts$t2el, sum)) * 0.75,
                    label = paste0("P ", ifelse(test$p.value < 0.001, "< .001", paste0("= ", format_p_rsna(test$p.value)))), size = 3.2) +
  theme_rsna(10)
save_publication_plot(figure4a, "figure4A_sac_status", io$output_dir, 6.5, 4.8)
print(survival_results)
