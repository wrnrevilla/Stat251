library(readr)
library(dplyr)
library(car)
library(MVN)
library(biotools)
library(ggplot2)
library(GGally)

file_path <- "C:/Users/William/Desktop/R code - Survey ops/Data for G10 Analysis.xlsx - final_data.csv"

df <- read_csv(file_path)

names(df)

g10_complete <- df %>%
  dplyr::select(
    Health_Deficit_Z,
    LeisureSocial_Deficit_Z,
    Household_Deficit_Z,
    Productivity_Deficit_Z,
    FA1_scaled,
    FA2_scaled,
    FA3_scaled,
    FA4_scaled,
    FA5_scaled,
    FA6_scaled,
    g1,
    a1,
    g6,
    g2,
    a2,
    a3,
    g4,
    category
  ) %>%
  dplyr::mutate(
    g1 = factor(g1),
    g6 = factor(g6),
    g2 = factor(g2),
    a2 = factor(a2),
    a3 = factor(a3),
    g4 = factor(g4),
    category = factor(category)
  ) %>%
  na.omit() %>%
  droplevels()

nrow(g10_complete)

sapply(
  g10_complete %>%
    dplyr::select(g1, g6, g2, a2, a3, g4, category),
  nlevels
)

model <- lm(
  cbind(
    Health_Deficit_Z,
    LeisureSocial_Deficit_Z,
    Household_Deficit_Z,
    Productivity_Deficit_Z
  ) ~
    FA1_scaled +
    FA2_scaled +
    FA3_scaled +
    FA4_scaled +
    FA5_scaled +
    FA6_scaled +
    g1 +
    a1 +
    g6 +
    g2 +
    a2 +
    a3 +
    g4 +
    category,
  data = g10_complete
)

summary(model)

manova_model <- manova(model)

pillai_result <- summary(manova_model, test = "Pillai")
wilks_result <- summary(manova_model, test = "Wilks")
hl_result <- summary(manova_model, test = "Hotelling-Lawley")
roy_result <- summary(manova_model, test = "Roy")
anova_result <- summary.aov(manova_model)

pillai_result
wilks_result
hl_result
roy_result
anova_result

vif_model <- lm(
  Health_Deficit_Z ~
    FA1_scaled +
    FA2_scaled +
    FA3_scaled +
    FA4_scaled +
    FA5_scaled +
    FA6_scaled +
    g1 +
    a1 +
    g6 +
    g2 +
    a2 +
    a3 +
    g4 +
    category,
  data = g10_complete
)

vif_result <- car::vif(vif_model)

vif_table <- as.data.frame(vif_result)

if ("GVIF" %in% names(vif_table)) {
  vif_table$Adjusted_GVIF <- vif_table$GVIF^(1 / (2 * vif_table$Df))
}

vif_table

resids <- as.data.frame(residuals(model))

mvn_result <- MVN::mvn(data = resids)

mvn_result$multivariate_normality
mvn_result$univariate_normality
mvn_result$descriptives

Y_boxm <- as.matrix(
  g10_complete %>%
    dplyr::select(
      Health_Deficit_Z,
      LeisureSocial_Deficit_Z,
      Household_Deficit_Z,
      Productivity_Deficit_Z
    )
)

group_boxm <- droplevels(g10_complete$a2)

boxm_result <- biotools::boxM(Y_boxm, group_boxm)

boxm_result

mahal <- mahalanobis(
  resids,
  colMeans(resids),
  cov(resids)
)

mahal_cutoff <- qchisq(0.999, df = 4)

outliers <- which(mahal > mahal_cutoff)

outlier_summary <- data.frame(
  cutoff = mahal_cutoff,
  number_of_outliers = length(outliers),
  percent_outliers = length(outliers) / nrow(g10_complete) * 100,
  max_mahalanobis = max(mahal),
  mean_mahalanobis = mean(mahal)
)

outlier_summary
outliers

par(mfrow = c(2, 2))

qqnorm(resids$Health_Deficit_Z, main = "Q-Q Plot: Health Deficit")
qqline(resids$Health_Deficit_Z)

qqnorm(resids$LeisureSocial_Deficit_Z, main = "Q-Q Plot: Leisure/Social Deficit")
qqline(resids$LeisureSocial_Deficit_Z)

qqnorm(resids$Household_Deficit_Z, main = "Q-Q Plot: Household Deficit")
qqline(resids$Household_Deficit_Z)

qqnorm(resids$Productivity_Deficit_Z, main = "Q-Q Plot: Productivity Deficit")
qqline(resids$Productivity_Deficit_Z)

par(mfrow = c(1, 1))

hist(resids$Health_Deficit_Z, main = "Histogram: Health Deficit Residuals", xlab = "Residuals")
hist(resids$LeisureSocial_Deficit_Z, main = "Histogram: Leisure/Social Deficit Residuals", xlab = "Residuals")
hist(resids$Household_Deficit_Z, main = "Histogram: Household Deficit Residuals", xlab = "Residuals")
hist(resids$Productivity_Deficit_Z, main = "Histogram: Productivity Deficit Residuals", xlab = "Residuals")

boxplot(
  resids$Health_Deficit_Z,
  resids$LeisureSocial_Deficit_Z,
  resids$Household_Deficit_Z,
  resids$Productivity_Deficit_Z,
  names = c("Health", "Leisure/Social", "Household", "Productivity"),
  main = "Boxplots of Residuals",
  ylab = "Residuals"
)

plot(
  mahal,
  main = "Mahalanobis Distance",
  ylab = "Mahalanobis Distance",
  xlab = "Observation"
)

abline(h = mahal_cutoff, lty = 2)

vif_plot_data <- vif_table
vif_plot_data$Variable <- rownames(vif_plot_data)

ggplot(vif_plot_data, aes(x = reorder(Variable, Adjusted_GVIF), y = Adjusted_GVIF)) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Multicollinearity Check: Adjusted GVIF",
    x = "Variable",
    y = "Adjusted GVIF"
  )

fitted_values <- as.data.frame(fitted(model))

par(mfrow = c(2, 2))

plot(
  fitted_values$Health_Deficit_Z,
  resids$Health_Deficit_Z,
  main = "Residuals vs Fitted: Health",
  xlab = "Fitted Values",
  ylab = "Residuals"
)
abline(h = 0, lty = 2)

plot(
  fitted_values$LeisureSocial_Deficit_Z,
  resids$LeisureSocial_Deficit_Z,
  main = "Residuals vs Fitted: Leisure/Social",
  xlab = "Fitted Values",
  ylab = "Residuals"
)
abline(h = 0, lty = 2)

plot(
  fitted_values$Household_Deficit_Z,
  resids$Household_Deficit_Z,
  main = "Residuals vs Fitted: Household",
  xlab = "Fitted Values",
  ylab = "Residuals"
)
abline(h = 0, lty = 2)

plot(
  fitted_values$Productivity_Deficit_Z,
  resids$Productivity_Deficit_Z,
  main = "Residuals vs Fitted: Productivity",
  xlab = "Fitted Values",
  ylab = "Residuals"
)
abline(h = 0, lty = 2)

par(mfrow = c(1, 1))

GGally::ggpairs(
  g10_complete %>%
    dplyr::select(
      Health_Deficit_Z,
      LeisureSocial_Deficit_Z,
      Household_Deficit_Z,
      Productivity_Deficit_Z
    )
)