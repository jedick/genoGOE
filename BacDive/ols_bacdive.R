# OLS on BacDive data with range parsing

# 1. Read data (keep original column names)
df <- read.csv(
  "cleaned_data.csv",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 2. Keep only rows with complete data for the model variables
vars_needed <- c(
  "GC_content.GC-content",
  "culture_temp.Temperature",
  "culture_pH.pH",
  "oxygen_tolerance.Oxygen tolerance",
  "nutrition_type.Nutrition type"
)

df_model <- df[complete.cases(df[, vars_needed]), vars_needed]

# 3. Fit OLS model
model <- lm(
  `GC_content.GC-content` ~
    `culture_temp.Temperature` +
    `culture_pH.pH` +
    `oxygen_tolerance.Oxygen tolerance` +
    `nutrition_type.Nutrition type`,
  data = df_model
)

# 4. Show results
summary(model)
