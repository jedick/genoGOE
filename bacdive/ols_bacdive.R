# OLS on BacDive data with range parsing

# Helper: convert "a-b" (optionally with spaces) to mean (a+b)/2
range_to_mean <- function(x) {
  # Treat empty strings as NA
  x <- ifelse(trimws(x) == "", NA, x)

  # If entry contains "-", treat as range; otherwise try simple numeric
  is_range <- grepl("-", x)

  out <- rep(NA_real_, length(x))

  # Handle ranges
  if (any(is_range, na.rm = TRUE)) {
    parts <- strsplit(gsub(" ", "", x[is_range]), "-", fixed = TRUE)
    nums <- lapply(parts, function(v) as.numeric(v[1:2]))
    means <- vapply(nums, function(v) mean(v, na.rm = TRUE), numeric(1))
    out[is_range] <- means
  }

  # Handle non-range numeric values
  if (any(!is_range & !is.na(x))) {
    out[!is_range & !is.na(x)] <- suppressWarnings(
      as.numeric(x[!is_range & !is.na(x)])
    )
  }

  out
}

# 1. Read data (keep original column names)
df <- read.csv(
  "cleaned_data.csv",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# 2. Parse continuous variables (overwrite as numeric means)
df$`culture_temp.Temperature`    <- range_to_mean(df$`culture_temp.Temperature`)
df$`culture_pH.pH`               <- range_to_mean(df$`culture_pH.pH`)
df$`GC_content.GC-content`       <- range_to_mean(df$`GC_content.GC-content`)

# 3. Convert categorical predictors to factors
df$`oxygen_tolerance.Oxygen tolerance` <-
  as.factor(df$`oxygen_tolerance.Oxygen tolerance`)

df$`nutrition_type.Nutrition type` <-
  as.factor(df$`nutrition_type.Nutrition type`)

# 4. Keep only rows with complete data for the model variables
vars_needed <- c(
  "GC_content.GC-content",
  "culture_temp.Temperature",
  "culture_pH.pH",
  "oxygen_tolerance.Oxygen tolerance",
  "nutrition_type.Nutrition type"
)

df_model <- df[complete.cases(df[, vars_needed]), vars_needed]

# 5. Fit OLS model
model <- lm(
  `GC_content.GC-content` ~
    `culture_temp.Temperature` +
    `culture_pH.pH` +
    `oxygen_tolerance.Oxygen tolerance` +
    `nutrition_type.Nutrition type`,
  data = df_model
)

# 6. Show results
summary(model)
