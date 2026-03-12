# Set input / output paths
input_file  <- "custom_download_bacdive_2026-03-11.csv"
output_file <- "cleaned_data.csv"

# Read full file as text
all_lines <- readLines(input_file, warn = FALSE)

# Find header line ("strains.ID_strains,...")
header_idx <- which(grepl("^strains.ID_strains,", all_lines))[1]
if (is.na(header_idx)) stop("Header line with 'strains.ID_strains' not found.")

# Work only on header + data below it
lines <- all_lines[header_idx:length(all_lines)]

# Remove all "References" blocks (from line 'References' through the next blank line)
clean_lines <- character()
i <- 1
n <- length(lines)

while (i <= n) {
  if (grepl("^References$", lines[i])) {
    # Skip "References" line and everything until the next completely empty line
    i <- i + 1
    while (i <= n && nzchar(lines[i])) {
      i <- i + 1
    }
    # Also skip that blank separator line (if present)
    if (i <= n && !nzchar(lines[i])) i <- i + 1
  } else {
    clean_lines <- c(clean_lines, lines[i])
    i <- i + 1
  }
}

# Parse cleaned CSV text into a data.frame
con <- textConnection(clean_lines)

dat <- read.csv(con, stringsAsFactors = FALSE, check.names = FALSE)
close(con)

# Safety: drop rows with missing / empty strain ID
id_col <- "strains.ID_strains"
if (!id_col %in% names(dat)) stop("Column 'strains.ID_strains' not found after import.")

dat <- dat[!is.na(dat[[id_col]]) & dat[[id_col]] != "", ]

# Keep only the first row per strain (first occurrence in file)
dat_first <- dat[!duplicated(dat[[id_col]]), ]

# Write result
write.csv(dat_first, output_file, row.names = FALSE)

# Optional: show how many strains/rows were kept
cat("Number of unique strains:", nrow(dat_first), "\n")
cat("Output written to:", output_file, "\n")
