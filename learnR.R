filter_table <- function(df, matching_values, key = "human.entrez") {
  if (!key %in% colnames(df)) {
    stop(paste("Column '", key, "' not found in the data frame."))
  }

  matching_rows <- df[df[, key] %in% matching_values, ]
  return(matching_rows)
}

# Example usage:
data <- data.frame(
  "human.entrez" = c(1, 2, 3),
  "psd" = c("A", "B", "C"),
  "pre" = c("X", "Y", "X"),
  "syn" = c("Y", "Z", "Z")
)

matching_values <- c(1, 3)
result <- filter_table(data, matching_values)

print(result)






