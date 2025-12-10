## Loading and combining BirdNET output .csv files

```
#By Jack A. Greenhalgh. December, 2025.
#Department of Biology, McGill University, 1205 Dr Penfield Ave, Montreal, Quebec, H3A 1B1, Canada.

library(dplyr)
library(readr)

##### Part 1. Loading and combining BirdNET output (Revised) #####

# Set the main directory containing all your data folders
# It's often better to use R Projects to manage working directories
setwd("D:/DATA - ALL COMBINED/BirdNET Output Combined")

# Define column types for consistent reading
column_types <- cols(
  `Start (s)` = col_double(),
  `End (s)` = col_double(),
  `Scientific name` = col_character(),
  `Common name` = col_character(),
  `Confidence` = col_double(),
  `File` = col_character(),
  .default = col_character()
)

# Get folder names in the working directory
folders <- list.dirs(path = ".", full.names = FALSE, recursive = FALSE)

# Best practice: Initialize a list to store the results
all_folder_data <- list()

# Loop through each folder
for (folder in folders) {
  folder_path <- file.path(getwd(), folder)
  
  # List all .csv files in the folder
  csv_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
  
  # Skip if no CSV files are found in the folder
  if (length(csv_files) == 0) {
    cat("No CSV files found in folder:", folder, "\n")
    next
  }
  
  # 1. Read and combine all .csv files in one step
  combined_data <- lapply(csv_files, function(file) {
    read_csv(file, col_types = column_types, show_col_types = FALSE) %>%
      mutate(source_file = basename(file))
  }) %>%
    bind_rows()
  
  # 2. Select only the columns you want to keep by NAME (more robust)
  #    This avoids the unnecessary write/read cycle and fragile indexing.
  #    Adjust the names in `select()` to match the columns you need.
  cleaned_data <- combined_data %>%
    select(
      `Start (s)`,
      `End (s)`,
      `Scientific name`,
      `Common name`,
      `Confidence`,
      `File`,
      source_file
    )
  
  # Create a filename and path for the final combined output
  output_filename <- paste0(folder, "_Combined_Results.csv")
  output_path <- file.path(folder_path, output_filename)
  
  # 3. Write the *cleaned* data to the final CSV file
  write_csv(cleaned_data, output_path)
  
  # 4. Store the cleaned data frame in the list (best practice)
  all_folder_data[[folder]] <- cleaned_data
  
  cat("Processed and saved data for folder:", folder, "\n")
}

# List all CSV files in the directory
csv_files <- list.files(pattern = "_Combined_Results.csv$", full.names = TRUE)

# Read each CSV into a list of data frames
data_list <- lapply(csv_files, read.csv, stringsAsFactors = FALSE)

# Name each list element by the file name (without extension)
names(data_list) <- tools::file_path_sans_ext(basename(csv_files))

# Combine all data frames into one (handles mismatched columns)
combined_data <- bind_rows(data_list)

combined_data <- combined_data %>%
  mutate(
    # Extract the leading datetime block: "20250513_151000"
    datetime_string = str_extract(source_file, "^\\d{8}_\\d{6}"),
    
    # Convert to a clean format: "2025-05-13 15:10:00"
    Datetime = ymd_hms(str_replace(datetime_string, "_", " "))
  )

combined_data <- combined_data %>%
  mutate(
    Recorder_ID = str_extract(File, "[A-Z]\\d{3}")
  )

combined_data <- combined_data %>%
  mutate(
    # Extract folder name containing the habitat (RecorderID_Habitat)
    Folder = str_extract(File, "[A-Z]\\d{3}_[A-Za-z0-9_]+"),
    
    # Remove the recorder ID and underscore → get only habitat
    Habitat_raw = str_replace(Folder, "^[A-Z]\\d{3}_", ""),
    
    # Replace underscores with spaces
    Habitat_clean = str_replace_all(Habitat_raw, "_", " "),
    
    # Apply override rule for J001
    Habitat = ifelse(Recorder_ID == "J001", "Maple_Beech", Habitat_clean)
  )

combined_data <- combined_data %>%
  select(-Folder, -Habitat_raw, -Habitat_clean)

combined_data <- combined_data %>%
  mutate(
    # Extract the YYYYMMDD_HHMMSS part using regex
    datetime_str = sub("(\\d{8}_\\d{6}).*", "\\1", source_file),
    # Parse it into POSIXct datetime
    Datetime = ymd_hms(gsub("_", " ", datetime_str))
  ) %>%
  select(-datetime_str)  # Remove helper column

# Check result
head(combined_data)

write.csv(combined_data, "combined_BirdNET_data.csv")

#### Part 2. Rarefraction ####

combined_data <- read.csv("combined_BirdNET_data.csv")

# ------------------------------------------------------------
# Load libraries
# ------------------------------------------------------------
library(dplyr)
library(vegan)

# ------------------------------------------------------------
# 1. Define each audio file as a "sample"
# ------------------------------------------------------------
BirdNET_Data <- BirdNET_Data %>%
  mutate(Sample = File)   # each WAV file = one sampling unit

# ------------------------------------------------------------
# 2. List unique recorder IDs
# ------------------------------------------------------------
recorders <- unique(BirdNET_Data$Recorder_ID)

# ------------------------------------------------------------
# 3. Build species × sample matrices for each recorder
# ------------------------------------------------------------
spec_matrices <- list()

for (rec in recorders) {
  
  df <- BirdNET_Data %>%
    filter(Recorder_ID == rec)
  
  # Create species-by-sample abundance matrix
  mat <- table(df$Sample, df$Common.name)
  
  spec_matrices[[rec]] <- mat
}

# ------------------------------------------------------------
# 4. Plot rarefaction curves for each recorder ID
# ------------------------------------------------------------
# Adjust the panel layout depending on how many recorders you have
n <- length(recorders)
rows <- ceiling(n / 3)

par(mfrow = c(rows, 3))   # 3 plots per row

for (rec in recorders) {
  
  mat <- spec_matrices[[rec]]
  
  rarecurve(mat,
            step = 20,
            col = "darkgreen",
            label = FALSE,
            main = paste("Recorder", rec),
            xlab = "Detections (N)",
            ylab = "Species Richness")
}


accum_list <- map2(spec_matrices, names(spec_matrices), function(mat, rec) {
  acc <- specaccum(mat, method = "random")
  tibble(
    Recorder_ID = rec,
    Samples = acc$sites,
    Richness = acc$richness
  )
})

accum_df <- bind_rows(accum_list)

Rarefaction <- ggplot(accum_df, aes(x = Samples, y = Richness, color = Recorder_ID)) +
  geom_line(size = 1) +
  theme_bw() +
  labs(
    title = "",
    x = "Number of audio files (Samples)",
    y = "Species richness"
  )

ggsave("BirdNET Rarefaction.jpeg", plot = Rarefaction, device = "jpeg",
       width = 5, height = 5, units = "in", dpi = 300)

```



