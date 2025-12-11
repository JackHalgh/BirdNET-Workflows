
```
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

#### Part 2. Rarefaction ####

combined_data <- read.csv("combined_BirdNET_data.csv")


# Load libraries

library(dplyr)
library(vegan)


# 1. Define each audio file as a "sample"

BirdNET_Data <- BirdNET_Data %>%
  mutate(Sample = File)   # each WAV file = one sampling unit


# 2. List unique recorder IDs

recorders <- unique(BirdNET_Data$Recorder_ID)


# 3. Build species × sample matrices for each recorder

spec_matrices <- list()

for (rec in recorders) {
  
  df <- BirdNET_Data %>%
    filter(Recorder_ID == rec)
  
  # Create species-by-sample abundance matrix
  mat <- table(df$Sample, df$Common.name)
  
  spec_matrices[[rec]] <- mat
}


# 4. Plot rarefaction curves for each recorder ID

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

#### Part 3. Heatmaps ####

library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(viridis)

# Load data
combined_data <- read.csv("combined_BirdNET_data.csv")

head(combined_data)

# Ensure Datetime is POSIXct and extract Date and Hour
combined_data <- combined_data %>%
  mutate(
    Datetime = parse_date_time(Datetime, orders = "Ymd_HMS"),  # parse YYYYMMDD_HHMMSS
    Date = as.Date(Datetime),
    Hour = hour(Datetime)
  )

# Call count

# Aggregate call count per Habitat, Date, and Hour
call_count_data <- combined_data %>%
  filter(!is.na(Hour)) %>%               # Remove rows with NA Hour
  group_by(Habitat, Date, Hour) %>%
  summarise(CallCount = n(), .groups = "drop") %>%  # Count calls
  mutate(HourLabel = sprintf("%02d:00", Hour))

# Create complete grid to fill missing combinations with 0
call_count_complete <- call_count_data %>%
  complete(Habitat, Date, Hour, fill = list(CallCount = 0)) %>%
  mutate(HourLabel = sprintf("%02d:00", Hour))

# Plot heatmap of call counts
p_callcount <- ggplot(call_count_complete, 
                      aes(y = Date, x = factor(HourLabel, levels = sprintf("%02d:00", 0:23)), fill = CallCount)) +
  geom_tile(color = NA) +  # remove tile borders
  facet_wrap(~ Habitat, scales = "free_y") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
  scale_x_discrete(breaks = sprintf("%02d:00", seq(0, 23, by = 2))) +  # show every 2 hours
  labs(y = "Date", x = "Hour of day", fill = "Call count") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"))

print(p_callcount)

ggsave("call_count_heatmap.jpeg", plot = p_callcount, device = "jpeg",
       width = 10, height = 10, units = "in", dpi = 300)

# Species richness

# Aggregate species richness per Habitat, Date, and Hour
richness_data <- combined_data %>%
  filter(!is.na(Hour)) %>%                     # remove rows with NA Hour
  group_by(Habitat, Date, Hour) %>%
  summarise(SpeciesRichness = n_distinct(Scientific.name), .groups = "drop") %>%  # unique species count
  mutate(HourLabel = sprintf("%02d:00", Hour))

# Create complete grid to fill missing combinations with 0
richness_complete <- richness_data %>%
  complete(Habitat, Date, Hour, fill = list(SpeciesRichness = 0)) %>%
  mutate(HourLabel = sprintf("%02d:00", Hour))

# Plot heatmap of species richness
p_richness <- ggplot(richness_complete, 
                     aes(y = Date, x = factor(HourLabel, levels = sprintf("%02d:00", 0:23)), fill = SpeciesRichness)) +
  geom_tile(color = NA) +  # remove tile borders
  facet_wrap(~ Habitat, scales = "free_y") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
  scale_x_discrete(breaks = sprintf("%02d:00", seq(0, 23, by = 2))) +  # show every 2 hours
  labs(y = "Date", x = "Hour of day", fill = "Species richness") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"))

print(p_richness)

ggsave("species_richness_heatmap.jpeg", plot = p_richness, device = "jpeg",
       width = 10, height = 10, units = "in", dpi = 300)

#### Part 4. Density time series plots ####

combined_data <- read.csv("combined_BirdNET_data.csv")
head(combined_data)

# Load required libraries
library(dplyr)
library(ggplot2)
library(lubridate)

# If your datetime column is not POSIXct, convert it
combined_data$Datetime <- as.POSIXct(combined_data$Datetime, format="%Y-%m-%d %H:%M:%S")

# 1. Aggregate data by time intervals (e.g., daily)
# You can also do hourly: floor_date(Datetime, "hour")
species_time <- combined_data %>%
  mutate(Date = as.Date(Datetime)) %>%  # Change to floor_date(Datetime, "hour") for hourly
  group_by(Date, Scientific.name) %>%
  summarise(Detections = n(), .groups = "drop")

# 2. Calculate relative composition per time interval
species_rel <- species_time %>%
  group_by(Date) %>%
  mutate(Relative = Detections / sum(Detections)) %>%
  ungroup()

# 3. Plot as stacked area plot
ggplot(species_rel, aes(x = Date, y = Relative, fill = Scientific.name)) +
  geom_area(alpha = 0.8 , color = "black", size = 0.2) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Relative Bird Species Composition Over Time",
    x = "Date",
    y = "Relative Composition (%)",
    fill = "Species"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )



# 3. Plot as stacked area plot without legend
ggplot(species_rel, aes(x = Date, y = Relative, fill = Scientific.name)) +
  geom_area(alpha = 0.8 , color = "black", size = 0.2) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "",
    x = "Date",
    y = "Relative Composition (%)",
    fill = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"  # removes the legend
  )

# top 10

library(dplyr)
library(lubridate)

# 1. Add a week column
species_rel <- species_rel %>%
  mutate(Week = floor_date(Date, "week"))

# 2. Aggregate relative composition per common name per week
weekly_species <- combined_data %>%
  mutate(Date = as.Date(Datetime),
         Week = floor_date(Date, "week")) %>%
  group_by(Week, Common.name) %>%
  summarise(WeeklyDetections = n(), .groups = "drop") %>%
  group_by(Week) %>%
  mutate(Relative = WeeklyDetections / sum(WeeklyDetections)) %>%  # relative contribution per week
  ungroup()

# 3. Calculate total relative contribution per common name across all weeks
total_contributions <- weekly_species %>%
  group_by(Common.name) %>%
  summarise(TotalRelative = sum(Relative), .groups = "drop") %>%
  arrange(desc(TotalRelative))

# 4. Select top 10 species by common name
top10_species <- total_contributions %>%
  slice_head(n = 10)

# 5. Print top 10 species with their total relative contributions
print(top10_species)

# Optional: show weekly breakdown for only top 10 common names
weekly_top10 <- weekly_species %>%
  filter(Common.name %in% top10_species$Common.name) %>%
  arrange(Week, desc(Relative))

print(weekly_top10)
head(weekly_top10)

library(ggplot2)
library(dplyr)
library(scales)

# Make sure Common.name is a factor ordered by total contribution
top_species_order <- top10_species$Common.name

weekly_top10 <- weekly_top10 %>%
  mutate(Common.name = factor(Common.name, levels = top_species_order))

# Plot weekly relative composition as a stacked area plot
ggplot(weekly_top10, aes(x = Week, y = Relative, fill = Common.name)) +
  geom_area(alpha = 0.8, color = "black", size = 0.2) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "",
    x = "Week",
    y = "Relative Composition (%)",
    fill = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

# top 20

library(dplyr)
library(lubridate)
library(ggplot2)
library(scales)

# 1. Add a week column
species_rel <- species_rel %>%
  mutate(Week = floor_date(Date, "week"))

# 2. Aggregate relative composition per common name per week
weekly_species <- combined_data %>%
  mutate(Date = as.Date(Datetime),
         Week = floor_date(Date, "week")) %>%
  group_by(Week, Common.name) %>%
  summarise(WeeklyDetections = n(), .groups = "drop") %>%
  group_by(Week) %>%
  mutate(Relative = WeeklyDetections / sum(WeeklyDetections)) %>%  # relative contribution per week
  ungroup()

# 3. Calculate total relative contribution per common name across all weeks
total_contributions <- weekly_species %>%
  group_by(Common.name) %>%
  summarise(TotalRelative = sum(Relative), .groups = "drop") %>%
  arrange(desc(TotalRelative))

# 4. Select top 20 species by common name
top20_species <- total_contributions %>%
  slice_head(n = 20)

# 5. Optional: show weekly breakdown for only top 20 common names
weekly_top20 <- weekly_species %>%
  filter(Common.name %in% top20_species$Common.name) %>%
  arrange(Week, desc(Relative))

# Make sure Common.name is a factor ordered by total contribution
weekly_top20 <- weekly_top20 %>%
  mutate(Common.name = factor(Common.name, levels = top20_species$Common.name))

# Plot weekly relative composition as a stacked area plot
ggplot(weekly_top20, aes(x = Week, y = Relative, fill = Common.name)) +
  geom_area(alpha = 0.8, color = "black", size = 0.2) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "",
    x = "Week",
    y = "Relative Composition (%)",
    fill = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

# top 30 #

library(dplyr)
library(lubridate)
library(ggplot2)
library(scales)

# 1. Add a week column
species_rel <- species_rel %>%
  mutate(Week = floor_date(Date, "week"))

# 2. Aggregate relative composition per common name per week
weekly_species <- combined_data %>%
  mutate(Date = as.Date(Datetime),
         Week = floor_date(Date, "week")) %>%
  group_by(Week, Common.name) %>%
  summarise(WeeklyDetections = n(), .groups = "drop") %>%
  group_by(Week) %>%
  mutate(Relative = WeeklyDetections / sum(WeeklyDetections)) %>%  # relative contribution per week
  ungroup()

# 3. Calculate total relative contribution per common name across all weeks
total_contributions <- weekly_species %>%
  group_by(Common.name) %>%
  summarise(TotalRelative = sum(Relative), .groups = "drop") %>%
  arrange(desc(TotalRelative))

# 4. Select top 30 species by common name
top30_species <- total_contributions %>%
  slice_head(n = 30)

# 5. Optional: show weekly breakdown for only top 30 common names
weekly_top30 <- weekly_species %>%
  filter(Common.name %in% top30_species$Common.name) %>%
  arrange(Week, desc(Relative))

# Make sure Common.name is a factor ordered by total contribution
weekly_top30 <- weekly_top30 %>%
  mutate(Common.name = factor(Common.name, levels = top30_species$Common.name))

# Plot weekly relative composition as a stacked area plot
density_plot <- ggplot(weekly_top30, aes(x = Week, y = Relative, fill = Common.name)) +
  geom_area(alpha = 0.8, color = "black", size = 0.2) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "",
    x = "Week",
    y = "Relative Composition (%)",
    fill = "Species"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave("density_plot.jpeg", plot = density_plot, width = 20, height = 10, units = "in", dpi = 300)

```



