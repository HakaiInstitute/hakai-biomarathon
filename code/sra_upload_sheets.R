# create SRA metadata import file from lab metadata files
# based this off the SRA MIMARKS.survey.water.6.0.xls file 

library(here)
library(readr)
library(tidyverse)
library(janitor)

# Read in the metadata
batch_metadata <- readr::read_csv(
  here("data", "batch_sheet.csv"),
  na = c("", "NA")
) |>
  janitor::clean_names() |>
  pivot_wider(names_from = column_header, values_from = example_s_for_column)

lab_metadata <- readr::read_csv(
  here("data", "lab_sheet.csv"),
  guess_max = Inf
) |>
  janitor::clean_names() |>
  janitor::remove_empty("cols") |>
  mutate(tech = gsub("/", " | ", tech, fixed = TRUE))

field_metadata <- readr::read_csv(
  here("data", "field_sheet.csv"),
  guess_max = Inf
) |>
  janitor::clean_names() |>
  janitor::remove_empty("cols")

# set the project name to append to file names
# eg: project <- "Biomarathon"
project <- ""

# set the project accession id, you get this after you submit the bioproject data and add it to the batch sheet manually or you could set it here as a string
bioproject_accession <- batch_metadata$bioproject_accession

# Create table for the biosample data upload
sra_biosample <- left_join(
  field_metadata,
  lab_metadata,
  by = c("sample_id", "sample_type", "field_blank")
) |>
  filter(sample_type == "Seawater") |>
  dplyr::rename(
    sample_name = sample_id,
    depth = sample_depth_m) |>
  mutate(
    bioproject_accession = bioproject_accession,
    lat_lon = paste0(latitude, " N ", longitude, " W"),
    `Sample description` = sample_name
  ) |>
  select(sample_name, bioproject_accession, organism, collection_date, 
         depth, env_broad_scale, env_local_scale, env_medium, 
         geo_loc_name, lat_lon, `Sample description`) 

# write out
write_tsv(
  sra_biosample,
  here("data", paste0(project, "_SRA_biosample_upload.tsv")),
  na = ""
)



# after the biosample data are submitted we can create the sra metadata acc file
# download the attributes.tsv file from the SRA submission page

attributes <- readr::read_tsv(
  here("data", "attributes.tsv"),
  na = c("", "NA")
) |>
  janitor::clean_names()

# set title based on the bioproject title 
# eg: bioproj_title <- "Quadra Island BioMarathon Raw sequence reads"
bioproj_title <- ""

# get columns we'll select from batch metadata tibble
batch_cols <- c("library_strategy", "library_source", "library_selection",
                "library_layout", "platform", "instrument_model",
                "design_description")

# combine attributes data with necessary batch data
sra_metadata_acc <- attributes |>
  select(biosample_accession = accession, library_ID = sample_name) |>
  mutate(title = bioproj_title) |>
  cross_join(batch_metadata |> select(all_of(batch_cols))) |>
  mutate(filetype = "fastq") |>
  left_join(
    lab_metadata |>
      select(library_ID = sample_id,
             filename = forward_file_name,
             filename2 = reverse_file_name),
    by = "library_ID"
  )

# write out
write_tsv(
  sra_metadata_acc,
  here("data", paste0(project, "_SRA_metadata_acc_upload.tsv")),
  na = ""
)
