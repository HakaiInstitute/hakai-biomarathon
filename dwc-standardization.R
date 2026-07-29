# Install packages
install.packages("devtools")
devtools::install_github("iobis/obistools", force = TRUE)
pak::pak("HakaiInstitute/biocleanr")

# Load packages
library(tidyr)
library(dplyr)
library(obistools)
library(janitor)
library(here)
library(readxl)
library(readr)
library(lubridate)
library(hms)

# Read in the metadata
batch_metadata <- readxl::read_excel(
  here("data", "Lab Metadata Template - Aquatic.xlsx"),
  sheet = "batch sheet (draft)"
) |>
  janitor::clean_names()
lab_metadata <- readxl::read_excel(
  here("data", "Lab Metadata Template - Aquatic.xlsx"),
  sheet = "lab sheet (draft)"
) |>
  janitor::clean_names()
field_metadata <- readxl::read_excel(
  here("data", "Lab Metadata Template - Aquatic.xlsx"),
  sheet = "field sheet (draft)"
) |>
  janitor::clean_names()

# Read in the data
asv <- readr::read_csv(here("data", "Biomarathon_ASV_table_MarineInverts.csv"))
taxa <- readr::read_csv(here("data", "Biomarathon_taxonomy_file.csv"))
fasta <- seqinr::read.fasta(
  file = here("data", "CO1.biomarathon_combined.ASV_sequences.fasta"),
  as.string = TRUE,
  forceDNAtolower = FALSE
)

# Create table for sampling event metadata
biom <- data.frame(
  eventID = "Hakai-Biomarathon",
  language = "en",
  license = "http://creativecommons.org/licenses/by/4.0/legalcode",
  rightsHolder = "Hakai Institute",
  modified = lubridate::today(),
  geodeticDatum = "WGS84",
  stateProvince = "British Columbia",
  institutionCode = "Hakai Institute",
  institutionID = "https://edmo.seadatanet.org/report/5148"
)

meta <- left_join(
  field_metadata,
  lab_metadata,
  by = c("sample_id", "sample_type")
) |>
  dplyr::rename(
    materialSampleID = sample_id,
    verbatimLocality = site,
    decimalLatitude = latitude,
    decimalLongitude = longitude,
    fieldNotes = field_notes,
    minimumDepthInMeters = sample_depth_m,
    maximumDepthInMeters = sample_depth_m,
    samplingProtocol = collection_method,
    country = geo_loc_name
  ) |>
  mutate(
    eventDate = biocleanr::create_event_date(collection_date),
    year = as.numeric(lubridate::year(eventDate)),
    month = as.numeric(lubridate::month(eventDate)),
    day = as.numeric(lubridate::day(eventDate)),
    countryCode = "CA"
  )

## Create Occurrence Core

# Read in the cleaned worms taxonomy table

# Pivot the ASV table
asv <- asv |>
  dplyr::rename(ASV = row_names) |>
  select(ASV, QU380001:QU380500) |>
  pivot_longer(
    QU380001:QU380500,
    names_to = "materialSampleID",
    values_to = "organismQuantity"
  ) |>
  filter(organismQuantity != 0) |>
  mutate(
    organismQuantityType = "DNA sequence reads",
    sampleSizeUnit = "DNA sequence reads",
    basisOfRecord = "MaterialSample",
    eventID = paste("Hakai-Biomarathon", materialSampleID, sep = "-")
  ) |>
  group_by(materialSampleID, ASV) |>
  mutate(
    occurrenceID = paste(eventID, ASV, "occ", sep = "-"),
    sampleSizeValue = sum(organismQuantity)
  ) |>
  ungroup()

# Match ASV to the sequence and taxa information:
asv_CO1 <- names(fasta)
DNA_sequence_CO1 <- unname(unlist(fasta))
fasta_CO1 <- cbind(asv_CO1, DNA_sequence_CO1) |>
  as.data.frame() |>
  dplyr::rename(ASV = asv_CO1, DNA_sequence = DNA_sequence_CO1)
