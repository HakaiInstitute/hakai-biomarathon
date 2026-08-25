# Install packages
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}
pak::pak("iobis/obistools")
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
  janitor::clean_names() |>
  pivot_wider(names_from = column_header, values_from = example_s_for_column)
lab_metadata <- readxl::read_excel(
  here("data", "Lab Metadata Template - Aquatic.xlsx"),
  sheet = "lab sheet (draft)"
) |>
  janitor::clean_names() |>
  janitor::remove_empty("cols") |>
  mutate(tech = gsub("/", " | ", tech, fixed = TRUE))
field_metadata <- readxl::read_excel(
  here("data", "Lab Metadata Template - Aquatic.xlsx"),
  sheet = "field sheet (draft)"
) |>
  janitor::clean_names() |>
  janitor::remove_empty("cols")

# Read in the data
asv <- readr::read_csv(here("data", "Biomarathon_ASV_table_MarineInverts.csv"))
taxa <- readr::read_csv(here("data", "Biomarathon_taxonomy_file.csv"))
fasta <- seqinr::read.fasta(
  file = here("data", "CO1.biomarathon_combined.ASV_sequences.fasta"),
  as.string = TRUE,
  forceDNAtolower = FALSE
)

# Create table for sampling event metadata
meta <- left_join(
  field_metadata,
  lab_metadata,
  by = c("sample_id", "sample_type")
) |>
  filter(sample_type == "Seawater") |>
  dplyr::rename(
    materialSampleID = sample_id,
    verbatimLocality = site,
    decimalLatitude = latitude,
    decimalLongitude = longitude,
    fieldNotes = field_notes,
    samplingProtocol = collection_method,
    country = geo_loc_name,
    recordedBy = tech
  ) |>
  mutate(
    datasetID = "Hakai-Biomarathon",
    language = "en",
    license = "http://creativecommons.org/licenses/by/4.0/legalcode",
    rightsHolder = "Hakai Institute",
    modified = lubridate::today(),
    geodeticDatum = "WGS84",
    stateProvince = "British Columbia",
    institutionCode = "Hakai Institute",
    institutionID = "https://edmo.seadatanet.org/report/5148",
    eventID = paste("Hakai-Biomarathon", materialSampleID, sep = "-"),
    eventDate = biocleanr::create_event_date(collection_date),
    year = as.numeric(lubridate::year(eventDate)),
    month = as.numeric(lubridate::month(eventDate)),
    day = as.numeric(lubridate::day(eventDate)),
    countryCode = "CA",
    minimumDepthInMeters = sample_depth_m,
    maximumDepthInMeters = sample_depth_m
  )

## Create Occurrence Core

# Read in the cleaned worms taxonomy table
taxonomy <- read_csv(
  here(
    "data",
    "Hakai_BioMarathon_taxonomy_COI_worms.csv"
  )
)

# Pivot the ASV table
asv <- asv |>
  dplyr::rename(ASV = row_names) |>
  select(ASV, starts_with("QU")) |>
  pivot_longer(
    starts_with("QU"),
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
  group_by(materialSampleID) |>
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

# Join ASV to samples and taxonomic information:
DNAtable <- dplyr::left_join(asv, taxonomy, by = "ASV")
DNAtable <- dplyr::inner_join(
  DNAtable,
  meta,
  by = c("materialSampleID", "eventID")
)
DNAtable <- dplyr::left_join(DNAtable, fasta_CO1, by = "ASV")

# Create occurrence core
# Columns shared between the occurrence-level table (HB_occ) and the final
# core (HB_occ_core, which adds the event-level fields joined from meta)
occurrence_fields <- c(
  "occurrenceID",
  "basisOfRecord",
  "scientificName",
  "scientificNameID",
  "scientificNameAuthorship",
  "kingdom",
  "phylum",
  "class",
  "order",
  "family",
  "genus",
  "taxonRank",
  "verbatimIdentification",
  "identificationQualifier",
  "organismQuantity",
  "organismQuantityType",
  "sampleSizeValue",
  "sampleSizeUnit",
  "materialSampleID"
)
event_fields <- c(
  "eventID",
  "datasetID",
  "language",
  "license",
  "rightsHolder",
  "modified",
  "geodeticDatum",
  "stateProvince",
  "institutionCode",
  "institutionID",
  "verbatimLocality",
  "decimalLatitude",
  "decimalLongitude",
  "fieldNotes",
  "minimumDepthInMeters",
  "maximumDepthInMeters",
  "samplingProtocol",
  "country",
  "eventDate",
  "year",
  "month",
  "day",
  "countryCode"
)

HB_occ <- DNAtable |>
  select(ASV, eventID, all_of(occurrence_fields))

# Make sure all NA cells are blank:
HB_occ[is.na(HB_occ)] <- ""

# Create core:
HB_occ_core <- HB_occ |>
  left_join(meta, by = c("eventID", "materialSampleID")) |>
  select(all_of(occurrence_fields), all_of(event_fields)) |>
  mutate(
    occurrenceStatus = "present",
    taxonRank = tolower(taxonRank)
  )

# Validate the Occurrence table
kit <- biocleanr::dwc_init(
  "hakai-biomarathon",
  description = "Hakai Bio Marathon data",
  standard = "obis"
) |>
  biocleanr::dwc_dataset(
    "Occurrence",
    HB_occ_core,
    name = "HB_occ_core"
  )

report <- biocleanr::dwc_validate(kit)
report

write_csv(
  HB_occ_core,
  here("obis", "Hakai_BM_occ.csv"),
  na = ""
)

## Create the DNA-derived data extension
# DNAtable already has meta joined in (see above), so we only need to add
# the batch-level (run-wide) metadata here via a cross join, since
# batch_metadata is a single row of constants that applies to every record
sequence_meta <- DNAtable |>
  dplyr::cross_join(batch_metadata) |>
  mutate(
    seq_meth = paste(platform, instrument_model, sep = " "),
    otu_class_appr = "",
    otu_seq_comp_appr = "",
    otu_db = ""
  ) |>
  select(
    occurrenceID,
    DNA_sequence,
    samp_collec_method = samplingProtocol,
    sop = SOP_link,
    target_gene,
    target_subfragment = target_subragment,
    samp_size = volume_filtered_ml,
    pcr_primer_forward,
    pcr_primer_name_forward,
    pcr_primer_reverse,
    pcr_primer_name_reverse,
    env_broad_scale,
    env_local_scale,
    seq_meth,
    lib_layout = library_layout,
    otu_class_appr,
    otu_seq_comp_appr,
    otu_db
  )

# Clean it up - omit rows that are completely NA
# or where cells are NA make them blank
sequence_meta <- sequence_meta |>
  filter(!if_all(everything(), ~ is.na(.) | . == ""))

kit <- kit |>
  biocleanr::dwc_dataset(
    "dnaDerivedData",
    sequence_meta,
    name = "sequence_meta"
  )

report <- biocleanr::dwc_validate(kit)
report

# Save this data table in the obis folder
write_csv(
  sequence_meta,
  here("obis", "HB_DNA_data.csv"),
  na = ""
)
