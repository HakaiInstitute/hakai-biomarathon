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
library(readr)
library(lubridate)
library(hms)

# Read in the metadata. guess_max = Inf so that column types are guessed from
# all rows, not just the first 1000 (some columns are sparse)
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
  by = c("sample_id", "sample_type", "field_blank")
) |>
  filter(sample_type == "Seawater") |>
  dplyr::rename(
    materialSampleID = sample_id,
    verbatimLocality = site,
    decimalLatitude = latitude,
    decimalLongitude = longitude,
    fieldNotes = field_notes,
    samplingProtocol = collection_method,
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
    country = "Canada",
    countryCode = "CA",
    minimumDepthInMeters = sample_depth_m,
    maximumDepthInMeters = sample_depth_m
  )

## Create Occurrence Core

# Read in the cleaned worms taxonomy table. Read all columns as character:
# sparse columns (e.g. taxonConceptID) are empty in the first rows, so readr
# would otherwise guess them as logical and drop their values
taxonomy <- read_csv(
  here(
    "data",
    "Hakai_BioMarathon_taxonomy_COI_worms.csv"
  ),
  col_types = cols(.default = col_character())
)

# Pivot the ASV table, keeping only the seawater sample columns. Field blanks
# also have QU IDs, so select by the seawater IDs in meta rather than prefix
seawater_samples <- meta$materialSampleID

asv <- asv |>
  dplyr::rename(ASV = row_names) |>
  select(ASV, all_of(seawater_samples)) |>
  pivot_longer(
    all_of(seawater_samples),
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

# ASVs to leave out of the submission: freshwater taxa with no WoRMS match
# (ASV4870 Margaritifera falcata, ASV10266 Sulcospira) as per LN and ML
excluded_asvs <- c("ASV4870", "ASV10266")

# Join ASV to samples and taxonomic information:
DNAtable <- dplyr::left_join(asv, taxonomy, by = "ASV") |>
  filter(!ASV %in% excluded_asvs)
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
  "identificationRemarks",
  "taxonConceptID",
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
  "recordedBy",
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
    otu_class_appr = "100% identity; DADA2 denoising algorithm error model",
    otu_seq_comp_appr = "blastn;2.14.0;e-value cutoff: 0.00001",
    otu_db = "NCBI NT;255"
  ) |>
  select(
    occurrenceID,
    DNA_sequence,
    samp_collec_method = samplingProtocol,
    sop = SOP_link,
    target_gene,
    target_subfragment,
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
# and where cells are NA make them blank
sequence_meta <- sequence_meta |>
  filter(!if_all(everything(), ~ is.na(.) | . == ""))

kit <- kit |>
  biocleanr::dwc_dataset(
    "dnaDerivedData",
    sequence_meta,
    name = "DNA Derived"
  )

report <- biocleanr::dwc_validate(kit)
report

# Save this data table in the obis folder
write_csv(
  sequence_meta,
  here("obis", "HB_DNA_data.csv"),
  na = ""
)

## Compare the ASV table against the OBIS submission
# One row per ASV in the ASV table, showing which taxon it links to and
# whether it ended up in the occurrence core. Only seawater sample columns
# are used, so an ASV is not included if all its reads come from blank
# samples (zero seawater reads), or if it is listed in excluded_asvs
asv_ids <- readr::read_csv(
  here("data", "Biomarathon_ASV_table_MarineInverts.csv"),
  col_select = 1
) |>
  dplyr::rename(ASV = row_names)

submitted <- HB_occ |>
  distinct(ASV, scientificName, scientificNameID, identificationQualifier)

asv_obis_comparison <- asv_ids |>
  left_join(
    taxa |> select(ASV = `#Query`, verbatimIdentification = `#lca taxon`),
    by = "ASV"
  ) |>
  left_join(submitted, by = "ASV") |>
  left_join(taxonomy |> select(ASV, taxonConceptID), by = "ASV") |>
  mutate(includedInOBIS = ASV %in% submitted$ASV)

write_csv(
  asv_obis_comparison,
  here("data", "Biomarathon_ASV_table_obis_comparison.csv"),
  na = ""
)
