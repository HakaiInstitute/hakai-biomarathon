# Install requirements for phyloseq
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")

# Functions
get_worms_lineage <- function(tax) {
  result <- tryCatch(classification(tax, db = "worms"), error = function(e) {
    NULL
  })
  worms_data <- result[[1]]

  if (!is.data.frame(worms_data)) {
    return(tibble(Species = tax)) # fallback: just return species name
  }

  # convert worms_data into wide form
  lineage <- as_tibble(t(worms_data$name), .name_repair = "minimal")
  colnames(lineage) <- worms_data$rank
  lineage$Species <- tax # retain for joining

  # get id number from worms, format it for AphiaID/LSID
  lineage$AphiaID <- paste0(
    "urn:lsid:marinespecies.org:taxname:",
    last(worms_data$id)
  )

  # get the rank of the lowest taxon
  lineage$Rank <- last(worms_data$rank)

  lineage
}

get_authority <- function(id) {
  tryCatch(worrms::wm_record(id = id)$authority, error = function(e) NA)
}

# Load libraries
library(tidyverse)
library(taxize)
library(data.table)
library(phyloseq)
library(here)

# Placeholder labels to treat as missing, and the identification qualifiers
# (cf./aff.) that should be pulled out of Number_species and dropped from
# scientificName
na_labels <- c(
  "unknown",
  "no identification",
  "unknown class",
  "unknown order",
  "unknown phylum",
  "unknown family",
  "unknown genus",
  "uncultured *"
)
qualifier_pattern <- "cf. promare|cf. depressum|cf. urtica|aff. granulata"

# Read in taxonomic data, filter for Eukaryota and Bacteria
taxtab_COI <- fread(
  here("data", "Biomarathon_taxonomy_file.csv"),
  sep = ",",
  header = TRUE,
  data.table = FALSE
) |>
  janitor::clean_names() |>
  select(1:number_species) |>
  filter(number_kingdom == c("Eukaryota", "Bacteria"))

colnames(taxtab_COI)[-1] <- str_to_title(colnames(taxtab_COI)[-1])
taxtab_COI <- taxtab_COI |> column_to_rownames("number_query")

# Convert to taxonomyTable class
tax_tab_COI <- as.data.frame(tax_table(as.matrix(taxtab_COI)))

# Make uncultured to unknown, then recode all placeholder labels as NA
tax_tab_COI <- tax_tab_COI |>
  mutate(
    Number_species = gsub("^uncultured.*", "unknown", Number_species),
    across(everything(), ~ replace(.x, .x %in% na_labels, NA))
  )

# Grab the most granular taxonomic data (i.e., the last column populated before 'NA' is provided)
tax_tab_COI$scientificName <- apply(tax_tab_COI, 1, function(row) {
  last_value <- tail(na.omit(row), 1)
  if (length(last_value) == 0) NA else last_value
})

# Build an ASV-level lookup table: species name (cleaned of qualifiers/codes),
# verbatim identification, and identification qualifier
tax_lookup_tbl <- tax_tab_COI |>
  rownames_to_column("ASV") |>
  transmute(
    ASV,
    verbatimIdentification = Number_lca_taxon,
    identificationQualifier = str_extract(Number_species, qualifier_pattern),
    Species = scientificName
  ) |>
  mutate(
    Species = gsub(qualifier_pattern, "", Species),
    Species = gsub("sp\\..*", "", Species),
    Species = gsub("2 ML-2018|HIMB59|RCA23|A6200|CMC02|NN016038", "", Species),
    Species = trimws(Species)
  )

### Convert to WoRMS ---
unique_species_COI <- unique(tax_lookup_tbl$Species)

# Map over unique species names
## 1 696531 Nitzschia dissipata var antartica     Peragallo, 1921                        unassessed
## 1 149142 Navicula                              J.B.M. Bory de Saint-Vincent, 1822     accepted
## 1 149045 Nitzschia                             A.H. Hassall, 1845                     accepted
## 1 110993 Alcyonidium                           Lamouroux, 1813                        accepted
## 1 120564 Hyalodiscus                           Hertwig & Lesser, 1874                 accepted
## 1 129619 Polydora                              Bosc, 1802                             accepted
## 1 345835 Spumellaria                           Ehrenberg, 1875                        accepted
## 1 415427 Glyptaesopus                          Pilsbry & Olsson, 1941                 accepted
## 2 137357 Nais                                  Muller, 1774                           accepted
## 2 129625 Spio                                  Fabricius, 1785                        accepted
## 1 110041 Gonyaulax spinifera                   (Claparede & Lachmann) Diesing, 1866   accepted
## 2 130932 Myxicola infundibulum                 (Montagu, 1808)                        accepted
worms_taxonomy_tbl <- purrr::map_dfr(unique_species_COI, get_worms_lineage)

# Join lineage info back to ASVs, add authority, and finalize DwC columns
worms_tax_tab_COI_full <- tax_lookup_tbl |>
  left_join(worms_taxonomy_tbl, by = "Species") |>
  rename(scientificNameID = AphiaID, taxonRank = Rank) |>
  mutate(
    aphiaID = as.numeric(sub(
      "urn:lsid:marinespecies.org:taxname:",
      "",
      scientificNameID
    )),
    scientificNameAuthorship = sapply(aphiaID, get_authority)
  ) |>
  column_to_rownames("ASV") |>
  select(
    kingdom = Kingdom,
    phylum = Phylum,
    class = Class,
    order = Order,
    family = Family,
    genus = Genus,
    scientificName = Species,
    scientificNameID,
    taxonRank,
    verbatimIdentification,
    scientificNameAuthorship,
    identificationQualifier
  )

# Print species that don't have an LSID
worms_tax_tab_COI_full[is.na(worms_tax_tab_COI_full$scientificNameID), ]

# There are 141 taxa that do not have an associated scientificNameID
# or taxonomic hierarchical information. Verify with Libby and Matt
# that these can be omitted from OBIS submission.
no_aphiaid <- worms_tax_tab_COI_full[
  is.na(worms_tax_tab_COI_full$scientificNameID),
]

# Create final WoRMS-cleaned taxonomy table
worms_tax_tab_COI_full <- rownames_to_column(
  worms_tax_tab_COI_full,
  var = "ASV"
)

# Write to csv
write.csv(
  worms_tax_tab_COI_full,
  here("data", "Hakai_BioMarathon_taxonomy_COI_worms.csv")
)
