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
    # fallback: just return the species name; map_dfr() will fill every
    # other lineage column with NA for this row when it row-binds results
    return(tibble(Species = tax))
  }

  # convert worms_data into wide form
  lineage <- as_tibble(t(worms_data$name), .name_repair = "minimal")
  colnames(lineage) <- worms_data$rank
  lineage$Species <- tax # retain for joining

  # get id number from worms, format it for AphiaID/LSID
  aphia_id <- last(worms_data$id)
  lineage$AphiaID <- paste0("urn:lsid:marinespecies.org:taxname:", aphia_id)

  # get the rank of the lowest taxon
  lineage$Rank <- last(worms_data$rank)

  # look up authority here (once per unique species) instead of later
  # (once per ASV), since many ASVs share the same species/AphiaID
  lineage$scientificNameAuthorship <- get_authority(aphia_id)

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
qualifier_pattern <- "cf. promare|cf. depressum|cf. pugettensis|aff. obscurus|cf. curdieana|var. palustre"

# Descriptive (morphological) notes that follow "sp." in Number_species, which
# would otherwise be lost when "sp\\..*" is stripped from scientificName
remarks_pattern <- "pink proboscis"

# ASVs present in the ASV table; only these need taxonomic cleaning
asv_ids <- fread(
  here("data", "Biomarathon_ASV_table_MarineInverts.csv"),
  select = 1,
  data.table = FALSE
)[[1]]

# Read in taxonomic data, filter for ASVs in the ASV table that are
# Eukaryota or Bacteria
taxtab_COI <- fread(
  here("data", "Biomarathon_taxonomy_file.csv"),
  sep = ",",
  header = TRUE,
  data.table = FALSE
) |>
  janitor::clean_names() |>
  select(1:number_species) |>
  filter(
    number_query %in% asv_ids
  )

colnames(taxtab_COI)[-1] <- str_to_title(colnames(taxtab_COI)[-1])
taxtab_COI <- taxtab_COI |> column_to_rownames("number_query")

# Route through phyloseq's tax_table() purely to force a uniform character
# matrix (rownames = ASV)
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
    identificationRemarks = str_extract(Number_species, remarks_pattern),
    # capture the BOLD process/BIN ID (e.g. "BOLD:AAO6949") before it gets
    # stripped out of Species below by the "sp\\..*" cleanup
    taxonConceptID = str_extract(scientificName, "BOLD:[A-Z0-9]+"),
    Species = scientificName
  ) |>
  mutate(
    Species = gsub(qualifier_pattern, "", Species),
    Species = gsub("sp\\..*", "", Species),
    # Strip non-taxonomic placeholder/clade codes (e.g., candidate lineage
    # names such as HIMB59, RCA23) that the classifier reports in place of
    # an actual species epithet
    Species = gsub(
      "2 ML-2018|HIMB59|RCA23|A6200|CMC02|NN016038|FS 5|HTCC1062|IMCC1322|1 AM-2017|CMC01|1 ML-2018",
      "",
      Species
    ),
    Species = trimws(Species)
  )

# There are a few taxa names that we have to change, based on conversations with Libby Natola
# We update those names below before getting their lineage below:
tax_lookup_tbl$Species <- gsub(
  "Trieres chinensis",
  "Odontella chinensis",
  tax_lookup_tbl$Species
)
tax_lookup_tbl$Species <- gsub(
  "Saprolegnia parasitica",
  "Saprolegnia",
  tax_lookup_tbl$Species
)
tax_lookup_tbl$Species <- gsub(
  "Clydonella sawyeri",
  "Clydonella",
  tax_lookup_tbl$Species
)
tax_lookup_tbl$Species <- gsub(
  "Nitzschia inconspicua",
  "Nitschia frustulum var. inconspicua",
  tax_lookup_tbl$Species
)
tax_lookup_tbl$Species <- gsub(
  "Ralfsioid",
  "Ralfsiaceae",
  tax_lookup_tbl$Species
)
tax_lookup_tbl$Species <- gsub(
  "Chiridota albatrossii complex",
  "Chiridota albatrossii",
  tax_lookup_tbl$Species
)
tax_lookup_tbl$Species <- gsub(
  "Neoporphyra perforata",
  "Pyropia perforata",
  tax_lookup_tbl$Species
)
tax_lookup_tbl$Species <- gsub(
  "Schizostauron trachyderma",
  "Achnanthes trachyderma",
  tax_lookup_tbl$Species
)
tax_lookup_tbl$Species <- gsub(
  "Picobiliphyte",
  "Picozoa",
  tax_lookup_tbl$Species
)
tax_lookup_tbl$Species <- gsub(
  "Chondria cf. curdieana",
  "Chondria",
  tax_lookup_tbl$Species
)
tax_lookup_tbl$Species <- gsub(
  "Candidatus Arcticimaribacter forsetii",
  "Flavobacteriaceae",
  tax_lookup_tbl$Species
)
tax_lookup_tbl$Species <- gsub(
  "[Olpidiopsis] porphyrae var. scotiae (nom. inval.)",
  "Olpidiopsis",
  tax_lookup_tbl$Species,
  fixed = TRUE
)
tax_lookup_tbl$Species <- gsub(
  "Olleya aquimaris",
  "Olleya",
  tax_lookup_tbl$Species
)

### Create list of unique taxa
unique_species_COI <- unique(tax_lookup_tbl$Species)

# Map over unique species names
worms_taxonomy_tbl <- purrr::map_dfr(unique_species_COI, get_worms_lineage)

## 1 129619 Polydora                              Bosc, 1802                                                                accepted
## 1 110993 Alcyonidium                           Lamouroux, 1813                                                           accepted
## 2 333595 Hermadion truncata                    (Czerniavsky, 1882)                                                       nomen dubium
## 2 129296 Glycera                               Lamarck, 1818                                                             accepted
## 1 115379 Paralaophonte                         Lang, 1948                                                                accepted
## 2 138019 Flabellina                            McMurtrie, 1831                                                           accepted
## 2 129625 Spio                                  Fabricius, 1785                                                           accepted
## 2 130932 Myxicola infundibulum                 (Montagu, 1808)                                                           accepted
## 2 338908 Lepidasthenia longicirrata            Berkeley, 1923                                                            accepted

# Join lineage info (including authorship, looked up per unique species in
# get_worms_lineage) back to ASVs, and finalize DwC columns
worms_tax_tab_COI_full <- tax_lookup_tbl |>
  left_join(worms_taxonomy_tbl, by = "Species") |>
  rename(scientificNameID = AphiaID, taxonRank = Rank) |>
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
    identificationQualifier,
    identificationRemarks,
    taxonConceptID
  )

# Join original taxonomy table with the WoRMS cleaned table to
# do a manual review, joining by ASV
worms_tax_tab_COI_full <- worms_tax_tab_COI_full |> rownames_to_column("ASV")
tax_tab_COI <- tax_tab_COI |> rownames_to_column("ASV")

# Where a taxa does not have any associated WoRMS taxonomic data, but that have been
# identified as (marine) bacteria (by LN) include "incertae sedis" as scientificName.
incertae_sedis_taxa <- c(
  "alpha proteobacterium HIMB59",
  "Planktomarina temperata RCA23",
  "Urechidicola croceus",
  "Polaribacter reichenbachii",
  "Serratia",
  "uncultured bacterium TAYNAYA-8",
  "Myroides odoratimimus",
  "Methylophilales bacterium",
  "Ralstonia insidiosa",
  "Candidatus Pelagibacter giovannonii",
  "Candidatus Nucleicultrix amoebiphila FS5"
)

worms_tax_tab_COI_full <- worms_tax_tab_COI_full |>
  mutate(
    scientificName = if_else(
      verbatimIdentification %in% incertae_sedis_taxa,
      "incertae sedis",
      scientificName
    ),
    scientificNameID = if_else(
      verbatimIdentification %in% incertae_sedis_taxa,
      "urn:lsid:marinespecies.org:taxname:12",
      scientificNameID
    )
  )

# Write to csv
write.csv(
  worms_tax_tab_COI_full,
  here("data", "Hakai_BioMarathon_taxonomy_COI_worms.csv"),
  row.names = FALSE
)
