# Code

The scripts in this folder turn the Hakai BioMarathon CO1 eDNA data into Darwin Core tables for publication to OBIS. Run them in numbered order, because `02` reads the cleaned taxonomy that `01` produces.

```
01_worms_taxonomic_cleaning.R
  in:  data/Biomarathon_taxonomy_file.csv
       data/Biomarathon_ASV_table_MarineInverts.csv
  out: data/Hakai_BioMarathon_taxonomy_COI_worms.csv  ──┐
                                                         │
02_dwc_standardization.R                                 │
  in:  data/Hakai_BioMarathon_taxonomy_COI_worms.csv  ◄──┘
       data/Biomarathon_ASV_table_MarineInverts.csv
       data/Biomarathon_taxonomy_file.csv
       data/CO1.biomarathon_combined.ASV_sequences.fasta
       data/field_sheet.csv, lab_sheet.csv, batch_sheet.csv
  out: obis/Hakai_BM_occ.csv            (Occurrence core)
       obis/HB_DNA_data.csv             (DNA-derived data extension)
       data/Biomarathon_ASV_table_obis_comparison.csv
```

## 01_worms_taxonomic_cleaning.R

Cleans the ASV taxonomy and matches it to the World Register of Marine Species (WoRMS).

- Reads the classifier output (`data/Biomarathon_taxonomy_file.csv`) and keeps only the ASVs that are in `data/Biomarathon_ASV_table_MarineInverts.csv`.
- Takes the most specific rank available for each ASV as its `scientificName`.
- Extracts identification qualifiers (cf./aff.), identification remarks and BOLD BIN IDs, and strips non-taxonomic clade or strain codes from the names.
- Replaces names that WoRMS doesn't recognize, based on review by a Subject Matter Expert.
- Queries WoRMS once for each unique name to get the full lineage, AphiaID (`scientificNameID`), rank and authorship.
- Writes `data/Hakai_BioMarathon_taxonomy_COI_worms.csv`, one row per ASV.

## 02_dwc_standardization.R

Builds the Darwin Core Occurrence core and the DNA-derived data extension.

- Reads the sample metadata from `data/field_sheet.csv`, `data/lab_sheet.csv` and `data/batch_sheet.csv`, and combines the field and lab sheets into event-level metadata for the seawater samples (location, date, depth, protocol, etc.).
- Converts the ASV read-count table to long format, dropping zero counts, and joins each ASV to its WoRMS taxonomy (from script `01`), its sample metadata and its DNA sequence from the FASTA file.
- Builds the Occurrence core, validates it with `biocleanr::dwc_validate()` and writes `obis/Hakai_BM_occ.csv`.
- Adds the run-wide batch metadata (primers, target gene, sequencing platform, etc.) to create the DNA-derived data extension, validates it and writes `obis/HB_DNA_data.csv`.
- Writes `data/Biomarathon_ASV_table_obis_comparison.csv`, which lists every ASV in the ASV table and whether it was included in the occurrence core.
