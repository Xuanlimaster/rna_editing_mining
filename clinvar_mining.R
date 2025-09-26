# Load necessary R packages
library(vcfR)
library(dplyr)
library(readr)

# Define ClinVar VCF file path
vcf_file <- "clinvar_20250923.vcf.gz"

# Read VCF file
vcf <- read.vcfR(vcf_file, verbose = FALSE)

# View basic information of the VCF file
# print(vcf)

# Examine INFO fields defined in the META section
# info_fields <- queryMETA(vcf, element = "INFO")
# print(info_fields)

# Convert fixed fields from VCF file to a data frame
fix_df <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)

# Extract and parse INFO field data into tidy format
vcf_info <- extract_info_tidy(vcf)

# Combine fixed fields
vcf_df <- cbind(fix_df, vcf_info) %>% select(-INFO)

# Write full VCF data to CSV
# write_csv(vcf_df, "clinvar_vcf_complete_parsed_data.csv")

# Filter variants theoretically correctable by ADAR (A-to-I) editing
adar_candidates <- vcf_df %>% filter(ALT == "A", REF == "G")

# Filter variants theoretically correctable by APOBEC (C-to-U) editing
apobec_candidates <- vcf_df %>% filter(ALT == "C", REF == "T")

# Merge both editing type candidates
all_candidates <- bind_rows(
  adar_candidates %>% mutate(RNA_Editing_Type = "ADAR_A_to_I"),
  apobec_candidates %>% mutate(RNA_Editing_Type = "APOBEC_C_to_U")
)

# Select final output columns to form comprehensive result table
candidates_table <- all_candidates %>% select(
  CHROM,            # Chromosome
  POS,              # Genomic position
  ID,               # Variant identifier
  REF,              # Reference allele
  ALT,              # Alternate allele
  RNA_Editing_Type, # Type of RNA editing mechanism applicable
  CLNSIG,           # Clinical significance
  CLNDN,            # Associated disease name
  CLNREVSTAT,       # ClinVar evidence level
  CLNDISDB,         # Disease database identifiers
  CLNHGVS           # HGVS nomenclature
)

# View first few rows of results
# print("Final candidate variants:")
# print(head(candidates_table, 10))

# Filter variants classified as Pathogenic
pathogenic_candidates <- candidates_table %>%
  filter(grepl("Pathogenic|Likely_pathogenic|conflicting", CLNSIG, ignore.case = TRUE))

# View the filtered pathogenic variants
# print(head(pathogenic_candidates, 10))

# Count the number of variants (#99356)
# print("Number of pathogenic/likely pathogenic candidates:")
# print(nrow(pathogenic_candidates))

# Define a reliability order for review status
# Ordered from most reliable to least reliable
review_status_order <- c("reviewed_by_expert_panel",
                         "practice_guideline",
                         "criteria_provided,_multiple_submitters,_no_conflicts",
                         "criteria_provided,_single_submitter")

# Convert review status to an ordered factor for sorting
pathogenic_candidates$CLNREVSTAT <- factor(pathogenic_candidates$CLNREVSTAT,
                                           levels = review_status_order,
                                           ordered = TRUE)

# Sort variants by review status reliability and clinical significance
prioritized_pathogenic <- pathogenic_candidates %>%
  arrange(CLNREVSTAT, CLNSIG)

# View sorted results and optionally save to file
# print(head(prioritized_pathogenic, 20))
# write_csv(prioritized_pathogenic,
#           "prioritized_pathogenic_variants_by_review_status.csv")

# Count and display diseases associated with pathogenic candidates
# disease_counts <- pathogenic_candidates %>% count(CLNDN, sort = TRUE)
# print("Top diseases associated with pathogenic RNA editing candidates:")
# print(head(disease_counts, n = 10))

# Load ClinVar variant summary data for additional annotations
variant_summary <- read_tsv("variant_summary.txt.gz")

# Examine variant summary data structure
# print("Variant Summary columns:")
# glimpse(variant_summary)

# Load allele-gene association data
allele_gene <- read_tsv("allele_gene.txt.gz", show_col_types = FALSE)

# Examine allele-gene data structure
# print("Allele-Gene columns:")
# glimpse(allele_gene)

# Integrate gene and clinical summaries with pathogenic variants data
pathogenic_gene_summary <- prioritized_pathogenic %>%
  left_join(
    allele_gene %>%
      group_by(ID = as.character(`#AlleleID`)) %>%
      summarise(
        
        # Concatenate unique GeneIDs
        GeneIDs = paste(unique(GeneID), collapse = ";"),
        
        # Concatenate unique gene symbols
        GeneSymbols = paste(unique(Symbol), collapse = ";"),
        
        # Count distinct genes affected
        Gene_Count = n_distinct(GeneID),
        .groups = "drop"
      ),
    by = "ID"
  ) %>%
  left_join(
    variant_summary %>%
      group_by(ID = as.character(`#AlleleID`)) %>%
      summarise(
        
        # All review statuses
        ReviewStatuses = paste(unique(ReviewStatus), collapse = "|"),
        
        # All associated phenotypes
        Phenotypes = paste(unique(PhenotypeList), collapse = "|"),
        
        # All origin types
        Origins = paste(unique(Origin), collapse = "|"),
        
        # GRCh38 availability flag
        Has_GRCh38 = any(Assembly == "GRCh38"),
        
        # GRCh37 availability flag
        Has_GRCh37 = any(Assembly == "GRCh37"),
        
        # GRCh38 chromosome
        GRCh38_Chromosome = first(Chromosome[Assembly == "GRCh38"]),
        
        # GRCh38 start position
        GRCh38_Start = first(Start[Assembly == "GRCh38"]),
        
        # GRCh38 stop position
        GRCh38_Stop = first(Stop[Assembly == "GRCh38"]),
        
        # Maximum submitter count
        Submitter_Count = max(NumberSubmitters, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "ID"
  ) %>%
  mutate(Is_MultiGene = Gene_Count > 1) # Flag variants affecting multiple genes

# Save pathogenic summary with comprehensive gene and clinical annotations
# write_csv(pathogenic_summary, 
#           "pathogenic_variants_full_gene_clinical_summary.csv")

# Apply final filters and priority sorting for high-confidence targets
high_priority_targets <- pathogenic_summary %>%
  filter(
    grepl("Pathogenic|Likely_pathogenic", CLNSIG), # Pathogenicity filter
    grepl("germline", Origins)                     # Germline origin filter
  ) %>%
  arrange(
    desc(Has_GRCh38),        # Priority to variants with GRCh38 coordinates
    desc(Submitter_Count),   # Then by evidence strength
    desc(Gene_Count),        # Then by number of genes affected
    CLNREVSTAT               # Finally by review status reliability
  )

# Save high-priority targets for RNA editing therapy
# write_csv(high_priority_targets, 
#           "high_priority_RNA_editing_therapeutic_targets.csv")

# Check how many candidate sites were found
# print(paste("Number of candidate targets found after correction:",
#             nrow(high_priority_targets)))

# Check the coverage of GRCh38 coordinates
# GRCh_stats <- high_priority_targets %>%
#   count(Has_GRCh38, Has_GRCh37)
# print("Genome assembly version coverage:")
# print(GRCh_stats)