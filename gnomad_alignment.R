# Load necessary R packages
library(purrr)
library(readr)
library(dplyr)
library(data.table)

# Read target list containing high-priority RNA editing therapeutic candidates
targets <- read_csv(
  "high_priority_RNA_editing_therapeutic_targets.csv",
  show_col_types = FALSE
)

# Extract variant information from gnomAD VCF files using bcftools
query_vcf <- function(chrom, pos_start, pos_stop) {
  
  # Construct VCF filename using chromosome number
  vcf_file <- sprintf("gnomad.joint.v4.1.sites.chr%s.vcf.bgz", chrom)
  
  # Build bcftools query command with format specifiers
  cmd <- sprintf(
    "bcftools query -r 'chr%s:%s-%s' -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\t%%FILTER\\t%%INFO\\n' '%s'",
    chrom, pos_start, pos_stop, vcf_file
  )
  
  # Execute bcftools command
  result <- fread(
    cmd = cmd,
    sep = "\t",
    col.names = c("GRCh38_Chromosome", "GRCH38_POS", "REF", "ALT", "FILTER", "INFO"),
    na.strings = c(".", "NA")
  )
  
  # Handle empty results when no variants found in specified region
  if (nrow(result) == 0) {
    return(
      data.frame(
        GRCh38_Chromosome = character(),
        GRCH38_POS = integer(),
        REF = character(),
        ALT = character(),
        FILTER = character(),
        INFO = character()
      )
    )
  }
  
  # Clean chromosome names and adjust coordinate system
  result <- result %>%
    mutate(
      GRCh38_Chromosome = sub("^chr", "", GRCh38_Chromosome),
      GRCH38_POS = GRCH38_POS + 1
    )
  
  return(result)
  
}

# Execute Batch Queries Against gnomAD Database
# Note: gnomAD VCF uses 0-based coordinates, adjust target positions accordingly
vcf_data <- pmap_dfr(
  list(targets$GRCh38_Chromosome,
       targets$GRCh38_Start - 1, # Convert to 0-based start
       targets$GRCh38_Stop - 1), # Convert to 0-based end
  query_vcf # Apply query function to each target region
)

# Test query with first target for debugging purposes
# vcf_data <- pmap_dfr(
#   list(targets$GRCh38_Chromosome[1],
#        targets$GRCh38_Start[1]- 1,
#        targets$GRCh38_Stop[1]- 1),
#   query_vcf
# )
# Examine structure of results for validation
# str(vcf_data) 

# Merge results
final_data <- targets %>%
  left_join(
    vcf_data,
    by = c("GRCh38_Chromosome",            # Chromosome match
           "GRCh38_Start" = "GRCH38_POS",  # Position match
           "REF",                          # Reference allele match
           "ALT")                          # Alternate allele match
  )

# Save final annotated dataset for downstream analysis
write_csv(final_data, "gnomad_annotated_targets.csv")