# ADAR-Specific Oligonucleotide Design Pipeline

# Description: 
# This script analyzes ADAR substrate potential and designs antisense
# oligonucleotide for RNA editing applications. It fetches genomic sequences,
# predicts secondary structures, and evaluates editing site characteristics.
# This R script demonstrates an idea for training multi-modal models using
# A-to-I editing databases to automatically construct optimal ASO designs,
# with results intended for reference purposes only. Therefore, this script only
# presents a design logic and will not be maintained or updated in the future.

# Author: Zhikai Xu
# Date: 2025-10-05

# Load required packages
library(httr)
library(xml2)
library(dplyr)
library(stringr)
library(purrr)
library(readr)

# Retrieves DNA sequence from UCSC DAS server for specified genomic coordinates
fetch_genomic_sequence <- function(chrom, start, end) {
  
  # Construct URL for UCSC DAS server (hg38 genome assembly)
  base_url <- "http://genome.ucsc.edu/cgi-bin/das/hg38/dna"
  url <- sprintf("%s?segment=chr%s:%d,%d", base_url, chrom, start, end)
  
  tryCatch({
    
    # Send HTTP GET request with proper headers for UTF-8 encoding
    response <- GET(
      url,
      add_headers(
        "Accept-Charset" = "UTF-8",
        "Content-Type" = "text/xml; charset=utf-8"
      ),
      encoding = "UTF-8"
    )
    
    # Check if request was successful
    if (status_code(response) == 200) {
      
      # Extract content from response
      content <- content(response, "text", encoding = "UTF-8")
      
      # Force UTF-8 encoding
      if (Encoding(content) != "UTF-8") {
        content <- enc2utf8(content)
      }
      
      # Extract DNA sequence
      xml <- read_xml(content)
      dna <- xml_text(xml_find_first(xml, "//DNA"))
      
      # Remove newlines and digits, convert to uppercase
      clean_seq <- gsub("[\\n\\d]", "", dna, perl = TRUE)
      clean_seq <- toupper(clean_seq)
      
      # Ensure UTF-8 encoding
      if (Encoding(clean_seq) != "UTF-8") {
        clean_seq <- enc2utf8(clean_seq)
      }
      
      return(clean_seq)
    }
    
    return(NA)
  }, error = function(e) {
    message(
      sprintf(
        "Failed to fetch sequence %s:%d-%d: %s", chrom, start, end, e$message
      )
    )
    return(NA)
  })
}

# Runs RNAfold to predict secondary structure and free energy
predict_with_vienna <- function(seq) {
  
  # Create temporary files for input and output
  temp_in <- tempfile(fileext = ".fa")
  temp_out <- tempfile(fileext = ".txt")
  on.exit(unlink(c(temp_in, temp_out)))
  
  # Write sequence to temporary FASTA file
  writeLines(paste0(">ASO_sequence\n", seq), temp_in)
  
  exit_code <- system2(
    "RNAfold",
    args = c("--noPS", "--noLP", "-T", "37", paste0("--infile=", temp_in)),
    stdout = temp_out,
    stderr = FALSE
  )
  
  if (exit_code != 0) {
    stop(sprintf("RNAfold execution failed, exit code: %d", exit_code))
  }
  
  if (!file.exists(temp_out)) {
    stop("RNAfold did not generate output file")
  }
  
  # Read and parse RNAfold output
  output <- readLines(temp_out)
  
  if (length(output) < 2) {
    stop("RNAfold output has insufficient lines")
  }
  
  # Extract structure line containing secondary structure notation
  structure_line <- output[grepl("\\(", output)][1]
  if (is.na(structure_line) || !grepl("\\(", structure_line)) {
    stop("RNAfold output format abnormal - no structure found")
  }
  
  # Extract pure structure notation
  structure <- gsub("\\s*\\(.*", "", structure_line)
  structure <- gsub("\\s+", "", structure)
  
  # Extract free energy value
  energy_match <- regmatches(
    structure_line,
    regexpr("\\(\\s*([-+]?\\d+\\.?\\d*)\\s*\\)", structure_line)
  )
  
  free_energy <- ifelse(
    length(energy_match) > 0,
    as.numeric(gsub("[()]", "", energy_match)),
    NA_real_
  )
  
  # Calculate GC content for the sequence
  gc_content <- sum(str_count(seq, "[GC]")) / nchar(seq) * 100
  
  # Return comprehensive structure prediction results
  return(
    list(
      sequence = seq,
      secondary_structure = structure,
      free_energy = free_energy,
      gc_content = gc_content,
      structure_stability = ifelse(
        free_energy < -10,
        "High",
        ifelse(free_energy < -5, "Medium", "Low")
      ),
      method = "ViennaRNA"
    )
  )
}

# Wrapper function that validates input before structure prediction
predict_secondary_structure <- function(seq, aso_length = 20) {
  
  # Validate input sequence
  if (!is.character(seq) || nchar(seq) < 10) {
    stop("Sequence length must be ≥10nt")
  }
  
  # Call ViennaRNA prediction function
  result <- predict_with_vienna(seq)
  return(result)
}

# Calculates various sequence features relevant for ADAR editing efficiency
analyze_adar_substrate <- function(seq, ed_pos) {
  
  n <- nchar(seq)
  
  # Extract sequence windows around editing site (5bp and 10bp windows)
  win_5bp <- substr(seq, max(1, ed_pos - 5), min(n, ed_pos + 5))
  win_10bp <- substr(seq, max(1, ed_pos - 10), min(n, ed_pos + 10))
  
  # Calculate comprehensive feature set for ADAR substrate assessment
  features <- list(
    sequence_length = n,
    gc_content_total = sum(str_count(seq, "[GC]")) / n * 100,
    gc_content_5bp = sum(str_count(win_5bp, "[GC]")) / nchar(win_5bp) * 100,
    gc_content_10bp = sum(str_count(win_10bp, "[GC]")) / nchar(win_10bp) * 100,
    has_5_neighbor = ifelse(ed_pos > 1, substr(seq, ed_pos - 1, ed_pos - 1), NA),
    has_3_neighbor = ifelse(ed_pos < n, substr(seq, ed_pos + 1, ed_pos + 1), NA),
    u_content = sum(str_count(seq, "U")) / n * 100
  )
  
  return(features)
}

# Creates reverse complement sequence for complementarity analysis
reverse_complement <- function(seq) {
  complement_map <- c("A" = "T", "T" = "A", "G" = "C", "C" = "G")
  seq_chars <- strsplit(seq, "")[[1]]
  rev_seq <- rev(seq_chars)
  rev_complement <- sapply(rev_seq, function(x) complement_map[x])
  return(paste(rev_complement, collapse = ""))
}

# Measures potential for self-hairpin formation by comparing ends
calculate_complementarity_score <- function(seq) {
  
  n <- nchar(seq)
  if (n < 6) return(0) # Sequence too short for meaningful assessment
  
  # Extract 5' and 3' ends
  end_len <- min(5, floor(n/2))
  five_prime_end <- substr(seq, 1, end_len)
  three_prime_end <- substr(seq, n - end_len + 1, n)
  
  # Generate reverse complement of 3' end for comparison
  rev_three_prime_end <- reverse_complement(three_prime_end)
  
  # Count matching positions between 5' end and reverse complement of 3' end
  match_count <- 0
  for (i in 1:end_len) {
    if (substr(five_prime_end, i, i) == substr(rev_three_prime_end, i, i)) {
      match_count <- match_count + 1
    }
  }
  
  return(match_count / end_len)
}

# Generates candidate oligonucleotides centered around editing site
design_oligo <- function(seq, ed_pos, length = 20, n_candidates = 3) {
  
  n <- nchar(seq)
  candidates <- list()
  
  # Generate multiple candidate starting positions around editing site
  start_pos1 <- max(1, ed_pos - floor(length/2))
  end_pos1 <- min(n, start_pos1 + length - 1)
  if (end_pos1 - start_pos1 + 1 < length) {
    start_pos1 <- max(1, end_pos1 - length + 1)
  }
  
  start_pos2 <- max(1, ed_pos - floor(length/2) - 2)
  end_pos2 <- min(n, start_pos2 + length - 1)
  if (end_pos2 - start_pos2 + 1 < length) {
    start_pos2 <- max(1, end_pos2 - length + 1)
  }
  
  start_pos3 <- max(1, ed_pos - floor(length/2) + 2)
  end_pos3 <- min(n, start_pos3 + length - 1)
  if (end_pos3 - start_pos3 + 1 < length) {
    start_pos3 <- max(1, end_pos3 - length + 1)
  }
  
  # Select unique starting positions and limit to requested number
  start_poses <- unique(c(start_pos1, start_pos2, start_pos3))
  start_poses <- start_poses[1:min(n_candidates, length(start_poses))]
  
  # Generate candidate oligonucleotides for each starting position
  for (i in seq_along(start_poses)) {
    start_pos <- start_poses[i]
    end_pos <- min(n, start_pos + length - 1)
    
    # Ensure oligonucleotide has correct length
    if (end_pos - start_pos + 1 < length) {
      start_pos <- max(1, end_pos - length + 1)
    }
    
    oligo_seq <- substr(seq, start_pos, end_pos)
    edit_in_oligo <- (ed_pos >= start_pos) && (ed_pos <= end_pos)
    edit_pos_in_oligo <- ifelse(edit_in_oligo, ed_pos - start_pos + 1, NA)
    
    # Calculate self-complementarity score for stability assessment
    complementarity_score <- calculate_complementarity_score(oligo_seq)
    
    # Compile candidate information
    candidates[[i]] <- list(
      candidate_id = i,
      oligonucleotide_sequence = oligo_seq,
      oligonucleotide_length = nchar(oligo_seq),
      edit_position_in_oligo = edit_pos_in_oligo,
      start_position = start_pos,
      end_position = end_pos,
      edit_site_included = edit_in_oligo,
      complementarity_score = complementarity_score,
      has_strong_complementarity = complementarity_score >= 0.6
    )
  }
  
  # Rank candidates by complementarity score (lower is better for specificity)
  complementarity_scores <- sapply(
    candidates,
    function(x) x$complementarity_score
  )
  
  candidate_order <- order(complementarity_scores, decreasing = TRUE)
  candidates <- candidates[candidate_order]
  
  return(list(
    all_candidates = candidates,
    best_candidate = candidates[[1]]
  ))
}

# Main analysis function that processes multiple editing sites
analyze_adar_sites <- function(targets, context_size = 50, oligo_length = 20) {
  
  results <- list()
  
  # Process each editing site in the input targets
  for (i in 1:nrow(targets)) {
    chrom <- targets$GRCh38_Chromosome[i]
    start <- targets$GRCh38_Start[i]
    ref <- targets$REF[i]
    alt <- targets$ALT[i]
    
    # Progress reporting
    cat(sprintf("Processing site %d/%d: %s:%d %s>%s\n", 
                i, nrow(targets), chrom, start, ref, alt))
    
    tryCatch({
      
      # Fetch genomic sequence with extended context
      geno_seq <- fetch_genomic_sequence(
        chrom,
        start - context_size,
        start + context_size
      )
      
      # Validate retrieved sequence
      if (is.na(geno_seq) || nchar(geno_seq) < (2 * context_size + 1)) {
        warning(sprintf("Cannot get valid sequence %s:%d", chrom, start))
        next
      }
      
      # Extract core sequence around editing site for structure prediction
      tag_pos <- context_size + 1
      core_start <- max(1, tag_pos - 25)
      core_end <- min(nchar(geno_seq), tag_pos + 25)
      core_seq <- substr(geno_seq, core_start, core_end)
      edit_pos_in_core <- tag_pos - core_start + 1
      
      # Perform comprehensive analysis: structure, features, and oligo design
      structure_prediction <- predict_secondary_structure(core_seq)
      adar_features <- analyze_adar_substrate(core_seq, edit_pos_in_core)
      oligo_design <- design_oligo(core_seq, edit_pos_in_core, oligo_length)
      
      # Compile results for this site
      results[[i]] <- list(
        site_id = sprintf("%s:%d", chrom, start),
        genomic_coordinates = sprintf("%s:%d-%d", chrom, start, start),
        ref_allele = ref,
        alt_allele = alt,
        core_sequence = core_seq,
        edit_position_in_core = edit_pos_in_core,
        structure_prediction = structure_prediction,
        adar_features = adar_features,
        oligonucleotide_design = oligo_design,
        retrieval_status = "Success"
      )
      
    }, error = function(e) {
      
      # Store error information for failed analyses
      results[[i]] <- list(
        site_id = sprintf("%s:%d", chrom, start),
        genomic_coordinates = sprintf("%s:%d-%d", chrom, start, start),
        ref_allele = ref,
        alt_allele = alt,
        core_sequence = NA,
        structure_prediction = NA,
        adar_features = NA,
        oligonucleotide_design = NA,
        retrieval_status = paste("Error:", e$message)
      )
    })
  }
  
  return(results)
}

# Converts nested results to flat dataframe and saves to CSV
summarize_aso_design <- function(results) {
  
  # Convert nested list structure to flat dataframe
  summary <- map_dfr(results, function(x) {
    
    if (x$retrieval_status == "Success" && !is.null(x$structure_prediction)) {
      
      best_candidate <- x$oligonucleotide_design$best_candidate
      
      data.frame(
        site_id = x$site_id,
        genomic_coordinates = x$genomic_coordinates,
        ref_allele = x$ref_allele,
        alt_allele = x$alt_allele,
        core_sequence = x$core_sequence,
        edit_position = x$edit_position_in_core,
        secondary_structure = x$structure_prediction$secondary_structure,
        free_energy = x$structure_prediction$free_energy,
        gc_content_total = x$adar_features$gc_content_total,
        gc_content_5bp = x$adar_features$gc_content_5bp,
        gc_content_10bp = x$adar_features$gc_content_10bp,
        u_content = x$adar_features$u_content,
        oligonucleotide_sequence = best_candidate$oligonucleotide_sequence,
        oligo_length = best_candidate$oligonucleotide_length,
        edit_in_oligo = best_candidate$edit_site_included,
        edit_pos_in_oligo = best_candidate$edit_position_in_oligo,
        complementarity_score = best_candidate$complementarity_score,
        has_strong_complementarity = best_candidate$has_strong_complementarity,
        structure_stability = x$structure_prediction$structure_stability,
        prediction_method = x$structure_prediction$method,
        retrieval_status = x$retrieval_status,
        stringsAsFactors = FALSE
      )
    } else {
      
      # Create placeholder row for failed analyses
      data.frame(
        site_id = x$site_id,
        genomic_coordinates = x$genomic_coordinates,
        ref_allele = x$ref_allele,
        alt_allele = x$alt_allele,
        core_sequence = NA,
        edit_position = NA,
        secondary_structure = NA,
        free_energy = NA,
        gc_content_total = NA,
        gc_content_5bp = NA,
        gc_content_10bp = NA,
        u_content = NA,
        oligonucleotide_sequence = NA,
        oligo_length = NA,
        edit_in_oligo = NA,
        edit_pos_in_oligo = NA,
        complementarity_score = NA,
        has_strong_complementarity = NA,
        structure_stability = NA,
        prediction_method = NA,
        retrieval_status = x$retrieval_status,
        stringsAsFactors = FALSE
      )
    }
  })
  
  # Save results to CSV file
  write_csv(summary, "adar_analysis_results.csv")
  cat(sprintf("Results saved to: %s\n", "adar_analysis_results.csv"))
  
  # Print analysis statistics
  success_count <- sum(summary$retrieval_status == "Success")
  cat(sprintf("\nADAR Analysis Statistics:\n"))
  cat(sprintf("Total sites analyzed: %d\n", nrow(summary)))
  cat(sprintf("Successful analyses: %d\n", success_count))
  
  return(summary)
}

# Orchestrates the entire ADAR analysis workflow from input to output
aso_design_pipeline <- function(context_size = 50, oligo_length = 20) {
  
  # Validate input file existence
  if (!file.exists("annotated_targets_with_vcf_info.csv")) {
    stop("Input file 'annotated_targets_with_vcf_info.csv' not found!")
  }
  
  # Read and filter input data for ADAR editing sites
  top_targets <- read_csv(
    "annotated_targets_with_vcf_info.csv",
    show_col_types = FALSE
  ) %>% filter(RNA_Editing_Type == "ADAR_A_to_I")
  
  # Pipeline initialization and progress reporting
  cat("Starting ADAR-specific RNA editing site analysis...\n")
  cat(
    sprintf(
      "Input data contains %d potential editing sites\n", nrow(top_targets)
    )
  )
  
  # Execute analysis and generate summary
  results <- analyze_adar_sites(top_targets, context_size, oligo_length)
  summary <- summarize_aso_design(results)
  
  cat("ADAR analysis completed!\n")
  return(list(design_results = summary))
}

# Execute the complete analysis pipeline
design_results <- aso_design_pipeline(context_size = 50, oligo_length = 20)
print(head(design_results$design_results))