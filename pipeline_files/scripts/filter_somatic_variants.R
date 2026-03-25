library("VariantAnnotation")
library("stringr")
library("dplyr")
library("tidyr")
library("openxlsx")

# ─── PARSE ARGUMENTS ───────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 8) {
  stop("Usage: Rscript script.R <input_vcf.gz> <output.xlsx> <gene_bed> <ras_file> <mtor_file> <mcd_panelapp_file> <austin_file> <g4e_file>")
}

vcf_file <- args[1]
output_file <- args[2]
gene_file <- args[3]
ras_file <- args[4]
mtor_file <- args[5]
mcd_panelapp_file <- args[6]
austin_file <- args[7]
g4e_file <- args[8]

# ─── LOAD VCF ──────────────────────────────────────────────────────────────────

vcf_files <- c(vcf_file)

#Load gene list info
gene_info <- read.table(file = gene_file, sep = '\t', header = FALSE)
colnames(gene_info) <- c('Chr', 'Start', 'End', 'Gene')

ras <- read.table(file = ras_file, sep = '\t', header = FALSE)
mtor <- read.table(file = mtor_file, sep = '\t', header = FALSE)
mcd_panelapp <- read.table(file = mcd_panelapp_file, sep = '\t', header = FALSE)
austin_solid_tumour <- read.table(file = austin_file, sep = '\t', header = FALSE)
g4e <- read.table(file = g4e_file, sep = '\t', header = TRUE)
colnames(g4e)[7] <- 'Phenotype'
gene_info$Pathway <- ''

for (k in 1:nrow(gene_info)){
  gene <- gene_info[k,'Gene']
  pathways <- character(0)  # Initialize an empty character vector to store pathways
  
  if (gene %in% ras$V1){
    pathways <- c(pathways, 'RAS')
  }
  if (gene %in% mtor$V1){
    pathways <- c(pathways, 'MTOR')
  }
  if (gene %in% austin_solid_tumour$V1){
    pathways <- c(pathways, 'AustinPanel')
  }
  if (gene %in% g4e$Gene){
    pathways <- c(pathways, 'G4E')
  }
  if (gene %in% mcd_panelapp$V1){
    pathways <- c(pathways, 'MCD_PanelApp')
  }
  
  # Concatenate the pathways if multiple conditions are true
  gene_info[k,'Pathway'] <- paste(pathways, collapse = ', ')
}

ref <- "hg38"

# Initialize an empty list to hold VCF info for the summary sheet
all_variants <- list()

# Create an Excel workbook
wb <- createWorkbook()

# Helper function to reformat variants for Franklin URL
reformat_variants <- function(variants) {
  base_url <- "https://franklin.genoox.com/clinical-db/variant/snp/"
  urls <- paste(base_url, gsub(":", "-", variants), "-hg38", sep = "")
  return(urls)
}

for (vcf_file in vcf_files) {
  # Load VCF file
  vcf <- readVcf(vcf_file, ref)
  sample_name <- colnames(geno(vcf)$GT)  # Get the sample name from VCF
  
  print(sample_name)
  
  if (length(sample_name) == 0) {
    cat("No sample found in VCF:", vcf_file, "\n")
    next
  }

  # Load VCF into a dataframe
  vcf_info <- data.frame(
    CHR = as.character(seqnames(vcf)),
    POS = start(vcf),
    REF = as.character(ref(vcf)),
    ALT = sapply(alt(vcf), function(x) as.character(x[1]))
  )
  
  if (nrow(vcf) == 0)  {
    vcf_info$Variant = start(vcf)
    vcf_info$VAF = start(vcf)
    
    # Add the data to the all_variants list for summary
    all_variants[[sample_name]] <- vcf_info
    
    # Add a worksheet to the workbook for this sample
    addWorksheet(wb, sample_name)
    #writeData(wb, sample_name, vcf_info)
    writeDataTable(wb, sheet = sample_name, x = vcf_info, tableStyle = "TableStyleLight9")
    next
  }
  
  # Construct variant string
  vcf_info$Variant <- paste(vcf_info$CHR, vcf_info$POS, vcf_info$REF, vcf_info$ALT, sep = ":")
  
  # Extract the INFO fields
  info_fields <- data.frame(info(vcf))
  vcf_info <- cbind(vcf_info, info_fields)
  
  # Process AD field for REF and ALT allele depths
  ad <- geno(vcf)$AD
  if (!is.null(ad)) {
    ref_alt_counts <- ad[, 1]
    ref_ad_values <- sapply(ref_alt_counts, function(counts) if (is.na(counts[1])) NA else counts[1])
    alt_ad_values <- sapply(ref_alt_counts, function(counts) if (is.na(counts[2])) NA else counts[2])
    vcf_info$REF_AD <- ref_ad_values
    vcf_info$ALT_AD <- alt_ad_values
    vcf_info$VAF <- alt_ad_values / (ref_ad_values + alt_ad_values)
  }
  
  #Replicate variant rows have multiple VEP entries (i.e. variants that impact multiple genes)
  vcf_info <- vcf_info %>% 
    unnest(CSQ) %>% 
    distinct(CSQ, .keep_all = TRUE)
  
  ##VEP CSQ##
  #VEP consequences are split into individual fields
  if ("CSQ" %in% colnames(vcf_info)){
    print('Processing Ensembl VEP CSQ')
    csq = info(header(vcf))['CSQ',]$Description
    csq_col = unlist(strsplit(csq, ':'))[2]
    csq_fields = unlist(strsplit(csq_col, '\\|'))
    csq_data <- data.frame(matrix(nrow = nrow(vcf_info), ncol = length(csq_fields)))
    colnames(csq_data) <- csq_fields
    for (i in 1:nrow(vcf_info)){
      transcript <- unlist(vcf_info[i,]$CSQ)
      fields <- unlist(strsplit(transcript, split = '\\|'))
      for (col in 1:length(csq_fields)){
        csq_data[i,col] <- fields[col]}
    }
    vcf_info <- cbind(csq_data, vcf_info)}
  
  
  ## gnomAD ##
  if ("AC_gnomad_exomes_4.0" %in% colnames(vcf_info)){
    print('Processing AC_gnomad_exomes_4.0')
    for (i in 1:length(vcf_info$AC_gnomad_exomes_4.0)){
      entry = vcf_info$AC_gnomad_exomes_4.0[i]
      if (length(unlist(entry)) > 1) {
        max = 0
        for (j in 1:length(unlist(entry))){
          if (!is.na(unlist(entry)[j])){
            count = as.numeric(unlist(entry)[j])
            if (count > max){
              max = count}
          }
          entry = max}
      } else if (is.na(entry)){
        entry = 0
      }
      vcf_info$AC_gnomad_exomes_4.0[i] = as.numeric(entry)
    }
    vcf_info$AC_gnomad_exomes_4.0 <- as.numeric(vcf_info$AC_gnomad_exomes_4.0)
  }
  
  vcf_info$AC_gnomad_exomes_4.0
  
  
  if ("AC_gnomad_genomes_4.0" %in% colnames(vcf_info)){
    print('Processing AC_gnomad_genomes_4.0')
    for (i in 1:length(vcf_info$AC_gnomad_genomes_4.0)){
      entry = vcf_info$AC_gnomad_genomes_4.0[i]
      if (length(unlist(entry)) > 1) {
        max = 0
        for (j in 1:length(unlist(entry))){
          if (!is.na(unlist(entry)[j])){
            count = as.numeric(unlist(entry)[j])
            if (count > max){
              max = count}
          }
          entry = max}
      } else if (is.na(entry)){
        entry = 0
      }
      vcf_info$AC_gnomad_genomes_4.0[i] = as.numeric(entry)
    }
    vcf_info$AC_gnomad_genomes_4.0 <- as.numeric(vcf_info$AC_gnomad_genomes_4.0)
  }
  
  vcf_info$AC_gnomad_total_4.0 <- as.numeric(vcf_info$AC_gnomad_exomes_4.0) + as.numeric(vcf_info$AC_gnomad_genomes_4.0)
  
  if ("nhomalt_gnomad_exomes_4.0" %in% colnames(vcf_info)){
    print('Processing nhomalt_gnomad_exomes_4.0')
    for (i in 1:length(vcf_info$nhomalt_gnomad_exomes_4.0)){
      entry = vcf_info$nhomalt_gnomad_exomes_4.0[i]
      if (length(unlist(entry)) > 1) {
        max = 0
        for (j in 1:length(unlist(entry))){
          if (!is.na(unlist(entry)[j])){
            count = as.numeric(unlist(entry)[j])
            if (count > max){
              max = count}
          }
          entry = max}
      } else if (is.na(entry)){
        entry = 0
      }
      vcf_info$nhomalt_gnomad_exomes_4.0[i] = as.numeric(entry)
    }
    vcf_info$nhomalt_gnomad_exomes_4.0 <- as.numeric(vcf_info$nhomalt_gnomad_exomes_4.0)
  }
  
  if ("nhomalt_gnomad_genomes_4.0" %in% colnames(vcf_info)){
    print('Processing nhomalt_gnomad_genomes_4.0')
    for (i in 1:length(vcf_info$nhomalt_gnomad_genomes_4.0)){
      entry = vcf_info$nhomalt_gnomad_genomes_4.0[i]
      if (length(unlist(entry)) > 1) {
        max = 0
        for (j in 1:length(unlist(entry))){
          if (!is.na(unlist(entry)[j])){
            count = as.numeric(unlist(entry)[j])
            if (count > max){
              max = count}
          }
          entry = max}
      } else if (is.na(entry)){
        entry = 0
      }
      vcf_info$nhomalt_gnomad_genomes_4.0[i] = as.numeric(entry)
    }
    vcf_info$nhomalt_gnomad_genomes_4.0 <- as.numeric(vcf_info$nhomalt_gnomad_genomes_4.0)
  }
  
  vcf_info$nhomalt_gnomad_total_4.0 <- as.numeric(vcf_info$nhomalt_gnomad_exomes_4.0) + as.numeric(vcf_info$nhomalt_gnomad_genomes_4.0)
  
  # Match VCF genes with Epilepsy Gene List genes
  matches <- match(vcf_info$SYMBOL, gene_info$Gene)
  
  # Subset gene_info based on the matching rows
  gene_summary <- gene_info[matches,]
  
  # Add the matched gene summary to vcf_info
  vcf_info <- cbind(vcf_info, gene_summary)
  
  # Merge with g4e to add Inheritance and Phenotype without duplicating columns
  vcf_info <- merge(vcf_info, g4e[, c("Gene", "Inheritance", "Phenotype")], by.x = "SYMBOL", by.y = "Gene", all.x = TRUE, suffixes = c("", "_g4e"))
  
  #Remove repeat entries (i.e. same variant across multiple transcripts). 
  # Find all duplicated variants and their corresponding rows
  rownames(vcf_info) <- seq(nrow(vcf_info))
  
  duplicate_variant_index_1 <- which(duplicated(vcf_info$Variant) | duplicated(vcf_info$Variant, fromLast = TRUE))

  #First, remove any duplicates that have a blank SYMBOL column (i.e. no gene name)
  vcf_info <- subset(vcf_info, !((rownames(vcf_info) %in% duplicate_variant_index_1) & (SYMBOL == "")))
  
  # Reformat variant URLs for Franklin lookup
  vcf_info$Franklin_URL <- reformat_variants(vcf_info$Variant)
  
  vcf_info <- vcf_info %>%
  dplyr::select(dplyr::any_of(c(
    "CHR", "POS", "REF", "ALT", "Variant", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature", 
    "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp", "HGVSg", "Existing_variation", "DISTANCE", 
    "FLAGS", "VARIANT_CLASS", "CADD_Score", "SIFT", "PolyPhen", "ClinVarDN", "ClinVarGene", "ClinVarSIG", "ClinVarSIGCONF", 
    "AC_gnomad_exomes_4.0", "AC_gnomad_genomes_4.0", "AC_gnomad_total_4.0", 
    "nhomalt_gnomad_exomes_4.0", "nhomalt_gnomad_genomes_4.0", "nhomalt_gnomad_total_4.0", 
    "REF_AD", "ALT_AD", "VAF", "Gene.1", "Tool_Count", "Mutect", "Freebayes", "Strelka",
    "Pathway", "Inheritance", "Phenotype", "Franklin_URL", "BQBZ", "DP", "DP4", "MQ", "MQBZ", "MQSBZ", "RPBZ", "SCBZ", "SGB", "VDB"
  )))

  
  # Add the data to the all_variants list for summary
  all_variants[[sample_name]] <- vcf_info
  
  # Add a worksheet to the workbook for this sample
  addWorksheet(wb, sample_name)
  #writeData(wb, sample_name, vcf_info)
  writeDataTable(wb, sheet = sample_name, x = vcf_info, tableStyle = "TableStyleLight9")
}

# Create a summary sheet
#summary_df <- all_variants[[1]][, c("Variant", "SYMBOL"), drop = FALSE]
#for (sample_name in names(all_variants)) {
#  sample_data <- all_variants[[sample_name]][, c("Variant", "SYMBOL", "VAF")]
#  colnames(sample_data)[-1] <- paste(sample_name, c("VAF"), sep = "_")
#  summary_df <- merge(summary_df, sample_data, by = c("Variant", "SYMBOL"), all = TRUE)
#}

# Add the summary sheet to the workbook
#addWorksheet(wb, "Summary")
#writeDataTable(wb, sheet = "Summary", summary_df, tableStyle = "TableStyleLight9")


saveWorkbook(wb, output_file, overwrite = TRUE)

cat("Workbook created successfully with additional INFO fields, VEP annotations, and gnomAD data.\n")
