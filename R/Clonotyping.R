# Clonotyping: From FASTQ files to clonotype tables
# To get the raw data (FASTQ files) please contact 
# Pastille, Eva <Eva.Pastille@uk-essen.de>

# load metadata
meta <- get(load("data/meta.RData"))

# create bash command to execute mixcr analysis
cmd <- paste0("mixcr analyze -t 20 -s mmu takara-mouse-rna-tcr-smarter ",
              meta$fullpath_fastq_R1, " ",
              meta$fullpath_fastq_R2, " ",
              "mixcr_output/", meta$sample_id)

# write bash command
writeLines(text = cmd, con = "step_mixcr.sh")

# Execute step_mixcr.sh with
system("bash step_mixcr.sh")

# Read data and format clonotype tables
fs <- list.files(path = "mixcr_output/", pattern = "clones_TR[AB]\\.tsv", full.names = T)

# Sort and rank clonotypes by expansion
mixcr_clones <- lapply(X = fs, FUN = function(x) {
  y <- read.csv(file = x, sep = "\t")
  y$file <- x
  y <- y[order(y$readCount, decreasing = T),]
  y$rank <- 1:nrow(y)
  return(y)
})
mixcr_clones <- do.call(rbind, mixcr_clones)

# Formatting, chain, experiment ID
mixcr_clones$locus <- ifelse(test = regexpr(pattern = "TRB", text = mixcr_clones$file)!=-1, yes = "TRB", no = "TRA")
mixcr_clones$experiment <- "2"
mixcr_clones$sample <- gsub(pattern = "mixcr_output\\/\\/|\\.clones|\\_TR[AB]\\.tsv", replacement = '', x = mixcr_clones$file)

# get meta data
meta <- get(load(file = "data/meta.RData"))
meta$sample <- meta$sample_id
meta$experiment <- 2
meta <- meta[, c("sample", "experiment", "group")]

# merge clonotypes and metadata
mixcr_clones <- merge(x = mixcr_clones, y = meta, by = c("experiment", "sample"), all.x = T)

# More formatting: productive vs. non productive, CD8 label, etc.
mixcr_clones$cell_type <- "CD8"
mixcr_clones$productive <- TRUE
mixcr_clones$productive[which(regexpr(pattern = "\\*|\\_|covered", text = mixcr_clones$aaSeqCDR3)!=-1)] <- FALSE
mixcr_clones$productive[mixcr_clones$aaSeqCDR3==""] <- FALSE
mixcr_clones <- mixcr_clones[mixcr_clones$productive==TRUE,]

# Save .RData
save(mixcr_clones, file = "data/mixcr_clones.RData")
