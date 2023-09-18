
# Version 2:

# Install dada2
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("dada2")

# import dada2
library(dada2)

# Set “path” to Where the Sequences Are
path <- "/home/your_username/your_datafiles_path"
list.files(path)

# Have R Find and Sort Your Sequencing Files
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), '[', 1)
fnFs
fnRs
sample.names

# plot
plotQualityProfile(fnFs[1:2])

# python3 figaro/figaro.py --ampliconLength 470 --forwardPrimerLength 16 --reversePrimerLength 24

# Prepare To Run FilterAndTrim
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs)
names(filtRs)


# Trim And Filter
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(305,225), trimLeft = c(16,24), maxN=0, maxEE=c(5,4), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)
head(out)

# Learn Forward And Reverse Errors
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
errF

# Visualizing the Error Models
plotErrors(errF, nominalq=TRUE)

# Dereplicate Reads
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)

# Apply the Denoising
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]
dadaRs[[1]]

# Attempt to Merge Read Pairs
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])

# Inspect Our Denoised Amplicons
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

# Remove Chimeras and Examine Effect
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))
sum(seqtab.nochim)/sum(seqtab),

# Track Our Reads
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Assign Taxa
taxa <- assignTaxonomy(seqtab.nochim, "/home/your_username/your_datafiles_path/silva_nr_v138_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "/home/your_username/your_datafiles_path/silva_species_assignment_v138.fa.gz", multithread=TRUE)
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
write.table(taxa.print, "taxaResults.txt", sep="\t", row.names=TRUE, col.names=TRUE)






