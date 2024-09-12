# Microbiome-16srRNA-analysis-using-DADA2-and-Phyloseq
# This is my pipeline for analysis of 16srRNA sequences of the gut microbiota of several species of eukaryotes\
# refer to https://benjjneb.github.io/dada2/tutorial_1_8.html for detailed reference


####Analysis Performed on Rstudio (R version 4.3.2)
###setting the scene
>library(c(dada2, csv, phyloseq))# need a lot more packages
>path <- file.path("C:", "Users", "rthekkoot", "Downloads", "real-seq")##set file path to file where results of metagenomic sequencing is downloaded to

## Extract sample names
>sample.names <- sapply(strsplit(fnFs, "_S"), `[`, 1)
>sample.names <- sapply(strsplit(fnFs, "_S"), `[`, 1)

## Specify the full path to the fnFs and FnRs
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
>fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
>fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
>fnFs <- file.path(path, fnFs)
>fnRs <- file.path(path, fnRs)
>all(file.exists(fnFs, fnRs))
[1] TRUE
>fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
>fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
>plotQualityProfile(fnFs, aggregate=T) ## sequences need to be trimmed (usually they are around 300bp) as both forward and reverse will have low sequencing quality at the end, though forward reads will be higher in quality than reverse reads, due to the nature of illumina sequencing
#splitting sample names so that all of them are not the same eventually
>sample.names <- sapply(strsplit(files_list, "_S"), `[`, 1)
>sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)  ## 1 --> 2

#filter and trim
>filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
>filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
>names(filtFs) <- sample.names
>names(filtRs) <- sample.names
>out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(290,210),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, multithread=TRUE, trimLeft = c(21,21)) # On Windows set multithread=FALSE, trunclen is where the reads are cut forward, reverse respectively

#estimating error rates
>errF <- learnErrors(filtFs, multithread=TRUE)
>errR <- learnErrors(filtRs, multithread=TRUE)

#plotting error profile
>plotErrors(errF, nominalQ=TRUE)
#dereplication
>derepFs <- derepFastq(filtFs, verbose=TRUE)
>derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
>names(derepFs) <- sample.names
>names(derepRs) <- sample.names

#Applying the core DADA2 algorithm

>dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
>dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

>dadaFs[[1]]
#changing the parameter pool=pseudo,(pooling information across samples can increase sensitivity to sequence variants that may be present at very low frequencies in multiple samples) 
>dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool="pseudo")
>dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool="pseudo")
>dadaFs[[1]]

#merging parired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
>head(mergers[[1]])

#constructing sequence table
>seqtab <- makeSequenceTable(mergers)
>dim(seqtab)

#examining the sequence table
table(nchar(getSequences(seqtab)))

#Removing Chimeras
>seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

#tracking reads through the pipeline
>getN <- function(x) sum(getUniques(x))
>track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
>colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
>rownames(track) <- sample.names
head(track)

#Assigning taxonomy using the SILVA database
>taxa <- assignTaxonomy(seqtab.nochim, "/Users/rthekkoot/Downloads/real-seq/silva_nr_v128_train_set.fa.gz", multithread=TRUE)
#download database from SILVA site and place in path folder, to note: the database needs to be downloaded from dada2 website






