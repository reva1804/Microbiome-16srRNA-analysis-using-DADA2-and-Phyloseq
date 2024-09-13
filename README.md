# Microbiome-16srRNA-analysis-using-DADA2-and-Phyloseq
# This is my pipeline for analysis of 16srRNA sequences of the gut microbiota of several species of eukaryotes, useful for systems with low RAM
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
> mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
> head(mergers[[1]])

#constructing sequence table
> seqtab <- makeSequenceTable(mergers)
> dim(seqtab)

#examining the sequence table
>table(nchar(getSequences(seqtab)))

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
#download database from SILVA site and place in path folder, to note: the database needs to be downloaded from dada2 website or it needs to be formatted

> taxa.print <- taxa
> rownames(taxa.print) <- NULL
> head(taxa.print)#removing row names to display the taxa

#Adding species level assignment, as it takes a lot of computational power, the species assignment is done in chunks
> chunk.size <- 4000
> chunks <- split(c(1:nrow(taxa)), sort(c(1:nrow(taxa))%%ceiling(nrow(taxa)/chunk.size)))
> chunks.species <- lapply(chunks,function(x){return(addSpecies(taxa[x,], refFasta = "/Users/rthekkoot/Downloads/real-seq/silva_species_assignment_v138.1.fa.gz", verbose = TRUE))})

#Displayingg species level assignment
> taxa.species <- do.call(rbind, chunks.species)
> taxa.species.print <- taxa.species
> rownames(taxa.species.print) <- NULL
> head(taxa.species.print)

###Moving on to phyloseq, creating the phyloseq object
>library(“phyloseq”,“vegan”,”picente”,”rstatix”,”ggplot2”,”ggpubr”,”ape”,”adephylo”,”phytools”,”tidyverse”,”csv”)

#saving the output of DADA2
>saveRDS(taxa.species, "C:/Users/rthekkoot/Downloads/practice-seq/taxa.species.rds")
>saveRDS(sequence_table,"C:/Users/rthekkoot/Downloads/practice-seq/seqtab.nochim.rds")

##assigning the metadata table, make sure column 1 is the sample names as per the output of DADA2 in the same order, ie. sample names need to match
>sdata <-as.csv("C:/Users/rthekkoot/Downloads/practice-seq/metadata-toby.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE)

#removing taxa with 0 counts
> seqtab.nochim <- as.matrix(seqtab.nochim)
> m <- (colSums(seqtab.nochim, na.rm=TRUE) != 0) 
> nonzero <- seqtab.nochim[, m]

#format for phyloseq
samdata = sample_data(sdata)
seqtab.nochim = otu_table(nonzero, taxa_are_rows = FALSE)
taxa.species = tax_table(taxa.species)

##Assembling the phyloseq object
> ps = phyloseq(otu_table(seqtab.nochim), tax_table(taxa.species), sample_data(samdata))
> ps

##Adding the phylogenetic tree
##align the sequences
>library(DECIPHER)
> seqs <- getSequences(sequence_table)
> names(seqs) <- seqs
> alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA))

##Making phylogenetic tree with Fasttree, first download fasttree and place into directory

##saving the alignment in fasta format
>Biostrings::writeXStringSet(alignment, "aligned_1.fasta")
>system2("fasttree", args = c("-gtr -nt", "-out dada2_nochim.tree", "alignment1.fasta"))
>ape::read.tree("dada2_nochim.tree")
>tree <- read_tree("dada2_nochim.tree")

#Assembling the phyloseq object again
>ps = phyloseq(otu_table(seqtab.nochim), tax_table(taxa.species), sample_data(samdata),phy_tree(tree))

#rooting the tree, tree needs to be rooted for calculating UNIFRAC/WUNIFRAC distances
>phy_tree(ps) <- ape::root(phy_tree(ps), sample(taxa_names(ps), 1), resolve.root=TRUE)
>is.rooted(phy_tree(ps))

##Running DECONTAM to remove contaminants based on prevalence
> library(decontam)
> sample_data(ps)$is.neg <- sample_data(ps)$Sample_or_Control == "Negative_control"
##name all the negative controls Negative_control under column Sample_or_Control in the metadata
> contamdf.prev <- isContaminant(ps, method="prevalence", neg="is.neg")
> table(contamdf.prev$contaminant)
> head(which(contamdf.prev$contaminant))
> ps.noncontam <- prune_taxa(!contamdf.prev$contaminant, ps)
> ps.noncontam

##removing mitochondria, eukaryotes and chloroplasts from the phyloseq object
> ps.clean <- subset_taxa(ps.noncontam, Kingdom == "Bacteria") %>%  subset_taxa(!is.na(Phylum)) %>%subset_taxa(!Class %in% c("Chloroplast")) %>% subset_taxa(!Family %in% c("mitochondria"))
> ps.clean

## The Phyloseq object is now ready for analysis, please contact me if you have any doubts or for downstream analysis, credit to Dr. Christin Manthey for helping me with this pipeline.











