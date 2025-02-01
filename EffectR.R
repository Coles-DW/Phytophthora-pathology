setwd("")
install.packages("seqinr")
install.packages("ggplot2")
install.packages("devtools")
devtools::install_github(repo = "grunwaldlab/effectR", build_vignettes = TRUE)
library("usethis")
library("devtools")
library("effectR")
library("seqinr")
library("ggplot2")

#Data generation
fasta.file <- system.file("extdata", "Phyt_med.fasta", package = "effectR")
ORF <- seqinr::read.fasta(fasta.file)
#lines 16 and 17 don't work
fasta.file_SP <- system.file("extdata", "Phyt_med_Secre_Pro.fasta", package = "effectR")
ORF_SP <- seqinr::read.fasta(fasta.file)
#line 19 works
SP <- read.fasta("Phyt_med_Secre_Pro.fasta")

#Motif identification
#Work
REGEX_RXLR <- regex.search(ORF, motif="RxLR")
REGEX_RXLR_S <- regex.search(ORF_SP, motif="RxLR")
REGEX_CRN <- regex.search(ORF, motif="CRN")
REGEX_RXLR_SP <- regex.search(SP, motif="RxLR")
REGEX_CRN_SP <- regex.search(SP, motif="CRN")
#Don't work
REGEX_PAAR_SP <- regex.search(SP, motif="PAAR")
#work
REGEX_PAAR_SP <- regex.search(SP, motif="custom", reg.pat = "PAAR")
#Don't work
REGEX_DNLxxP_SP <- regex.search(SP, motif="custom", reg.pat = "DNLxxP")
REGEX_DNLxxP_SP <- regex.search(SP, motif="custom", reg.pat = "RxLR")
REGEX_EAR_SP <- regex.search(SP, motif="custom", reg.pat = "EAR")

#Hmm step
candidate.rxlr <- hmm.search(original.seq = fasta.file, regex.seq = REGEX_RXLR,
                                alignment.file=NULL, save.alignment=T, mafft.path = "", num.threads = 2,
                                hmm.path = "")

candidate.rxlr_SP <- hmm.search(original.seq = fasta.file_SP, regex.seq = REGEX_RXLR_SP,
                             alignment.file=NULL, save.alignment=T, mafft.path = "", num.threads = 2,
                             hmm.path = "")

#Summary and write to file
#RxLR
RxLR.effectors <- effector.summary(candidate.rxlr, motif = "RxLR")
length(RxLR.effectors$consensus.sequences)
head(RxLR.effectors$motif.table, n=5)
hmm.logo(candidate.rxlr$HMM_Table)
write.csv(RxLR.effectors$motif.table, "RxLR_effectors_motif_table.csv")
con.seq <- RxLR.effectors$consensus.sequences
write.csv(con.seq, "RxLR_effectors_consensus_sequences.csv")
#RxLR_SP
RxLR.effectors_SP <- effector.summary(candidate.rxlr_SP, motif = "RxLR")
length(RxLR.effectors_SP$consensus.sequences)
head(RxLR.effectors_SP$motif.table, n=5)
hmm.logo(candidate.rxlr_SP$HMM_Table)
write.csv(RxLR.effectors_SP$motif.table, "RxLR_effectors_motif_table_SP.csv")
con.seq_SP <- RxLR.effectors_SP$consensus.sequences
write.csv(con.seq_SP, "RxLR_effectors_consensus_sequences.csv")

#CRN

#CRN
REGEX_CRN <- regex.search(ORF, motif="CRN")
candidate.crn <- hmm.search(original.seq = fasta.file, regex.seq = REGEX_CRN,
                                alignment.file=NULL, save.alignment=T, mafft.path = "", num.threads = 2,
                                hmm.path = "")
#Error in hmm.search(original.seq = my_fasta3.0, regex.seq = my_REGEX_CRN,  : Not enough sequences for HMM step. At least 4 sequences are required.
CRN.effectors <- effector.summary(my_REGEX_CRN, motif = "CRN")
length(CRN.effectors$consensus.sequences)
head(CRN.effectors$motif.table, n=5)
hmm.logo(my_REGEX_CRN$HMM_Table)
#Error in apply(hmm.sums, 2, function(x) x/sum(x)) : dim(X) must have a positive length 
#In addition: Warning message: In max(x) : no non-missing arguments to max; returning -Inf
write.csv(RxLR.effectors$motif.table, "RxLR_effectors_motif_table.csv")
write.csv(RxLR.effectors$consensus.sequences, "RxLR_effectors_consensus_sequences.csv")