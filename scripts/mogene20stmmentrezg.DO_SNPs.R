# Assess how MoGene 2.0 Gene ST Brainarray Entrez Gene probeset mapping
# is affected by SNPs in Diversity Outbred (DO) mouse strain

# Adam Gower, for Gillian Beamer and Dan Gatti

# Global variables
brainarray.version <- "17.0.0"

# Define paths
bumsr.path <- "/restricted/projectnb/bumsr"
annotations.path <- file.path(bumsr.path, "Annotations")
working.path <- file.path(annotations.path, "DO_SNPs")
brainarray.path <- file.path(
  bumsr.path, "Brainarray", paste0("R-", getRversion()), brainarray.version
)

# Load and list packages
.libPaths(brainarray.path)
library(mogene20stmmentrezgprobe)
library(mogene20stmmentrezg.db)
library(Biostrings)
sessionInfo()

# File obtained from Thermo Fisher at:
# https://sec-assets.thermofisher.com/TFS-Assets/LSG/Support-Files/MoGene-2_0-st-v1-mm10-probe-tab.zip
# Note: cannot download file directly in this script,
#       as this requires authentication from Thermo Fisher
probetab.filename <- file.path(
  annotations.path, "MoGene-2_0-st-v1.mm10.probe.tab"
)
# Parse Affymetrix probe table
probetab <- read.delim(probetab.filename, stringsAsFactors=FALSE)
# Reduce to description of unique probe IDs/sequences
# (some probes map to more than one location/SNP)
probetab <- unique(
  probetab[c("Probe.ID", "probe.sequence", "target.strandedness")]
)
# Ensure all probes are in sense orientation
if (any(probetab$target.strandedness == "Antisense")) {
  i <- which(probetab$target.strandedness == "Antisense")
  probetab$probe.sequence[i] <- as.character(
    reverseComplement(DNAStringSet(probetab$probe.sequence[i]))
  )
  probetab$target.strandedness[i] <- "Sense"
}

# File received by email from Dan Gatti
snp.csv.filename <- file.path(working.path, "snps_in_probes_within.csv")
# Parse file generated by Dan Gatti
snp.table <- read.csv(snp.csv.filename, stringsAsFactors=FALSE)
# Merge SNP table with Affymetrix probe annotation on Probe.ID
snp.table <- merge(snp.table, probetab)

# Initialize a probe-level table from the Brainarray probe annotation table
probes <- mogene20stmmentrezgprobe
# Coerce columns to character as needed
probes$sequence <- as.character(probes$sequence)
probes$Probe.Set.Name <- as.character(probes$Probe.Set.Name)
probes$Target.Strandedness <- as.character(probes$Target.Strandedness)
# Ensure all probes are in sense orientation
if (any(probes$Target.Strandedness == "Antisense")) {
  i <- which(probes$Target.Strandedness == "Antisense")
  probes$sequence[i] <- as.character(
    reverseComplement(DNAStringSet(probes$sequence[i]))
  )
  probes$Target.Strandedness[i] <- "Sense"
}
# Indicate whether each probe is affected by a SNP
probes$SNP <- is.element(probes$sequence, snp.table$probe.sequence)
# Store table as RDS file
rds.filename <- file.path(working.path, "mogene20stmmentrezg.probes.RDS")
saveRDS(probes, rds.filename)

# Initialize a gene-specific table, collapsing probes by probeset ID
genes <- aggregate(probes["SNP"], by=probes["Probe.Set.Name"], FUN=c)
# Add gene-level annotation
genes <- cbind(
  genes,
  Symbol = unlist(mget(genes$Probe.Set.Name, mogene20stmmentrezgSYMBOL)),
  Description = unlist(mget(genes$Probe.Set.Name, mogene20stmmentrezgGENENAME)),
  stringsAsFactors = FALSE
)
# Tabulate total number of probes and total number of SNP-affected probes
genes$n.probes <- sapply(genes$SNP, length)
genes$n.SNP.probes <- sapply(genes$SNP, sum)
# Remove temporary 'SNP' column
genes$SNP <- NULL
# Store table as RDS file
rds.filename <- file.path(working.path, "mogene20stmmentrezg.genes.RDS")
saveRDS(genes, rds.filename)
# Reformat gene symbols as formulas (to avoid conversion, e.g., DEC1 -> 12/01)
options(useFancyQuotes=FALSE)
genes$Symbol <- paste0("=", dQuote(genes$Symbol))
# Store table as TSV file
tsv.filename <- file.path(working.path, "mogene20stmmentrezg.genes.txt")
write.table(
  genes, tsv.filename, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE
)

# Create a histogram of the fraction of probes affected by SNPs per probeset
pdf.filename <- file.path(working.path, "mogene20stmmentrezg.histogram.pdf")
pdf(pdf.filename, width=8.5, height=8.5)
hist(
  genes$n.SNP.probes / genes$n.probes, breaks=100,
  main="Effect of DO SNPs on mogene20stmmentrezg probesets",
  xlim=c(0,1), xlab="Fraction of probes affected by SNPs",
  ylab="Number of Entrez Genes"  
)
dev.off()
