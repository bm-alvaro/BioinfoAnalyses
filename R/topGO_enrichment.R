# Alvaro Benitez, 2020

# User folder
setwd("~/GO_analysis")

# Mandatory library
library(topGO)

# Input files and (gene IDs and background annotation with Sma3s format)
file_genes <- "cmos_whole_list.txt"
file_background <- "new_go_grouped_cmos.txt"
Nodes <- 100 # number of processes to show
Ontology <- "GO.P.ID" #PFC (BP MF CC)

# Get gene IDs for the enrichment
genes <- read.csv(file_genes, header=F)$V1

# Get background annotation
GOesByID <- readMappings(file = "new_go_grouped_cmos.txt")
bg_genes <- names(GOesByID)

# Compare genes vs bg_genes
compared_genes <- factor(as.integer(bg_genes %in% genes))
names(compared_genes) <- bg_genes

# Create topGO object
GOdata <- new("topGOdata", ontology = "CC", allGenes = compared_genes,
              annot = annFUN.gene2GO, gene2GO = GOesByID)

# Run Fisher test
resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

# Create and print table with enrichment result
allRes <- GenTable(GOdata, classicFisher = resultFisher, topNodes = Nodes)
allRes

write.table(allRes, file = "CC_cmos_full_list.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
