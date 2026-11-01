
library(AnnotationDbi)
library(org.Hs.eg.db)





# Gene is UC or MC or NA (no data)


# PhyData cols - "GeneID"      "Entrez"      "Phylostrata"
# annot cols -  "ensg"   "entrez" "symbol"
get_phylostratum <- function(phyDat, annot) {
  age_from_entrez <- phyDat$Phylostrata[match(annot$entrez, phyDat$Entrez)]
  age_from_symbol <- phyDat$Phylostrata[match(annot$symbol, phyDat$GeneID)]

  phylostratum <- age_from_entrez                    # try Entrez first
  missing_age  <- is.na(phylostratum)
  phylostratum[missing_age] <- age_from_symbol[missing_age]   # fall back to the symbol

  phylostratum
}


# Label each gene UC, MC, or NA (no age)
classify_genes <- function(phyDat, annot) {
  phylostratum <- get_phylostratum(phyDat, annot)

  gene_class <- rep(NA_character_, length(phylostratum))       # start with NA for every gene

  has_age <- !is.na(phylostratum)
  is_old  <- has_age & phylostratum <= 3
  is_new  <- has_age & phylostratum >= 4

  gene_class[is_old] <- "UC"
  gene_class[is_new] <- "MC"

  gene_class
}


# 2. Test one module: build the 2x2 table, run the two one-sided tests
test_one_module <- function(gene_module, is_uc, module_name) {
  in_module <- gene_module == module_name            # which genes are in this module

  n_uc_in  <- sum(in_module  &  is_uc)               # in module, UC
  n_mc_in  <- sum(in_module  & !is_uc)               # in module, MC
  n_uc_out <- sum(!in_module &  is_uc)               # outside module, UC
  n_mc_out <- sum(!in_module & !is_uc)               # outside module, MC

  tab <- matrix(c(n_uc_in,  n_mc_in,
                  n_uc_out, n_mc_out),
                nrow = 2, byrow = TRUE)

  test_uc <- fisher.test(tab, alternative = "greater")   # more UC than expected?
  test_mc <- fisher.test(tab, alternative = "less")      # more MC than expected?

  module_size <- n_uc_in + n_mc_in
  odds_ratio  <- test_uc$estimate
  p_uc        <- test_uc$p.value
  p_mc        <- test_mc$p.value

  out <- data.frame(module = module_name,
                    size   = module_size,
                    n_uc   = n_uc_in,
                    OR     = odds_ratio,
                    p_uc   = p_uc,
                    p_mc   = p_mc)
  out
}


# 3. Test every module in one grade
test_modules <- function(phyDat, annot, moduleDat, g) {
  # length(annot)
  # length(moduleDat)
  # stopifnot(length(annot) == length(moduleDat))      # one module label per gene

  gene_class <- classify_genes(phyDat, annot)

  # keep only genes that have an age: these are the "bag"
  has_age     <- !is.na(gene_class)
  gene_module <- moduleDat[has_age]                  # module of each aged gene
  is_uc       <- gene_class[has_age] == "UC"

  modules <- unique(gene_module)
  modules <- modules[modules != "grey"]              # drop the unassigned genes

  per_module <- list()
  for (m in modules) {
    per_module[[m]] <- test_one_module(gene_module, is_uc, m)
  }
  module_results <- do.call(rbind, per_module)
  module_results$group <- g

  # multiple-testing correction across the modules in this grade
  module_results$fdr_uc <- p.adjust(module_results$p_uc, method = "BH")
  module_results$fdr_mc <- p.adjust(module_results$p_mc, method = "BH")

  # classify each module
  module_results$class <- "Mixed"
  module_results$class[module_results$fdr_mc < 0.05] <- "MC"
  module_results$class[module_results$fdr_uc < 0.05] <- "UC"

  module_results
}




# Load in the module data 
# Load in GG1 
mods <- readRDS(file.path("../networks", paste0("GG5_modules.rds")))

mods$mergedColors
names(mods$geneTree)             # slots of an hclust object: merge, height, order, labels, method, ...
is.null(mods$geneTree$labels)    # TRUE means no gene names were stored
head(mods$geneTree$labels)
genes <- mods$genes
#genes

# We have Gene ID, Enterz, and gene age 
# Everything >= 4 is considered MC 


# now per module, perform the fisher exact test to see if UC or MC enriched 
# Make sure to check how sensitive we are to p value cutoffs 

# Read in the phylostrata data
phyDat <- read.table("../gene_phylostrata_table.txt",
                          sep = ",",
                          stringsAsFactors = F,
                          header = T)
phyDat[,1] <- as.character(phyDat[,1])
phyDat[,2] <- as.character(phyDat[,2])
phyDat[,3] <- as.numeric(as.character(phyDat[,3]))
# Check the class of the columns
sapply(phyDat, class)
#phyDat
colnames(phyDat)
genes
phyDat$GeneID


# Convert ENSG to gene symbols
genes_ensg <- sub("\\..*$", "", mods$genes)                # drop the version suffix

entrez <- mapIds(org.Hs.eg.db, keys = genes_ensg, keytype = "ENSEMBL",
                 column = "ENTREZID", multiVals = "first")
symbol <- mapIds(org.Hs.eg.db, keys = genes_ensg, keytype = "ENSEMBL",
                 column = "SYMBOL",   multiVals = "first")

annot <- data.frame(ensg = genes_ensg, entrez = unname(entrez), symbol = unname(symbol))
saveRDS(annot, "gene_annotation.rds")
# annot     
annot <- readRDS("gene_annotation.rds")  
annot 
colnames(annot)


phylostratum <- get_phylostratum(phyDat, annot)
phylostratum
mean(!is.na(phylostratum))                           # fraction of network genes with an age

# For starters just try to classify all the genes 
clfy <- classify_genes(phyDat, annot) # Note we're using the correctly annotated genes

# Do we get junk??
#clfy


 
mods$mergedColors[1:3]          # should print three colours
modules <- mods$mergedColors
length(modules)


# Drop NA genes without ages 
has_age     <- !is.na(clfy)
gene_module <- modules[has_age]                  # module of each aged gene
is_uc       <- clfy[has_age] == "UC"

modules <- unique(gene_module)
modules <- modules[modules != "grey"]              # drop the unassigned genes

# Now pick one module and run the fisher test 
module_test = modules["turquoise"]
module_test

# test_mod_class <- test_one_module(gene_module,is_uc,"turquoise")
# test_mod_class
modules <- mods$mergedColors
test_modules <- test_modules(phyDat, annot, modules, "GG1")
test_modules
table(test_modules$class)

# # Run fisher test on one module for one grade (lets go with GG1)
# modules_gg1 = mergedColors[["GG1"]]
# # Think turquoise is a safe bet
# module_test = modules_gg1[["Turquoise"]]

# # Now fisher test 
# has_age     <- !is.na(gene_class)
# gene_module <- moduleDat[has_age]
