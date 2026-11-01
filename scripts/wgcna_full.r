library(edgeR)
library(SummarizedExperiment)
library(WGCNA)

# Setup input data, output dir etc.
rds_file <- "../prad_se.rds"
outdir <- "results"


# Not going to bother changing this but just know that its here
MIN_MOD <- 20       # minimum module size
MERGE_H <- 0.35     # module merge height
PRI_COL <- "primary_gleason_grade"
SEC_COL <- "secondary_gleason_grade"
nthr = 4
# Lets try to read in and then print our data 

samples_df <- readRDS(rds_file)

# Function to work out gleason grade 
get_gleason_grade <- function(samples_df){
    cd  <- as.data.frame(colData(samples_df))
    stopifnot(!is.null(cd$primary_gleason_grade), !is.null(cd$secondary_gleason_grade))

    num <- function(x) as.integer(gsub("\\D", "", as.character(x)))   # "Pattern 3" -> 3
    pri <- num(cd$primary_gleason_grade)
    sec <- num(cd$secondary_gleason_grade)
    stopifnot(length(pri) == ncol(samples_df), length(sec) == ncol(samples_df))
    tot <- pri + sec

    # cross-check against the curated total: this should be a clean diagonal
    table(tot, cd$gleason_score, useNA = "ifany")

    gg <- rep(NA_integer_, length(tot))
    gg[which(tot <= 6)]            <- 1L
    gg[which(tot == 7 & pri == 3)] <- 2L    # 3+4
    gg[which(tot == 7 & pri == 4)] <- 3L    # 4+3
    gg[which(tot == 8)]            <- 4L
    gg[which(tot >= 9)]            <- 5L
    gg <- factor(gg, levels = 1:5, labels = paste0("GG", 1:5))

    return(gg)
}
plot_density <- function(lcpm, title, lcpm_cutoff) {
  d <- apply(lcpm, 2, density)
  plot(NA, xlim = range(sapply(d, `[[`, "x")), ylim = range(sapply(d, `[[`, "y")),
       xlab = "log-CPM", ylab = "Density", main = title)
  for (i in seq_along(d)) lines(d[[i]], col = adjustcolor("steelblue", alpha.f = 0.25))
  abline(v = lcpm_cutoff, lty = 3)
}


# Filter and normalise function
filter_normalise <- function(samples_df, samples,g) {

  y <- DGEList(counts = assay(samples_df[, samples], "unstranded"),
               genes  = as.data.frame(rowData(samples_df)))

  # Edge R filtering code 
  # L <- mean(y$samples$lib.size)   * 1e-6                 # mean library size, millions
  # M <- median(y$samples$lib.size) * 1e-6                 # median library size, millions
  # lcpm_cutoff <- log2(10 / M + 2 / L)                    # the min.count = 10 cutoff on the log-CPM scale
  min_cpm = 1
  L <- mean(y$samples$lib.size) * 1e-6
  lcpm_cutoff <- log2(min_cpm + 2 / L)

  lcpm_before <- cpm(y, log = TRUE)

  # Keep genes with CPM >= min_cpm in every sample of this group
  keep <- rowSums(cpm(y) >= min_cpm) == ncol(y)
  message(g, ": ", sum(keep), " of ", nrow(y), " genes kept")
  y <- y[keep, , keep.lib.sizes = FALSE]

  # edgeR filtering code 
  # lcpm_before <- cpm(y, log = TRUE)
  # y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]  

  y <- calcNormFactors(y)
  lcpm_after <- cpm(y, log = TRUE)
  file <- paste0("../figures/filtering_density_", g, ".pdf")
  # Call the plotting function
  pdf(file, width = 11, height = 5)
  par(mfrow = c(1, 2))
  plot_density(lcpm_before, sprintf("Before filtering (%d genes)", nrow(lcpm_before)),lcpm_cutoff)
  plot_density(lcpm_after,  sprintf("After filtering (%d genes)",  nrow(lcpm_after)),lcpm_cutoff)
  dev.off()

  # Bit sloppy with the variables but don't log for power fitting
  lcpm_after <- cpm(y)
  return(lcpm_after)
}

# fit_power_thres <- function(dat, g){
  
#   powers <- 1:20                                     # candidate powers to test 20 is probably excessive
#   sft <- pickSoftThreshold(
#   dat,
#   powerVector = powers,
#   networkType = "signed",                          # negative correlation is NOT a connection
#   verbose = 5                                      # print progress
#   )
#   fit <- sft$fitIndices                              # one row per power tested


#   # SFT.R.sq is always positive. If the slope is positive the network is not scale-free,
#   # so flip the sign for those rows to make bad fits show up as negative.
#   fit$signed_R2 <- -sign(fit$slope) * fit$SFT.R.sq


#   fit[, c("Power", "signed_R2", "slope", "mean.k.", "median.k.", "max.k.")]

#   #Pick the lowest power that clears the threshold
#   r2_cut <- 0.80
#   passing <- fit$Power[fit$signed_R2 >= r2_cut]
#   power <- if (length(passing) > 0) min(passing) else NA
#   power                                              # NA means nothing reached the cutoff

#   file <- paste0("threshold_", g, ".pdf")
#   ## Plot both diagnostics
#   pdf(file, width = 10, height = 5)
#   par(mfrow = c(1, 2))
#   plot(fit$Power, fit$signed_R2, type = "n",
#      xlab = "Soft threshold (power)", ylab = "Scale-free fit (signed R^2)")
#   text(fit$Power, fit$signed_R2, labels = fit$Power, col = "red")
#   abline(h = r2_cut, col = "red")

#   plot(fit$Power, fit$mean.k., type = "n",
#      xlab = "Soft threshold (power)", ylab = "Mean connectivity")
#   text(fit$Power, fit$mean.k., labels = fit$Power, col = "red")
#   dev.off()

#   power
# }

fit_power_thres <- function(dat, g) {

  powers = c(c(1:10), seq(from = 10, to = 20, by = 2))

  sft <- pickSoftThreshold(
    dat,
    blockSize= 13000,
    powerVector = powers,
    networkType = "unsigned",
    verbose = 5
    
  )

  fit <- sft$fitIndices
  fit$signed_R2 <- -sign(fit$slope) * fit$SFT.R.sq

  print(fit[, c("Power", "signed_R2", "slope", "mean.k.")])

  r2_lo <- 0.8
  r2_hi <- 0.9

  in_window <- fit$Power[fit$signed_R2 > r2_lo & fit$signed_R2 < r2_hi]

  if (length(in_window) > 0) {
    power <- min(in_window)
  } else {
    above <- fit$Power[fit$signed_R2 > r2_lo]
    if (length(above) > 0) {
      power <- min(above)
    } else {
      power <- NA
    }
    message(g, ": no power in ", r2_lo, "-", r2_hi, "; fallback power = ", power)
  }

  pdf(paste0("../figures/threshold_", g, ".pdf"), width = 10, height = 5)
  par(mfrow = c(1, 2))
  plot(fit$Power, fit$signed_R2, type = "n",
       xlab = "Soft threshold (power)", ylab = "Scale-free fit (signed R^2)",
       main = g)
  text(fit$Power, fit$signed_R2, labels = fit$Power, col = "red")
  abline(h = c(r2_lo, r2_hi), col = "red")

  plot(fit$Power, fit$mean.k., type = "n",
       xlab = "Soft threshold (power)", ylab = "Mean connectivity", main = g)
  text(fit$Power, fit$mean.k., labels = fit$Power, col = "red")
  dev.off()

  return(power)
}

make_wgcna_dendrogram <- function(dat, sft_power, g){


  adjacency = adjacency(dat, power = sft_power)



  # Turn adjacency into topological overlap
  TOM = TOMsimilarity(adjacency)
  dissTOM = 1-TOM

  # Call the hierarchical clustering function
  geneTree = hclust(as.dist(dissTOM), method = "average")
  file = paste0(g,"Genetree_Plot.png")
  # Plot the resulting clustering tree (dendrogram)
  png(filename = file,
      width = 300,
      height = 200,
      units = "mm",
      res = 100)

  label = paste0("Gene clustering on TOM-based dissimilarity\n", g)
  
  plot(geneTree,
     xlab="", 
     sub="", 
     main = label,
     labels = FALSE, 
     hang = 0.04)

  dev.off()

# Want module sizes in the 20s 
minModuleSize = 30

# We'll do this sort of manually rather than using the automated WGCNA mergerging stuff
dynamicMods = cutreeDynamic(dendro = geneTree, 
                            distM = dissTOM,
                            deepSplit = 2, 
                            pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)

  # how many genes in each module?
  table(dynamicMods)

  # lables into colors
  dynamicColors = labels2colors(dynamicMods)

  # same as above but just for colours
  table(dynamicColors)

  file =  paste0(g,"_Genetree_Plot_With_Modules.png")
  # Plot the dendrogram w/ colours 
  png(filename = file,
      width = 300,
      height = 200,
      units = "mm",
      res = 100)

  label = paste0("Gene dendrogram and module colors\n",g)
  plotDendroAndColors(geneTree, 
                    dynamicColors, 
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = label)

  dev.off()

# Calculate eigengenes
  MEList = moduleEigengenes(dat, colors = dynamicColors)
  MEs = MEList$eigengenes

  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs)

  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average")

  file = paste0("./",g, "_Eigengenes_Dendogram.png")
  # Plot the eigengene dendogram
  png(filename = file,
    width = 300,
    height = 200,
    units = "mm",
    res = 100)
  label = paste0("Clustering of module eigengenes\n", g)
  plot(METree, 
     main = label,
     cex.main = 2,
     xlab = "", sub = "")                        

  dev.off()

  # Merge some modules 
  merge <- mergeCloseModules(dat, 
                           dynamicColors,         
                           verbose = 3)

  # The merged module colors
  mergedColors = merge$colors

  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs

  file = paste0("./",g,"Dendogram_And_Colours.png")
  png(filename = file,
      width = 300,
      height = 200,
      units = "mm",
      res = 100)

  label = paste0("TCGA-PRAD \n" ,g)
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    main = label,
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

  dev.off()



}


# Adjacency + TOM 
build_tom <- function(dat, power, network_type = "unsigned") {
  adj <- adjacency(dat, power = power, type = network_type)
  TOM <- TOMsimilarity(adj, TOMType = network_type)
  list(adjacency = adj, TOM = TOM)
}

find_modules <- function(dat, TOM, min_size = 30, deep_split = 2, merge_h = 0.25) {
  dissTOM  <- 1 - TOM
  geneTree <- hclust(as.dist(dissTOM), method = "average")

  moduleLabels <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                deepSplit = deep_split, pamRespectsDendro = FALSE,
                                minClusterSize = min_size)
  moduleColors <- labels2colors(moduleLabels)
  MEs          <- moduleEigengenes(dat, colors = moduleColors)$eigengenes

  merged <- mergeCloseModules(dat, moduleColors, verbose = 0)

  list(genes        = colnames(dat),
       geneTree     = geneTree,
       moduleLabels = moduleLabels,
       moduleColors = moduleColors,
       MEs          = MEs,
       mergedColors = merged$colors,
       mergedMEs    = merged$newMEs)
}



save_network <- function(tom, mods, g, outdir = ".") {
  saveRDS(tom,  file.path(outdir, paste0(g, "_tom.rds")))       # adjacency + TOM 
  saveRDS(mods, file.path(outdir, paste0(g, "_modules.rds")))   # tree, labels, colours, MEs 

}


plot_network <- function(mods, g, outdir = "../figures") {
  save_png <- function(name, expr) {
    png(file.path(outdir, paste0(g, "_", name, ".png")),
        width = 300, height = 200, units = "mm", res = 100)
    expr
    dev.off()
  }

  save_png("Genetree",
    plot(mods$geneTree, xlab = "", sub = "", labels = FALSE, hang = 0.04,
         main = paste0("Gene clustering on TOM-based dissimilarity\n", g)))

  save_png("Genetree_With_Modules",
    plotDendroAndColors(mods$geneTree, mods$moduleColors, "Dynamic Tree Cut",
                        dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
                        guideHang = 0.05,
                        main = paste0("Gene dendrogram and module colors\n", g)))

  save_png("Eigengene_Dendrogram",
    plot(hclust(as.dist(1 - cor(mods$MEs)), method = "average"),
         main = paste0("Clustering of module eigengenes\n", g),
         cex.main = 2, xlab = "", sub = ""))

  save_png("Dendrogram_And_Merged",
    plotDendroAndColors(mods$geneTree, cbind(mods$moduleColors, mods$mergedColors),
                        c("Dynamic Tree Cut", "Merged"),
                        main = paste0("TCGA-PRAD\n", g),
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05))
  invisible(NULL)
}


tum <- samples_df[, substr(colnames(samples_df), 14, 15) == "01"]
tum <- tum[, !duplicated(substr(colnames(tum), 1, 12))]    # Remove duplicates
nor <- samples_df[, substr(colnames(samples_df), 14, 15) == "11"]
samples_df  <- cbind(tum, nor) # combine gg samples and normal sample

gg <- get_gleason_grade(samples_df)
grp <- as.character(gg)
grp[substr(colnames(samples_df), 14, 15) == "11"] <- "Normal"

ok  <- !is.na(grp)
samples_df  <- samples_df[, ok]
grp <- factor(grp[ok], levels = c("Normal", paste0("GG", 1:5)))
table(grp)

lcpm_by_group <- list()


# Want to calculate lcpm before so that I can make the diagonstic plots

# Filter and normalise by group
# Not that levels here is needed to iterate over stuff
# probably closest to iters or enenmerate in python
for (g in levels(grp)) {
  samples <- which(grp == g)
  lcpm_by_group[[g]] <- filter_normalise(samples_df, samples, g)
}


# Transform the data to the format WGCNA expects
datExpr_list <- lapply(lcpm_by_group, t)
print("Transformed data")

# This is a sanity check for GG3 which is a bit shakey, do we need to remove outliers?  
dat3 <- datExpr_list[['GG3']]
colnames(dat3)
# Check data for excessive missing values and identify outlier samples
gsg = goodSamplesGenes(dat3)

# Check if its ok, if its true we can skip the next step:
gsg$allOK

#### Remove outlier samples ####
# Cluster the genes together using a distance matrix computation function and then clustering using hclust(). 
sampleTree = hclust(dist(dat3), method = "average")

# We want to export 
png(file = "./Gleason_Group_3_TCGA_PRAD_Clustering.png",
    width = 1400,
    height = 800,
    units = "px")

# Set text size
par(cex = 0.6)

# Set marginv: c(bottom, left, top, right)
par(mar = c(9,10,7,2))

# Create the plot:
plot(sampleTree,
     main = "TCGA-PRAD Clustering to Detect Outliers\n Gleason Group 3 (4+3)",
     labels = rownames(sampleTree),
     sub = "",
     xlab = "",
     cex.lab = 3,
     cex.axis = 3,
     cex.main = 3,
     cex = 2,
     hang = 0.05)

# Plot a line to show the cut
# This will need to change depending on the samples/genes used. 
# THIS NEEDS TO BE CHECKED IN EACH WORKFLOW
abline(h = 10000000, col = "red")
# I don't believe we need to remove outliers because it looks like the samples are quite closely clustered. 

dev.off()

# Sample-network connectivity (Oldham et al.): z-scored connectivity per sample
A   <- adjacency(t(dat3), type = "distance")     # samples x samples
k   <- rowSums(A) - 1
Z.k <- (k - mean(k)) / sd(k)
plot(Z.k, ylab = "Standardised connectivity (Z.k)", xlab = "Sample")
abline(h = -5, col = "red")                      # common cutoff; -2.5 is a looser one
names(Z.k)[Z.k < -5] 



# Now we'll make the power plots for each group, I really want to see that 
# match Mik's thesis in terms of ideal power that can be chosen. otherwise I'll do some fiddling with the counts

# Uncomment if you want to do it again. 
sft_powers <- numeric(length(datExpr_list))
names(sft_powers) <- names(datExpr_list)

for (g in names(datExpr_list)) {

  message("Fitting soft threshold for ", g)

  dat <- datExpr_list[[g]]
  power <- fit_power_thres(dat, g)

  sft_powers[g] <- power
}

# sft_powers

sft_powers <- list()
# Based on the threshold plots its fine to use the powers Mik used 
sft_powers$Normal <-8
sft_powers$GG1 <- 8
sft_powers$GG2 <- 6
sft_powers$GG3 <- 8 
sft_powers$GG4 <- 7 
sft_powers$GG5 <- 6 

# Now lets start making some networks 

# for (g in names(datExpr_list)){
#   dat = datExpr_list[[g]]
#   sft_power <- sft_powers[[g]]
#   sft_power 
#   make_wgcna_dendrogram(dat, sft_power, g)
# }


# SPLIT THE FUNCTION UP SO IT CAN BE USED to save eigengenes, TOM, etc. 
# Want to save 
# MEs, moduleLabels, moduleColors, geneTree, TOM, Adjaceny

network_dir = "../networks"
dir.create(network_dir, showWarnings = FALSE)

for (g in names(datExpr_list)) {
  dat <- datExpr_list[[g]]
  tom  <- build_tom(dat, sft_powers[[g]])
  mods <- find_modules(dat, tom$TOM)
  plot_network(mods, g)
  save_network(tom, mods, g, outdir = network_dir)
  rm(tom); gc()                       # garbage collect after saving 
  colnames(dat)
}



# TODO Phylostriagrphy 

# I'm going to start with GG1 as a test

# adjacency = adjacency(datExpr_list[["GG1"]], power = sft_powers$GG1)

# ### Topology Adjacency Into Topology Overlap ###

# # Turn adjacency into topological overlap
# TOM = TOMsimilarity(adjacency)
# dissTOM = 1-TOM

# # Call the hierarchical clustering function
# geneTree = hclust(as.dist(dissTOM), method = "average")

# # Plot the resulting clustering tree (dendrogram)
# png(filename = "./GG1_Genetree_Plot.png",
#     width = 300,
#     height = 200,
#     units = "mm",
#     res = 100)

# # Define the plot
# plot(geneTree,
#      xlab="", 
#      sub="", 
#      main = "Gene clustering on TOM-based dissimilarity\nGleason Group 1",
#      labels = FALSE, 
#      hang = 0.04)

# dev.off()

# # We would like large modules, so we set the minimum module size relatively high:
# minModuleSize = 30

# # Module identification using dynamic tree cut:
# dynamicMods = cutreeDynamic(dendro = geneTree, 
#                             distM = dissTOM,
#                             deepSplit = 2, 
#                             pamRespectsDendro = FALSE,
#                             minClusterSize = minModuleSize)

# # Visualise the number of genes in each module
# table(dynamicMods)

# # Convert numeric lables into colors
# dynamicColors = labels2colors(dynamicMods)

# # Visualise the number of genes in each colour
# table(dynamicColors)

# # Plot the dendrogram and colors underneath
# png(filename = "./GG1_Genetree_Plot_With_Modules.png",
#     width = 300,
#     height = 200,
#     units = "mm",
#     res = 100)

# # Define the plot
# plotDendroAndColors(geneTree, 
#                     dynamicColors, 
#                     "Dynamic Tree Cut",
#                     dendroLabels = FALSE,
#                     hang = 0.03,
#                     addGuide = TRUE,
#                     guideHang = 0.05,
#                     main = "Gene dendrogram and module colors\nGleason Group 1")

# dev.off()

# # Calculate eigengenes
# MEList = moduleEigengenes(datExpr_list[["GG1"]], colors = dynamicColors)
# MEs = MEList$eigengenes

# # Calculate dissimilarity of module eigengenes
# MEDiss = 1-cor(MEs)

# # Cluster module eigengenes
# METree = hclust(as.dist(MEDiss), method = "average")

# # Plot the eigengene dendogram
# png(filename = "./GG1_Eigengenes_Dendogram.png",
#     width = 300,
#     height = 200,
#     units = "mm",
#     res = 100)

# plot(METree, 
#      main = "Clustering of module eigengenes\nGleason Group 1",
#      cex.main = 2,
#      xlab = "", sub = "")                        

# dev.off()

# # Call an automatic merging function
# merge <- mergeCloseModules(datExpr_list[["GG1"]], 
#                            dynamicColors,         
#                            verbose = 3)

# # The merged module colors
# mergedColors = merge$colors

# # Eigengenes of the new merged modules:
# mergedMEs = merge$newMEs

# # We can write a script to export the script
# png(filename = "./GG1_Dendogram_And_Colours.png",
#     width = 300,
#     height = 200,
#     units = "mm",
#     res = 100)

# # Define the plot
# plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
#                     c("Dynamic Tree Cut", "Merged dynamic"),
#                     main = "TCGA-PRAD \n Gleason Group 1",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)

# dev.off()


