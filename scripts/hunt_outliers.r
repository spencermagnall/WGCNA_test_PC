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
#   pdf(file, width = 11, height = 5)
#   par(mfrow = c(1, 2))
#   plot_density(lcpm_before, sprintf("Before filtering (%d genes)", nrow(lcpm_before)),lcpm_cutoff)
#   plot_density(lcpm_after,  sprintf("After filtering (%d genes)",  nrow(lcpm_after)),lcpm_cutoff)
#   dev.off()

  # Bit sloppy with the variables but don't log for power fitting
  lcpm_after <- cpm(y, log=FALSE)
  return(lcpm_after)
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

plot_sample_tree <- function(dat, g, cut_height = NULL, outdir = "../figures") {
  sample_tree <- hclust(dist(dat), method = "average")

  png(file.path(outdir, paste0("sample_tree_", g, ".png")), width = 1400, height = 800)
  par(cex = 0.6, mar = c(9, 10, 7, 2))
  plot(sample_tree, main = paste("Sample clustering,", g),
       sub = "", xlab = "", hang = 0.05)

  if (!is.null(cut_height)) {
    abline(h = cut_height, col = "red")
  }
  dev.off()

  return(sample_tree)
}

# Transform the data to the format WGCNA expects
datExpr_list <- lapply(lcpm_by_group, t)
print("Transformed data")



for (g in names(datExpr_list)) {
  sample_tree <- plot_sample_tree(datExpr_list[[g]], g)
  cat("\n", g, "- top merge heights:\n")
  print(round(sort(sample_tree$height, decreasing = TRUE)[1:5], 1))
}

for (g in names(datExpr_list)) {
  dat_g      <- datExpr_list[[g]]
  sample_adj <- adjacency(t(dat_g), type = "distance")
  sample_k   <- rowSums(sample_adj) - 1
  z_k        <- (sample_k - mean(sample_k)) / sd(sample_k)

  cat("\n", g, "- lowest z_k:\n")
  print(round(sort(z_k)[1:3], 2))
}

mods <- readRDS("../networks/GG4_modules_sft7.rds")
length(setdiff(unique(mods$moduleColors), "grey"))
length(setdiff(unique(mods$mergedColors), "grey"))