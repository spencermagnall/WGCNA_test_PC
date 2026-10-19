library(edgeR)
library(SummarizedExperiment)
library(WGCNA)

# Setup input data, output dir etc.
rds_file <- "prad_se.rds"
outdir <- "results"


# Not going to bother changing this but just know that its here
N_GENES <- 5000   # top most-variable genes (MAD) used for the network
MIN_MOD <- 20       # minimum module size
MERGE_H <- 0.35     # module merge height
PRI_COL <- "primary_gleason_grade"
SEC_COL <- "secondary_gleason_grade"
nthr = 4
# Lets try to read in and then print our data 

samples_df <- readRDS(rds_file)

# Just going to do the head for obvious reasons 
#head(samples_df, n = 6)

# Print the entire thing why not?
#samples_df


# Remove genes that are poorly expressed 

# Use edgeR I guess
# Just get the primary tumor
samples_df <- samples_df[, samples_df$sample_type == "Primary Tumor"]
# remove duplicates
samples_df <- samples_df[, !duplicated(substr(colnames(samples_df), 1, 12))]

counts <- assay(samples_df, "unstranded")
# Check that this works
#class(counts[1, 1])        # should be integer-like
#head(colSums(counts))      # library sizes, tens of millions per sample




cd <- as.data.frame(colData(samples_df))      # after the tumour filtering step
for (v in c("primary_gleason_grade", "secondary_gleason_grade", "gleason_score",
            "paper_Reviewed_Gleason", "paper_Reviewed_Gleason_category",
            "paper_Reviewed_Gleason_sum")) {
  cat("\n==", v, "==\n")
  print(table(cd[[v]], useNA = "ifany"))
}

# claude slop
# Basic gist is to find the primary and secodary gleason grades and sum them to get a group
# The way its doing this is opaque to me to say the least
#  Grade groups (set your two column names)
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
table(gg, useNA = "ifany")

ok <- !is.na(gg)
samples_df<- samples_df[, ok]; gg <- droplevels(gg[ok])

y <- DGEList(counts = assay(samples_df, "unstranded"), genes = as.data.frame(rowData(samples_df)))
length(gg); ncol(y)   



# Claude wrote this 
# --- BEFORE filtering: record what we need to plot ---
lcpm_before <- cpm(y, log = TRUE)                      # prior.count = 2 by default
L <- mean(y$samples$lib.size)   * 1e-6                 # mean library size, millions
M <- median(y$samples$lib.size) * 1e-6                 # median library size, millions
lcpm_cutoff <- log2(10 / M + 2 / L)                    # the min.count = 10 cutoff on the log-CPM scale

print("Before gg stuff")

# Filter, include gleason grade in the filtering 
keep <- filterByExpr(y, group=gg)
y <- y[keep, , keep.lib.sizes = FALSE]
lcpm_after <- cpm(y, log = TRUE)

plot_density <- function(lcpm, title) {
  d <- apply(lcpm, 2, density)
  plot(NA, xlim = range(sapply(d, `[[`, "x")), ylim = range(sapply(d, `[[`, "y")),
       xlab = "log-CPM", ylab = "Density", main = title)
  for (i in seq_along(d)) lines(d[[i]], col = adjustcolor("steelblue", alpha.f = 0.25))
  abline(v = lcpm_cutoff, lty = 3)
}



pdf("filtering_density.pdf", width = 11, height = 5)
par(mfrow = c(1, 2))
plot_density(lcpm_before, sprintf("Before filtering (%d genes)", nrow(lcpm_before)))
plot_density(lcpm_after,  sprintf("After filtering (%d genes)",  nrow(lcpm_after)))
dev.off()

# Now normalise the data 
y <- calcNormFactors(y, method = "TMM")


# pick some with different behaviour 
genes <- c("KLK3", "AR", "MKI67", "TOP2A", "GAPDH")
idx <- match(genes, y$genes$gene_name)               # y is the filtered, TMM-normalised DGEList
if (anyNA(idx)) stop("Not in filtered data: ", paste(genes[is.na(idx)], collapse = ", "))

# before = same data with the normalisation factors reset to 1
y_before <- y
y_before$samples$norm.factors <- 1


lcpm5_before <- cpm(y_before, log = TRUE)[idx, ]     # 5 genes x samples
lcpm5_after  <- cpm(y,        log = TRUE)[idx, ]
rownames(lcpm5_before) <- rownames(lcpm5_after) <- genes

# same y-axis so the panels are comparable
yl <- range(lcpm5_before, lcpm5_after)

pdf("normalisation_genes.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
boxplot(t(lcpm5_before), ylim = yl, col = "grey80", las = 2,
        main = "Before TMM", ylab = "log-CPM")
boxplot(t(lcpm5_after),  ylim = yl, col = "steelblue", las = 2,
        main = "After TMM")
dev.off()
# Bit sloppy but whateves 
lcpm <- cpm(y, log = TRUE)

top <- order(apply(lcpm, 1, mad), decreasing = TRUE)[1:N_GENES]
dat <- t(lcpm[top, ])
dim(dat)            # check: ~468 x 5000


dim(dat)                                  # expect 468 x 5000
sum(is.na(dat)); sum(is.infinite(dat))    # expect 0 and 0
sum(apply(dat, 2, sd) == 0)               # expect 0
summary(as.vector(dat))                   # expect roughly -3 to 15 (log2 CPM), no NaN


powers <- 1:20                                     # candidate powers to test 20 is probably excessive
sft <- pickSoftThreshold(
  dat,
  powerVector = powers,
  networkType = "signed",                          # negative correlation is NOT a connection
  verbose = 5                                      # print progress
)
fit <- sft$fitIndices                              # one row per power tested


# SFT.R.sq is always positive. If the slope is positive the network is not scale-free,
# so flip the sign for those rows to make bad fits show up as negative.
fit$signed_R2 <- -sign(fit$slope) * fit$SFT.R.sq


fit[, c("Power", "signed_R2", "slope", "mean.k.", "median.k.", "max.k.")]

#Pick the lowest power that clears the threshold
r2_cut <- 0.85
passing <- fit$Power[fit$signed_R2 >= r2_cut]
power <- if (length(passing) > 0) min(passing) else NA
power                                              # NA means nothing reached the cutoff

## Plot both diagnostics
pdf("threshold.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot(fit$Power, fit$signed_R2, type = "n",
     xlab = "Soft threshold (power)", ylab = "Scale-free fit (signed R^2)")
text(fit$Power, fit$signed_R2, labels = fit$Power, col = "red")
abline(h = r2_cut, col = "red")

plot(fit$Power, fit$mean.k., type = "n",
     xlab = "Soft threshold (power)", ylab = "Mean connectivity")
text(fit$Power, fit$mean.k., labels = fit$Power, col = "red")




cor <- WGCNA::cor
net <- blockwiseModules(dat, power = power, networkType = "signed", TOMType = "signed",
                        minModuleSize = MIN_MOD, mergeCutHeight = MERGE_H,
                        deepsplit=0,
                        numericLabels = TRUE, maxBlockSize = 6000,
                        nThreads = nthr, verbose = 3)
cor <- stats::cor


# Grey fraction and stats
table(net$colors)                               # module sizes 0 = grey (unassigned)
mean(net$colors == 0)                           # grey fraction
length(unique(net$colors)) - 1                  # number of real modules

MEs <- net$MEs                                  # module eigengenes 468 x n_modules
mod_cols <- WGCNA::labels2colors(net$colors)    # numeric label -> colour name

# Plot the dendogram:
mod_cols <- WGCNA::labels2colors(net$colors)

pdf("gene_dendrogram.pdf", width = 12, height = 6)
plotDendroAndColors(net$dendrograms[[1]], mod_cols[net$blockGenes[[1]]],
                    "Module", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = paste("Gene dendrogram, power", power))
dev.off()


MElist <- moduleEigengenes(dat, colors = mod_cols)
ME.dissimilarity <- 1 - cor(MElist$eigengenes, use = "complete")
METree <- hclust(as.dist(ME.dissimilarity), method = "average")
pdf("eigengene_tree.pdf", width = 12, height = 6)
plot(METree, main = "Module eigengene clustering")
abline(h = 0.25, col = "red")
dev.off()

MEs_wo_grey <- MEs[, colnames(MEs) != "ME0"]

kME  <- cor(dat, MEs_wo_grey)                      # genes x modules (MEs without grey)
best <- colnames(kME)[max.col(kME)]
bestr <- apply(kME, 1, max)

reassign <- net$colors == 0 & bestr > 0.5
sum(reassign)                              # how many grey genes would move
table(best[reassign])                      # where they would go

MEs <- MEs_wo_grey
stopifnot(nrow(MEs) == length(gg))

# Kruskal-Wallis per module, FDR across modules
kw  <- sapply(MEs, function(me) kruskal.test(me ~ gg)$p.value)
fdr <- p.adjust(kw, method = "BH")

# Box colour = module colour (column names like "MEturquoise" -> "turquoise")
mod_names <- sub("^ME", "", colnames(MEs))
# if your columns are numeric (ME1, ME2...), use this instead:
mod_names <- WGCNA::labels2colors(as.integer(mod_names))

n <- ncol(MEs)
pdf("ME_by_grade_group.pdf", width = 12, height = 3 * ceiling(n / 4))
par(mfrow = c(ceiling(n / 4), 4), mar = c(3, 4, 3, 1))
for (i in seq_len(n)) {
  boxplot(MEs[[i]] ~ gg, col = mod_names[i], las = 1,
          xlab = "", ylab = "Eigengene",
          main = sprintf("%s\nFDR = %.2g", colnames(MEs)[i], fdr[i]))
  abline(h = 0, lty = 3, col = "grey50")
}
dev.off()

res <- data.frame(module = colnames(MEs), p = kw, fdr = fdr)
res[order(res$fdr), ]
