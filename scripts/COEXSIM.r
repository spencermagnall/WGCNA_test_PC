

normal_mods <- readRDS("../networks/Normal_modules.rds")
gg1_mods    <- readRDS("../networks/GG1_modules.rds")

normal_tom <- readRDS("../networks/Normal_tom.rds")$TOM
gg1_tom <- readRDS("../networks/GG1_tom.rds")$TOM

get_genes_per_mod <- function(mod){
    genes = mod$genes
    gene_colours = mod$mergedColors
    mod_colours = unique(gene_colours)
    mod_colours = mod_colours[mod_colours != "grey"]

    gene_per_mod <- list()

    for (m in mod_colours){
        gene_per_mod[[m]] = genes[gene_colours == m]
    }
    gene_per_mod
}



mean_edge_weight <- function(tom_matrix, gene_index) {
  sub_matrix <- tom_matrix[gene_index, gene_index]
  edges      <- sub_matrix[upper.tri(sub_matrix)]
  mean(edges)
}

# coexsim_pair <- function(genes_a, genes_b, tom_a, tom_b, all_genes_a, all_genes_b) {

#   shared   <- intersect(genes_a, genes_b)
#   n_shared <- length(shared)
#   size_rel <- n_shared / sqrt(length(genes_a) * length(genes_b))

#   # need at least 2 shared genes to have any edge
#   if (n_shared < 2) {
#     return(data.frame(n_shared = n_shared, size_rel = size_rel,
#                       sig_coex = NA, coexsim = 0))
#   }

#   w_a <- mean_edge_weight(tom_a, match(genes_a, all_genes_a))
#   w_b <- mean_edge_weight(tom_b, match(genes_b, all_genes_b))

#   w_shared_in_a <- mean_edge_weight(tom_a, match(shared, all_genes_a))
#   w_shared_in_b <- mean_edge_weight(tom_b, match(shared, all_genes_b))
#   w_shared      <- (w_shared_in_a + w_shared_in_b) / 2

#   sig_coex <- w_shared / (sqrt(w_a) * sqrt(w_b))
#   coexsim  <- size_rel * sig_coex

#   data.frame(n_shared = n_shared, size_rel = size_rel,
#              sig_coex = sig_coex, coexsim = coexsim)
# }

coexsim_all_pairs <- function(normal_gene_per_mod, tumour_gene_per_mod,
normal_tom, tumour_tom, normal_genes, tumour_genes){

   tumour_pairs <- NULL

  for (n in names(normal_gene_per_mod)) {
    for (t in names(tumour_gene_per_mod)) {
      row <- coexsim_pair(genes_a     = normal_gene_per_mod[[n]],
                          genes_b     = tumour_gene_per_mod[[t]],
                          tom_a       = normal_tom,
                          tom_b       = tumour_tom,
                          all_genes_a = normal_genes,
                          all_genes_b = tumour_genes)
      row$normal_module <- n
      row$tumour_module <- t
      tumour_pairs <- rbind(tumour_pairs, row)
    }
  }

  return(tumour_pairs)
}


# best_match_per_normal_module <- function(pairs) {
#   best_rows <- NULL

#   for (n in unique(pairs$normal_module)) {
#     rows_for_n <- pairs[pairs$normal_module == n, ]
#     best_index <- which.max(rows_for_n$coexsim)
#     best_rows  <- rbind(best_rows, rows_for_n[best_index, ])
#   }

#   return(best_rows)
# }


coexsim_pair <- function(genes_a, genes_b, tom_a, tom_b, all_genes_a, all_genes_b) {

  shared   <- intersect(genes_a, genes_b)
  n_shared <- length(shared)

  # gene-based (matches the paper's text)
  size_rel <- n_shared / sqrt(length(genes_a) * length(genes_b))

  # edge-based (matches the pasted code), assuming every pair of genes is an edge
  edges_a          <- choose(length(genes_a), 2)
  edges_b          <- choose(length(genes_b), 2)
  edges_shared     <- choose(n_shared, 2)
  size_rel_edges   <- edges_shared / sqrt(edges_a * edges_b)

  # need at least 2 shared genes to have any edge
  if (n_shared < 2) {
    return(data.frame(n_shared = n_shared,
                      size_rel = size_rel, size_rel_edges = size_rel_edges,
                      sig_coex = NA,
                      coexsim = 0, coexsim_edges = 0))
  }

  w_a <- mean_edge_weight(tom_a, match(genes_a, all_genes_a))
  w_b <- mean_edge_weight(tom_b, match(genes_b, all_genes_b))

  w_shared_in_a <- mean_edge_weight(tom_a, match(shared, all_genes_a))
  w_shared_in_b <- mean_edge_weight(tom_b, match(shared, all_genes_b))
  w_shared      <- (w_shared_in_a + w_shared_in_b) / 2

  sig_coex <- w_shared / (sqrt(w_a) * sqrt(w_b))

  data.frame(n_shared = n_shared,
             size_rel = size_rel, size_rel_edges = size_rel_edges,
             sig_coex = sig_coex,
             coexsim = size_rel * sig_coex,
             coexsim_edges = size_rel_edges * sig_coex)
}

best_match_per_normal_module <- function(pairs, score_column = "coexsim") {
  best_rows <- NULL

  for (n in unique(pairs$normal_module)) {
    rows_for_n <- pairs[pairs$normal_module == n, ]
    best_index <- which.max(rows_for_n[[score_column]])
    best_rows  <- rbind(best_rows, rows_for_n[best_index, ])
  }

  return(best_rows)
}


normal_gene_per_mod= get_genes_per_mod(normal_mods)
gg1_gene_per_mod = get_genes_per_mod(gg1_mods)

length(gg1_gene_per_mod)

# Testing modules 
genes_a = normal_gene_per_mod[["blue"]]
genes_b = gg1_gene_per_mod[["turquoise"]]
genes_a <- normal_gene_per_mod[["lightcyan"]]
genes_b <- gg1_gene_per_mod[["lightcyan"]]
genes_b

shared   <- intersect(genes_a, genes_b)
n_shared <- length(shared)

size_rel <- n_shared / sqrt(length(genes_a) * length(genes_b))
size_rel

w_a <- mean_edge_weight(normal_tom, match(genes_a, normal_mods$genes))
w_b <- mean_edge_weight(gg1_tom,    match(genes_b, gg1_mods$genes))

w_shared_in_a <- mean_edge_weight(normal_tom, match(shared, normal_mods$genes))
w_shared_in_b <- mean_edge_weight(gg1_tom,    match(shared, gg1_mods$genes))
w_shared      <- (w_shared_in_a + w_shared_in_b) / 2

sig_coex <- w_shared / (sqrt(w_a) * sqrt(w_b))
sig_coex
coexsim  <- size_rel * sig_coex
coexsim

# coexsim_pair(normal_gene_per_mod[["blue"]], gg1_gene_per_mod[["turquoise"]],
#              normal_tom, normal_mods$genes, gg1_tom, gg1_mods$genes)

# coexsim_pair(normal_gene_per_mod[["lightcyan"]], gg1_gene_per_mod[["lightcyan"]],
#              normal_tom, normal_mods$genes, gg1_tom, gg1_mods$genes)

# mod_pairs <- coexsim_all_pairs(normal_gene_per_mod, gg1_gene_per_mod,
# normal_tom, gg1_tom, normal_mods$genes, tumour_mods$genes)
# mod_pairs
# mod_pairs <- mod_pairs[order(-mod_pairs$coexsim), ]
# head(mod_pairs, 15)


#all_best <- NULL
all_best_edges <- NULL
all_best_nodes <- NULL

for (g in paste0("GG", 1:5)) {
  tumour_mods <- readRDS(file.path("../networks_plate_corrected", paste0(g, "_modules.rds")))
  tumour_tom  <- readRDS(file.path("../networks_plate_corrected", paste0(g, "_tom.rds")))$TOM
  tumour_gene_per_mod <- get_genes_per_mod(tumour_mods)

  pairs <- coexsim_all_pairs(normal_gene_per_mod, tumour_gene_per_mod,
                             normal_tom, tumour_tom,
                             normal_mods$genes, tumour_mods$genes)

  #best <- best_match_per_normal_module(pairs)
  best_nodes <- best_match_per_normal_module(pairs, "coexsim")
  best_edges <- best_match_per_normal_module(pairs, "coexsim_edges")
  best_nodes$group <- g
  best_edges$group <- g
  #best$group <- g
  #all_best <- rbind(all_best, best)
  all_best_edges <- rbind(all_best_edges, best_edges)
  all_best_nodes <- rbind(all_best_nodes, best_nodes)
}
#all_best

tapply(all_best_nodes$coexsim       > 0.25, all_best_nodes$group, sum)
tapply(all_best_edges$coexsim_edges > 0.25, all_best_edges$group, sum)

worst_nodes <- tapply(all_best_nodes$coexsim,       all_best_nodes$normal_module, min)
worst_edges <- tapply(all_best_edges$coexsim_edges, all_best_edges$normal_module, min)
sum(worst_nodes > 0.25)
sum(worst_edges > 0.25)

sort(worst_edges, decreasing = TRUE)[1:8]

# core <- c("coral2", "skyblue", "lightgreen", "blue")

# for (g in paste0("GG", 1:5)) {
#   m <- readRDS(file.path("../networks", paste0(g, "_modules.rds")))
#   for (mod in core) {
#     present <- mean(normal_mods[[mod]] %in% m$genes)
#     cat(g, mod, "fraction of genes in network:", round(present, 2), "\n")
#   }
# }

# all_best_edges[all_best_edges$group == "GG4" &
#                all_best_edges$normal_module %in% core, ]

# normal_sets <- get_genes_per_mod(normal_mods)
# length(normal_sets)              # expect 40
# length(normal_sets[["coral2"]])  # expect 46
# names(normal_sets)[1:5]

# worst_edges <- tapply(all_best_edges$coexsim_edges, all_best_edges$normal_module, min)
# sort(worst_edges, decreasing = TRUE)[1:8]



# # How do our core modules match across GGs?
# core <- c("darkslateblue", "white", "lightyellow", "bisque4")

# core_matches <- all_best_edges[all_best_edges$normal_module %in% core,
#                                c("normal_module", "group", "tumour_module",
#                                  "n_shared", "coexsim_edges")]

# core_matches <- core_matches[order(core_matches$normal_module, core_matches$group), ]
# core_matches$coexsim_edges <- round(core_matches$coexsim_edges, 2)
# core_matches

# worst_score <- tapply(all_best$coexsim, all_best$normal_module, min)
# sort(worst_score, decreasing = TRUE)
# tapply(all_best$coexsim, all_best$group, max)
# tapply(all_best$n_shared, all_best$group, max)

# for (g in paste0("GG", 1:5)) {
#   m    <- readRDS(file.path("../networks", paste0(g, "_modules.rds")))
#   sets <- get_genes_per_mod(m)
#   cat(g,
#       " genes:", length(m$genes),
#       " shared with Normal:", length(intersect(normal_mods$genes, m$genes)),
#       " modules:", length(sets),
#       " largest module:", max(lengths(sets)), "\n")
# }

# tumour_tom <- readRDS(file.path("../networks", "GG4_tom.rds"))$TOM
# dim(tumour_tom)
# length(readRDS(file.path("../networks", "GG4_modules.rds"))$genes)
# all_best[all_best$normal_module == "lightcyan", ]
# tapply(all_best$coexsim > 0.25, all_best$group, sum)


# worst_score <- tapply(all_best$coexsim, all_best$normal_module, min)
# sort(worst_score, decreasing = TRUE)


# tapply(all_best$coexsim  > 0.25, all_best$group, sum)
# tapply(all_best$size_rel > 0.25, all_best$group, sum)

# for (g in paste0("GG", 1:5)) {
#   m <- readRDS(file.path("../networks", paste0(g, "_modules.rds")))
#   cat(g, "grey fraction:", round(mean(m$mergedColors == "grey"), 3), "\n")
# }

# gg3_mods <- readRDS("../networks/GG3_modules.rds")
# gg4_mods <- readRDS("../networks/GG4_modules.rds")
# gg3_tom  <- readRDS("../networks/GG3_tom.rds")$TOM
# gg4_tom  <- readRDS("../networks/GG4_tom.rds")$TOM

# pairs_34 <- coexsim_all_pairs(get_genes_per_mod(gg3_mods), get_genes_per_mod(gg4_mods),
#                               gg3_tom, gg4_tom,
#                               gg3_mods$genes, gg4_mods$genes)

# best_34 <- best_match_per_normal_module(pairs_34)
# sum(best_34$coexsim > 0.25)
# nrow(best_34)

# sets3 <- get_genes_per_mod(gg3_mods)
# sets4 <- get_genes_per_mod(gg4_mods)
# length(sets3)     # expect 21
# length(sets4)     # expect 19

# big3 <- names(which.max(lengths(sets3)))
# big4 <- names(which.max(lengths(sets4)))

# observed <- length(intersect(sets3[[big3]], sets4[[big4]]))
# expected <- length(sets3[[big3]]) * length(sets4[[big4]]) /
#             length(intersect(gg3_mods$genes, gg4_mods$genes))
# c(observed = observed, expected = expected)

# coexsim_pair(genes_a = sets3[[big3]], genes_b = sets4[[big4]],
#              tom_a = gg3_tom, all_genes_a = gg3_mods$genes,
#              tom_b = gg4_tom, all_genes_b = gg4_mods$genes)

# nrow(pairs_34)                                  # expect 21 * 19 = 399
# length(unique(pairs_34$normal_module))          # expect 21
# pairs_34[pairs_34$normal_module == big3 & pairs_34$tumour_module == big4, ]

# sum(is.na(gg3_tom))
# sum(is.na(gg4_tom))

# bad_genes_gg4 <- which(colSums(is.na(gg4_tom)) > 0)
# length(bad_genes_gg4)
# head(gg4_mods$genes[bad_genes_gg4])


# a <- sets3[[big3]]
# b <- sets4[[big4]]
# shared <- intersect(a, b)

# sum(is.na(match(a, gg3_mods$genes)))        # expect 0
# sum(is.na(match(b, gg4_mods$genes)))        # expect 0
# sum(is.na(match(shared, gg3_mods$genes)))   # expect 0
# sum(is.na(match(shared, gg4_mods$genes)))   # expect 0

# mean_edge_weight(gg3_tom, match(a, gg3_mods$genes))
# mean_edge_weight(gg4_tom, match(b, gg4_mods$genes))

# codetools::findGlobals(coexsim_pair)
# codetools::findGlobals(coexsim_all_pairs)