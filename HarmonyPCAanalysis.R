#!/usr/bin/env Rscript
# =============================================================================
# Neighbourhood-matched comparison of clonal histories within the My-HSC state
# Stage 1 of 2 : analysis  ->  writes the tables that pipeline_02_figures.py plots
#
# Inputs
#   C02_2.rds  Seurat object, 11,717 cells; meta.data must carry
#              Sample, BC (clone barcode, "N" = unbarcoded), Cluster (c1..c4, "N")
#              reductions: harmony (>=30 dims), umap2
#   MolO.txt   tab-delimited, one row per cell IN THE SAME ORDER as the object,
#              columns include MolO plus nCount_RNA / BC / Cluster / Sample for checking
#
# Label convention (fixed by the study)
#   c1 = P-HSC, c2 = PRG-HSC, c3 = PRGB-HSC   -> collectively My-HSCs
#   c4 = balanced-HSC
#
# Outputs (out/)
#   panelA_nhood_composition.csv   2,781 neighbourhoods: DA stats + subclass fractions
#   panelB_clone_fractions.csv     per-clone % of cells inside the My-HSC-enriched region
#   panelCD_percell.csv            per-cell MolO + region membership + subclass + sample
#   panelE_percell_weights.csv     per-cell composition-adjustment weights
#   stats_*.csv                    the numbers quoted in the text
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(RANN); library(edgeR)
})

SEURAT_RDS <- Sys.getenv("SEURAT_RDS", "C02_2.rds")
MOLO_TXT   <- Sys.getenv("MOLO_TXT",   "MolO.txt")
OUT        <- "out"; dir.create(OUT, showWarnings = FALSE)

N_DIMS    <- 30    # harmony dimensions defining the neighbourhood space
K_NN      <- 50    # neighbours per neighbourhood (size = K_NN + 1 index cell)
IDX_PROP  <- 0.45  # fraction of cells sampled as candidate index cells
FDR_CUT   <- 0.05  # neighbourhood-level DA threshold
MIN_CELLS <- 10    # clones kept in panel B
SEED      <- 42

MY <- c("c1","c2","c3"); SUB <- c("c1","c2","c3","c4")

# ---- 1. load, and verify the MolO file is aligned cell-for-cell ---------------
obj <- readRDS(SEURAT_RDS)
md  <- obj@meta.data
CL  <- md$Cluster; BC <- md$BC; SAMPLE <- md$Sample

mol <- read.delim(MOLO_TXT, stringsAsFactors = FALSE)
stopifnot(nrow(mol) == ncol(obj))
stopifnot(all(mol$BC == BC), all(mol$Cluster == CL), all(mol$Sample == SAMPLE))
MOLO <- mol$MolO
cat(sprintf("loaded %d cells; MolO aligned (BC/Cluster/Sample all match)\n", ncol(obj)))

# Cluster is essentially a CLONE-level attribute. In this dataset 3 of 434 clones
# carry two labels (the minority side being 1 cell in 2 of them), so take the
# majority label per clone and record how often that had to be done.
ok     <- BC != "N" & CL != "N"
cl_raw <- table(BC[ok], CL[ok])
cl_tab <- data.frame(BC = rownames(cl_raw),
                     CL = colnames(cl_raw)[apply(cl_raw, 1, which.max)],
                     stringsAsFactors = FALSE)
n_mixed <- sum(rowSums(cl_raw > 0) > 1)
cat(sprintf("clones: %d | %d carried >1 label -> majority label used\n",
            nrow(cl_tab), n_mixed))
CL <- ifelse(ok, cl_tab$CL[match(BC, cl_tab$BC)], CL)   # harmonise cell labels

# ---- 2. clone-blind neighbourhoods on the harmony embedding ------------------
# Milo-style: sample candidate index cells, refine each to the medoid of its
# own kNN, keep the unique set. Neighbourhoods OVERLAP - this is not a clustering.
set.seed(SEED)
E  <- Embeddings(obj, "harmony")[, 1:N_DIMS]
Nc <- nrow(E)
nn.idx <- nn2(E, k = K_NN + 1)$nn.idx

cand   <- sample(seq_len(Nc), size = floor(IDX_PROP * Nc))
refine <- function(i) {
  ix  <- nn.idx[i, ]
  med <- colMeans(E[ix, , drop = FALSE])
  ix[which.min(colSums((t(E[ix, , drop = FALSE]) - med)^2))]
}
idx <- unique(vapply(cand, refine, integer(1)))

nh_list <- lapply(idx, function(i) nn.idx[i, ])
M <- sparseMatrix(i = unlist(nh_list),
                  j = rep(seq_along(nh_list), each = K_NN + 1),
                  x = 1, dims = c(Nc, length(nh_list)))
cat(sprintf("neighbourhoods: %d | size %d | cells covered %d / %d\n",
            ncol(M), K_NN + 1, sum(rowSums(M) > 0), Nc))

# ---- 3. differential abundance, clones as replicates -------------------------
# Features = neighbourhoods, columns = clones, offset = clone size, Sample as covariate.
lab_ok <- BC != "N" & CL %in% SUB
bcl    <- factor(BC[lab_ok])
bc_ind <- sparse.model.matrix(~ 0 + bcl)                 # cells x clones
cnt    <- as.matrix(t(M[lab_ok, , drop = FALSE]) %*% bc_ind)   # nhood x clone
colnames(cnt) <- levels(bcl)

clone_meta <- data.frame(BC = colnames(cnt), stringsAsFactors = FALSE)
clone_meta$CL     <- cl_tab$CL[match(clone_meta$BC, cl_tab$BC)]
clone_meta$hist   <- ifelse(clone_meta$CL %in% MY, "My", "Bal")
clone_meta$Sample <- vapply(clone_meta$BC, function(b)
  names(sort(table(SAMPLE[BC == b]), decreasing = TRUE))[1], character(1))
# A clone qualifies as a replicate on its CELL count, not on its footprint across
# the overlapping neighbourhoods (the latter counts the same cell many times and
# would admit ~2x as many clones).
clone_meta$n_cells <- as.integer(table(factor(BC[lab_ok], levels = colnames(cnt))))

keep_clone <- clone_meta$n_cells >= 10       # clones with enough cells to be a replicate
Y  <- cnt[, keep_clone, drop = FALSE]
cm <- clone_meta[keep_clone, ]
# carry the neighbourhood id explicitly: filtering rows below would otherwise
# leave the DA results joined to the wrong neighbourhoods
rownames(Y) <- as.character(seq_len(nrow(Y)))
Y  <- Y[rowSums(Y) > 0, , drop = FALSE]
cat(sprintf("DA design: %d neighbourhoods x %d clones (My %d / Bal %d)\n",
            nrow(Y), ncol(Y), sum(cm$hist == "My"), sum(cm$hist == "Bal")))

run_da <- function(Y, meta, grp) {
  d   <- DGEList(Y, lib.size = colSums(Y))
  d   <- calcNormFactors(d, method = "TMM")
  des <- model.matrix(~ Sample + grp, data.frame(Sample = meta$Sample, grp = grp))
  d   <- estimateDisp(d, des)
  fit <- glmQLFit(d, des, robust = TRUE)
  cf  <- tail(colnames(des), 1)
  tt  <- topTags(glmQLFTest(fit, coef = cf), n = Inf, sort.by = "none")$table
  tt$nhood <- as.integer(rownames(Y))
  tt
}
# positive logFC = enriched for My-HSC clones
tt <- run_da(Y, cm, factor(cm$hist, levels = c("Bal","My")))

sig_my  <- tt$nhood[tt$FDR < FDR_CUT & tt$logFC > 0]
sig_bal <- tt$nhood[tt$FDR < FDR_CUT & tt$logFC < 0]
cat(sprintf("My-HSC-enriched: %d | balanced-enriched: %d (FDR < %.2f)\n",
            length(sig_my), length(sig_bal), FDR_CUT))

# region = union of the significant My-enriched neighbourhoods,
# minus any cell that also sits in a balanced-enriched one
in_my  <- rowSums(M[, sig_my,  drop = FALSE]) > 0
in_bal <- rowSums(M[, sig_bal, drop = FALSE]) > 0
REGION <- ifelse(in_my & !in_bal, "MyRich",
          ifelse(in_bal & !in_my, "BalRich", "Other"))
cat(sprintf("region: MyRich %d cells (of which barcoded+labelled %d)\n",
            sum(REGION == "MyRich"), sum(REGION == "MyRich" & lab_ok)))

# ---- 4. panel A : composition along the DA axis ------------------------------
frac <- vapply(SUB, function(k)
  as.numeric((t(M) %*% (CL == k)) / pmax(as.numeric(t(M) %*% (CL %in% SUB)), 1)),
  numeric(ncol(M)))
colnames(frac) <- paste0("f_", SUB)

U2 <- Embeddings(obj, "umap2")
A <- data.frame(nhood = seq_len(ncol(M)),
                u2_1 = U2[idx, 1], u2_2 = U2[idx, 2],
                molo = as.numeric((t(M) %*% MOLO) / pmax(rowSums(t(M)), 1)),
                frac, check.names = FALSE)
A <- merge(A, tt[, c("nhood","logFC","F","PValue","FDR")], by = "nhood", all.x = TRUE)
A$sig  <- ifelse(A$nhood %in% sig_my, "My",
          ifelse(A$nhood %in% sig_bal, "Bal", "ns"))
A      <- A[order(A$logFC), ]
A$rank <- seq_len(nrow(A)) - 1L
write.csv(A, file.path(OUT, "panelA_nhood_composition.csv"), row.names = FALSE)

# ---- 5. panel B : per-clone occupancy of the region --------------------------
cb <- data.frame(BC = BC, CL = CL, REGION = REGION,
                 Sample = SAMPLE, stringsAsFactors = FALSE)
cb <- subset(cb, BC != "N" & CL %in% SUB)
B  <- aggregate(cbind(n_cells = rep(1, nrow(cb)),
                      n_in    = as.integer(cb$REGION == "MyRich")) ~ BC, cb, sum)
B$CL     <- cl_tab$CL[match(B$BC, cl_tab$BC)]
B$Sample <- clone_meta$Sample[match(B$BC, clone_meta$BC)]
B$pct    <- 100 * B$n_in / B$n_cells
write.csv(B, file.path(OUT, "panelB_clone_fractions.csv"), row.names = FALSE)

Bk <- subset(B, n_cells >= MIN_CELLS)
sB <- do.call(rbind, lapply(SUB, function(k) {
  v <- Bk$pct[Bk$CL == k]
  data.frame(Cluster = k, n_clones = length(v), mean = mean(v),
             sem = sd(v) / sqrt(length(v)))
}))
write.csv(sB, file.path(OUT, "stats_panelB.csv"), row.names = FALSE)
cat("\npanel B (clones with >=", MIN_CELLS, "cells):\n"); print(sB, row.names = FALSE, digits = 4)

# ---- 6. panels C/D : per-cell MolO ------------------------------------------
CD <- data.frame(cell = seq_len(Nc), Sample = SAMPLE, Cluster = CL, BC = BC,
                 region = REGION, MolO = MOLO,
                 inout = ifelse(REGION == "MyRich", "inside", "outside"),
                 u2_1 = U2[, 1], u2_2 = U2[, 2], stringsAsFactors = FALSE)
write.csv(CD, file.path(OUT, "panelCD_percell.csv"), row.names = FALSE)

sC <- do.call(rbind, lapply(c("inside","outside"), function(z) {
  v <- CD$MolO[CD$inout == z]
  data.frame(region = z, n = length(v), mean = mean(v),
             sem = sd(v) / sqrt(length(v)), sd = sd(v))
}))
sC$p_mannwhitney <- wilcox.test(CD$MolO[CD$inout == "inside"],
                                CD$MolO[CD$inout == "outside"])$p.value
write.csv(sC, file.path(OUT, "stats_panelC.csv"), row.names = FALSE)
cat("\npanel C (inside vs outside, cells as n):\n"); print(sC, row.names = FALSE, digits = 4)

inreg <- subset(CD, region == "MyRich" & Cluster %in% SUB)
sD <- do.call(rbind, lapply(SUB, function(k) {
  v <- inreg$MolO[inreg$Cluster == k]
  data.frame(Cluster = k, n = length(v), mean = mean(v),
             sem = sd(v) / sqrt(length(v)), sd = sd(v))
}))
pw <- t(combn(SUB, 2)); praw <- apply(pw, 1, function(p)
  wilcox.test(inreg$MolO[inreg$Cluster == p[1]],
              inreg$MolO[inreg$Cluster == p[2]])$p.value)
sDp <- data.frame(pair = paste(pw[,1], pw[,2], sep = "-"),
                  p = praw, q = p.adjust(praw, "BH"))
write.csv(sD,  file.path(OUT, "stats_panelD.csv"),      row.names = FALSE)
write.csv(sDp, file.path(OUT, "stats_panelD_pairs.csv"), row.names = FALSE)
cat("\npanel D (within region, cells as n):\n");  print(sD,  row.names = FALSE, digits = 4)
cat("pairwise (BH):\n");                          print(sDp, row.names = FALSE, digits = 3)

# ---- 7. panel E : composition-adjusted MolO ---------------------------------
# Each cell gets w = sum over its My-enriched neighbourhoods of n_c4(j)/n_k(j),
# so every subclass carries the same total weight as c4 within each neighbourhood.
# One row per CELL: neighbourhoods overlap, so per-(cell,neighbourhood) rows
# would count the same cell many times and shrink the SEM artificially.
w <- numeric(Nc)
for (j in sig_my) {
  mem <- which(M[, j] > 0); mem <- mem[CL[mem] %in% SUB]
  n4  <- sum(CL[mem] == "c4"); if (!length(mem) || n4 == 0) next
  for (k in SUB) {
    ik <- mem[CL[mem] == k]
    if (length(ik)) w[ik] <- w[ik] + n4 / length(ik)
  }
}
Ew <- data.frame(cell = which(w > 0 & CL %in% SUB))
Ew$Cluster <- CL[Ew$cell]; Ew$Sample <- SAMPLE[Ew$cell]
Ew$MolO    <- MOLO[Ew$cell]; Ew$w <- w[Ew$cell]
write.csv(Ew, file.path(OUT, "panelE_percell_weights.csv"), row.names = FALSE)

sE <- do.call(rbind, lapply(SUB, function(k) {
  g <- Ew[Ew$Cluster == k, ]
  m <- weighted.mean(g$MolO, g$w)
  v <- weighted.mean((g$MolO - m)^2, g$w)
  ne <- sum(g$w)^2 / sum(g$w^2)
  data.frame(Cluster = k, n_cells = nrow(g), sum_w = sum(g$w), n_eff = ne,
             wmean = m, wsd = sqrt(v), wsem = sqrt(v / ne))
}))
write.csv(sE, file.path(OUT, "stats_panelE.csv"), row.names = FALSE)
cat("\npanel E (composition-adjusted):\n"); print(sE, row.names = FALSE, digits = 4)

# clone-level label permutation: the null that respects the sampling unit
set.seed(SEED + 1)
wdiff <- function(g1, g2) weighted.mean(g1$MolO, g1$w) - weighted.mean(g2$MolO, g2$w)
perm <- do.call(rbind, lapply(MY, function(k) {
  s   <- Ew[Ew$Cluster %in% c(k, "c4"), ]
  obs <- wdiff(s[s$Cluster == k, ], s[s$Cluster == "c4", ])
  ct  <- unique(data.frame(BC = BC[s$cell], lab = s$Cluster, stringsAsFactors = FALSE))
  nul <- replicate(2000, {
    lt <- ct; lt$lab <- sample(lt$lab)
    nl <- lt$lab[match(BC[s$cell], lt$BC)]
    if (sum(nl == k) < 5 || sum(nl == "c4") < 5) return(NA_real_)
    wdiff(s[nl == k, ], s[nl == "c4", ])
  })
  data.frame(pair = paste0(k, "-c4"), diff = obs,
             null_mean = mean(nul, na.rm = TRUE), null_sd = sd(nul, na.rm = TRUE),
             p = mean(abs(nul) >= abs(obs), na.rm = TRUE))
}))
write.csv(perm, file.path(OUT, "stats_panelE_perm.csv"), row.names = FALSE)
cat("\nclone-permuted null vs c4:\n"); print(perm, row.names = FALSE, digits = 3)

# ---- 8. weighted differential expression, calibrated by permutation ---------
# The composition weights of section 7 are reused on the EXPRESSION matrix.
#
# WARNING, and the reason this section exists: handing these weights to a
# parametric fit (lmFit(..., weights = w)) inflates the gene count by more than
# an order of magnitude. Weighting changes the effective sample size but not the
# nominal residual degrees of freedom, so the variance is underestimated and the
# moderated t is too large. Running the SAME cells through the same fit with the
# weights removed gives 5-32 genes at FDR<0.01 (10-57 at FDR<0.05) where the
# weighted fit gives 79-1208 at FDR<0.01, across the six criterion x comparison
# combinations in weighted_DE_threshold_summary.csv.
# Never report the weighted parametric counts. Calibrate against a permuted null
# built with the same weighting instead, which is what is done below.
#
# Test statistic: weighted mean difference in log-normalised expression.
# Two nulls, because they answer different questions:
#   cell  permutation -> appropriate for an EXISTENCE claim ("some cells differ")
#   clone permutation -> appropriate for a NO-DIFFERENCE claim; it respects the
#                        sampling unit, and widens specifically where the signal
#                        is, because a whole clone moves as a block.
NPERM_DE <- 200
MIN_CELLS_DE <- 5

Xn  <- GetAssayData(obj, layer = "data")
sel <- Ew$cell                                   # cells carrying a weight
expressed <- Matrix::rowSums(Xn[, sel, drop = FALSE] > 0) >= 10
Xw  <- Xn[expressed, sel, drop = FALSE]
wv  <- Ew$w; clv <- Ew$Cluster; bcv <- BC[sel]
cat(sprintf("\nweighted DE: %d genes x %d cells | clones %d\n",
            nrow(Xw), ncol(Xw), length(unique(bcv[bcv != "N"]))))

# weighted mean difference, vectorised over all genes at once
wmd <- function(Y, ww, is_k) {
  wk <- ww[is_k]; w4 <- ww[!is_k]
  as.numeric((Y[, is_k, drop = FALSE] %*% wk) / sum(wk) -
             (Y[, !is_k, drop = FALSE] %*% w4) / sum(w4))
}

set.seed(SEED + 2)
de_all <- list()
for (k in MY) {
  keep <- clv %in% c("c4", k)
  Y <- Xw[, keep, drop = FALSE]; ww <- wv[keep]
  is_k <- clv[keep] == k; bck <- bcv[keep]
  obs  <- wmd(Y, ww, is_k)

  # null 1: shuffle the subclass label across CELLS (weights stay on the cell)
  nc <- replicate(NPERM_DE, wmd(Y, ww, sample(is_k)))

  # null 2: shuffle the label across CLONES, so every cell of a clone moves together
  ct <- unique(data.frame(BC = bck, lab = ifelse(is_k, k, "c4"),
                          stringsAsFactors = FALSE))
  nl <- replicate(NPERM_DE, {
    lt <- ct; lt$lab <- sample(lt$lab)
    lab <- lt$lab[match(bck, lt$BC)] == k
    if (sum(lab) < MIN_CELLS_DE || sum(!lab) < MIN_CELLS_DE) return(rep(NA_real_, nrow(Y)))
    wmd(Y, ww, lab)
  })

  sd_cell  <- apply(nc, 1, sd)
  sd_clone <- apply(nl, 1, sd, na.rm = TRUE)
  z_cell   <- obs / pmax(sd_cell,  1e-8)
  z_clone  <- obs / pmax(sd_clone, 1e-8)
  de_all[[k]] <- data.frame(
    comparison = paste0(k, "_vs_c4"), gene = rownames(Y), diff = obs,
    sd_cell = sd_cell,   z_cell  = z_cell,
    FDR_cell  = p.adjust(2 * pnorm(-abs(z_cell)),  "BH"),
    sd_clone = sd_clone, z_clone = z_clone,
    FDR_clone = p.adjust(2 * pnorm(-abs(z_clone)), "BH"),
    stringsAsFactors = FALSE)
  cat(sprintf("  %s vs c4: n = %d vs %d | FDR<0.05  cell-perm %4d  clone-perm %4d\n",
              k, sum(is_k), sum(!is_k),
              sum(de_all[[k]]$FDR_cell < 0.05), sum(de_all[[k]]$FDR_clone < 0.05)))
}
DE <- do.call(rbind, de_all)
write.csv(DE, file.path(OUT, "weighted_DE_percell.csv"), row.names = FALSE)

sDE <- do.call(rbind, lapply(MY, function(k) {
  d <- DE[DE$comparison == paste0(k, "_vs_c4"), ]
  data.frame(comparison = paste0(k, "_vs_c4"), n_genes = nrow(d),
             cell_perm_05  = sum(d$FDR_cell  < 0.05), cell_perm_01  = sum(d$FDR_cell  < 0.01),
             clone_perm_05 = sum(d$FDR_clone < 0.05), clone_perm_01 = sum(d$FDR_clone < 0.01),
             median_sd_ratio = median(d$sd_clone / d$sd_cell, na.rm = TRUE))
}))
write.csv(sDE, file.path(OUT, "stats_weighted_DE.csv"), row.names = FALSE)
cat("\nweighted DE summary:\n"); print(sDE, row.names = FALSE, digits = 4)

cat("\ndone -> ", normalizePath(OUT), "\n")
