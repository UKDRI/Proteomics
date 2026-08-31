library(readxl)
library(readr)
library(reshape2)
library(ggplot2)
library(limma)
library(pheatmap)
library(mice)
library(preprocessCore)
library(corrplot)
library(odbc)
library(DBI)
library(RPostgres)
library(dplyr)
library(Rtsne)
library(rrcovNA)
library(future)
library(furrr)

set.seed(500)

beth.colours.reverse <- colorRampPalette(c("yellow", "darkgrey", "blue"), space = "rgb")(100)
beth.colours.colour  <- colorRampPalette(c("#D52D00", "#EF7627", "#FF9A56", "white", "#D162A4", "#B55690", "#A30262"), space = "rgb")(100)
source("~/Scripts/BETH_functions.r")

# -------------------------------------------------------------
# Data Ingestion & Preprocessing
# -------------------------------------------------------------
Aonly <- read_delim("report.pg_matrix.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Clean column headers
colnames(Aonly) <- gsub(".*DRI", "DRI", colnames(Aonly))
colnames(Aonly) <- gsub("_S2.*", "", colnames(Aonly))
colnames(Aonly) <- gsub(".raw", "", colnames(Aonly))
colnames(Aonly) <- gsub(".202.*", "", colnames(Aonly))
colnames(Aonly) <- gsub("_2026.*", "", colnames(Aonly))

m.Aonly <- data.matrix(Aonly[, -c(1:6)])
rownames(m.Aonly) <- Aonly$Protein.Group

# Mapping metadata
Condition <- mapper$Condition[match(colnames(m.Aonly), mapper$Raw)]
colnames(m.Aonly) <- mapper$Unique[match(colnames(m.Aonly), mapper$Raw)]
Outlier <- mapper$Outlier[match(colnames(m.Aonly), mapper$Unique)]

m.Aonly[m.Aonly == 0] <- NA

# Exclude outliers
m.Aonly <- m.Aonly[, grep(1, Outlier, invert = TRUE), drop = FALSE]
Condition <- Condition[grep(1, Outlier, invert = TRUE)]

saved.con <- Condition
saved.mAonly <- m.Aonly

# -------------------------------------------------------------
# Parallel Processing Setup
# -------------------------------------------------------------
plan(multisession, workers = availableCores() - 1)
mainDir <- getwd()

# Non-redundant pairwise condition list
cond_pairs <- combn(unique(saved.con), 2, simplify = FALSE)

future_walk(cond_pairs, function(pair) {
  
  cond1 <- pair[1]
  cond2 <- pair[2]
  pair_name <- paste(cond1, cond2, sep = "_vs_")
  
  sub_dir <- file.path(mainDir, pair_name)
  if (!dir.exists(sub_dir)) dir.create(sub_dir, recursive = TRUE)
  
  # Filter data for current pair
  pair_mask <- saved.con %in% c(cond1, cond2)
  m_sub <- saved.mAonly[, pair_mask, drop = FALSE]
  cond_sub <- saved.con[pair_mask]
  
  # Drop empty rows
  m_sub <- m_sub[rowSums(is.na(m_sub)) != ncol(m_sub), , drop = FALSE]
  
  # -------------------------------------------------------------
  # 1. QC Plots
  # -------------------------------------------------------------
  na_count <- colSums(is.na(m_sub))
  na.df <- data.frame(
    ProteinIDs = nrow(m_sub) - na_count,
    Sample = colnames(m_sub),
    stringsAsFactors = FALSE
  )
  
  raw_df <- as.data.frame(m_sub)
  raw_df$Protein.Group <- rownames(m_sub)
  
  me.data <- melt(raw_df, id.vars = "Protein.Group", variable.name = "Sample", value.name = "value")
  me.data <- me.data[complete.cases(me.data), ]
  me.data$Condition <- gsub(".*_", "", me.data$Sample)
  
  # Proteins per fraction plot
  pdf(file.path(sub_dir, "proteins_per_fraction.pdf"))
  print(ggplot(me.data, aes(Sample, fill = Condition)) +
          geom_bar() +
          ylab("Proteins per fraction") +
          ggtitle("Amount of proteins in each sample") +
          xlab("Fraction name") +
          theme_bw() + 
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)))
  dev.off()
  
  # Identification overlap plot
  rownas <- data.frame(Na.Count = rowSums(is.na(m_sub)), Protein.Group = rownames(m_sub))
  rownas$IDin <- ncol(m_sub) - rownas$Na.Count
  
  pdf(file.path(sub_dir, "ID_overlap.pdf"))
  print(ggplot(rownas, aes(IDin)) +
          geom_bar(fill = "#ffe119") +
          ylab("Number of proteins") +
          xlab("Identified in number of samples") +
          ggtitle("Protein identifications overlap") +
          theme_bw())
  dev.off()
  
  # Missing values density plot
  me.data$missing <- rownas$Na.Count[match(me.data$Protein.Group, rownas$Protein.Group)] > 1
  me.data$log2_val <- log2(me.data$value)
  
  pdf(file.path(sub_dir, "Missing_values_density.pdf"))
  print(ggplot(me.data, aes(log2_val, col = missing)) +
          geom_density() +
          scale_color_manual(values = c('#4363d8', '#f58231', '#dcbeff', '#800000', '#000075', '#a9a9a9', '#000000')) +
          xlab("Log2 Intensity") + 
          ylab("Density") + 
          theme_bw())
  dev.off()
  
  # Sample Density plot
  pdf(file.path(sub_dir, "density_plots.pdf"))
  print(ggplot(me.data, aes(log2_val, col = Sample)) +
          geom_density(show.legend = FALSE) +
          xlab("Log2 Intensity") + 
          facet_wrap(vars(Sample)) + 
          ylab("Density") + 
          theme_bw())
  dev.off()
  
  # -------------------------------------------------------------
  # 2. Normalisation
  # -------------------------------------------------------------
  norm_mat <- normalizeMedianValues(m_sub)
  rownames(norm_mat) <- rownames(m_sub)
  colnames(norm_mat) <- colnames(m_sub)
  
  norm.dat <- as.data.frame(norm_mat)
  norm.dat$Protein.Group <- rownames(norm_mat)
  
  me.nor.dat  <- melt(norm.dat, id.vars = "Protein.Group", variable.name = "Sample", value.name = "value")
  me.data_raw <- melt(raw_df, id.vars = "Protein.Group", variable.name = "Sample", value.name = "value")
  
  me.nor.dat$Norm  <- "Normalised"
  me.data_raw$Norm <- "Not Normalised"
  
  comb.me.data <- rbind(me.nor.dat, me.data_raw)
  comb.me.data$value <- log2(comb.me.data$value)
  
  pdf(file.path(sub_dir, "Normalisation_plot.pdf"))
  print(ggplot(comb.me.data, aes(value, Sample, fill = Norm)) +
          facet_grid(cols = vars(Norm)) +
          scale_fill_manual(values = c('#4363d8', '#f58231', '#dcbeff', '#800000', '#000075', '#a9a9a9', '#000000')) +
          xlab("Log2 Intensity") +
          ylab("") +
          geom_boxplot() + 
          theme_bw())
  dev.off()
  
  # Update m_sub with median-normalised matrix
  m_sub <- norm_mat
  
  # -------------------------------------------------------------
  # 3. Robust Imputation per Condition Split
  # -------------------------------------------------------------
  zeroProp <- 0.35
  
  impute_group_safe <- function(sub_mat, target_prop) {
    n_samples <- ncol(sub_mat)
    zero_prop <- apply(sub_mat, 1, function(r) sum(is.na(r) | r == 0) / n_samples)
    eligible_mask <- zero_prop <= target_prop
    
    out_mat <- sub_mat
    eligible_rows <- which(eligible_mask)
    
    if (length(eligible_rows) > 0) {
      for (r in eligible_rows) {
        row_vals <- out_mat[r, ]
        na_idx   <- is.na(row_vals) | row_vals == 0
        if (any(na_idx)) {
          med_val <- median(row_vals[!na_idx], na.rm = TRUE)
          if (is.na(med_val) || med_val <= 0) med_val <- 1
          out_mat[r, na_idx] <- med_val
        }
      }
    }
    return(out_mat)
  }
  
  s1_cols <- colnames(m_sub)[cond_sub == cond1]
  s2_cols <- colnames(m_sub)[cond_sub == cond2]
  
  imp_s1 <- impute_group_safe(m_sub[, s1_cols, drop = FALSE], zeroProp)
  imp_s2 <- impute_group_safe(m_sub[, s2_cols, drop = FALSE], zeroProp)
  
  comb.AO <- cbind(imp_s1, imp_s2)
  comb.AO <- comb.AO[order(rownames(comb.AO)), ]
  
  # Maintain exact metadata ordering matching comb.AO column structure
  ordered_cols <- colnames(comb.AO)
  cond_ordered <- saved.con[match(ordered_cols, colnames(saved.mAonly))]
  
  comb.AO[comb.AO < 0] <- 1
  comb.AO[is.na(comb.AO)] <- 1
  
  pdf(file.path(sub_dir, paste0(pair_name, "_cor_plots.pdf")))
  corrplot(cor(comb.AO, method = "spearman"), method = 'color', order = 'AOE')
  dev.off()
  
  # Export processed matrix
  write.csv(comb.AO, file = file.path(sub_dir, paste0(pair_name, "_processed_normalised_filtered.csv")), quote = FALSE)
  
  # -------------------------------------------------------------
  # 4. Statistical Comparisons & Limma Loops (Unadjusted with df check)
  # -------------------------------------------------------------
  m_log <- log2(comb.AO)
  zeroProps <- c(0, 0.5, 1)
  pros <- c("100", "50", "0")
  
  for (b in 1:3) {
    z_prop <- zeroProps[b]
    n_samp <- ncol(m_log)
    
    v_zero <- apply(m_log, 1, function(col) sum(col == 0 | is.na(col)) / n_samp)
    m_limma <- m_log[v_zero <= z_prop, , drop = FALSE]
    
    if (nrow(m_limma) == 0 || ncol(m_limma) < 3) next
    
    m_limma[is.na(m_limma) | m_limma == -Inf] <- 1
    m_limma <- m_limma[rowSums(m_limma) != ncol(m_limma), , drop = FALSE]
    
    if (nrow(m_limma) == 0 || ncol(m_limma) < 3) next
    
    # Set up single-factor model matrix
    new_cond <- factor(ifelse(cond_ordered == cond1, "Control", "Case"), levels = c("Control", "Case"))
    
    if (length(unique(new_cond)) < 2) next
    
    design <- model.matrix(~ new_cond)
    
    # Ensure residual degrees of freedom > 0 (samples > parameters)
    if (ncol(m_limma) <= ncol(design)) next
    
    fit  <- lmFit(m_limma, design)
    efit <- eBayes(fit)
    
    res <- topTreat(efit, coef = "new_condCase", adjust = "fdr", sort.by = "p", number = 50000)
    res <- as.data.frame(res)
    res <- res[order(res$logFC, decreasing = TRUE), ]
    res$Comparison <- pair_name
    
    write.csv(res, file = file.path(sub_dir, paste0(pair_name, "_", pros[b], "pct_Limma_Unadjusted.csv")), quote = FALSE)
    
    # -------------------------------------------------------------
    # Heatmaps & PCA Diagnostics (Safe Ellipse Evaluation)
    # -------------------------------------------------------------
    pdf(file.path(sub_dir, paste0(pair_name, "_", pros[b], "pct_HeatmapsPCA.pdf")))
    
    if (ncol(m_limma) >= 2 && nrow(m_limma) >= 2) {
      pca_res <- prcomp(t(m_limma), scale = FALSE)
      
      tryCatch({
        print(ggbiplot(pca_res, scale = 1, obs.scale = 0, var.scale = 0, 
                       var.axes = FALSE, varname.size = 0.75, ellipse.prob = 0.45, 
                       groups = cond_ordered, ellipse = TRUE, circle = FALSE))
      }, error = function(e) {
        print(ggbiplot(pca_res, scale = 1, obs.scale = 0, var.scale = 0, 
                       var.axes = FALSE, varname.size = 0.75, 
                       groups = cond_ordered, ellipse = FALSE, circle = FALSE))
      })
      
      tryCatch({
        print(ggbiplot(pca_res, scale = 1, obs.scale = 0, var.scale = 0, 
                       labels = colnames(m_limma), labels.size = 2, var.axes = FALSE, 
                       varname.size = 0.75, ellipse.prob = 0.45, groups = cond_ordered, 
                       ellipse = TRUE, circle = FALSE))
      }, error = function(e) {
        print(ggbiplot(pca_res, scale = 1, obs.scale = 0, var.scale = 0, 
                       labels = colnames(m_limma), labels.size = 2, var.axes = FALSE, 
                       varname.size = 0.75, groups = cond_ordered, 
                       ellipse = FALSE, circle = FALSE))
      })
      
      ann_col <- data.frame(
        Condition = factor(cond_ordered)
      )
      rownames(ann_col) <- colnames(m_limma)
      
      u_cond <- unique(as.character(ann_col$Condition))
      cond_palette <- c("#999999", "#E69F00")[1:length(u_cond)]
      names(cond_palette) <- u_cond
      
      ann_colors <- list(
        Condition = cond_palette
      )
      
      pheatmap(m_limma, annotation_colors = ann_colors, color = beth.colours.reverse, 
               show_colnames = TRUE, show_rownames = FALSE, annotation_col = ann_col, scale = "row")
    }
    dev.off()
  }
}, .options = furrr_options(seed = TRUE))

plan(sequential)