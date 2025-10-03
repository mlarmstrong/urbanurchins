# Urban Urchins Window Association Tests
# Author: Katie Erickson
# Last Edit: September 19, 2025

# Running on a cluster with bcftools loaded as a dependancy
# Load packages and set parameters --------------------------------------------
packages <- c("dplyr", "data.table", "lostruct", "tidyr", "ggplot2", "RColorBrewer")
lapply(packages, library, character.only = TRUE)

window_size <- 100  # number of SNPs in each windows
k_kept <- 2         # number of PCs to keep
num_sims <- 100     # number of simulations to run

# set paths to VCF and metadata
vcf_path <- "/group/rbaygrp/katie/local_pca/uban_urchins/uban_urchins.vcf.gz"
meta_data_csv <- "/group/rbaygrp/katie/local_pca/uban_urchins/urbanurchins_metadata_thinned.csv"

# read in and reformat metadata,  (make sure variables of interest are factor)
meta_data <- read.csv(meta_data_csv)
meta_data$Dev <- factor(meta_data$Dev, levels = c("nonurban", "urban"))
meta_data$region <- factor(meta_data$region, levels = c("SD", "LA", "Vic"))
num_samples <- length(meta_data$Sample_ID)

#create positions matrix from vcf
vcf <- read.vcfR(vcf_path, verbose = FALSE)
vcf_fix <- as.data.frame(vcf@fix)
chroms <- vcf_fix$CHROM
snp_positions <- as.numeric(vcf_fix$POS)

# use lostruct functions to create windows and run PCAs -------------------
# create windows from vcf
sites <- vcf_positions(vcf_path)
win.fn.snp <- vcf_windower(vcf_path, size = window_size, type = "snp", sites = sites)
max_windows <- attr(win.fn.snp, "max.n")

# run PCA across all windows
all_windows_pca <- eigen_windows(win.fn.snp, k = k_kept)

# extract and reformat pc and sample names for ease of use in looping
pc_labels <- names(all_windows_pca[1,])
pc_labels <- grep("^PC_\\d+_", pc_labels, value = TRUE) #remove total lam_1, lam_2
pcs <- sort(unique(sub("^((PC_\\d+))_.*$", "\\1", pc_labels)))
sample_names <- sort(unique(sub("^PC_\\d+_", "", pc_labels)))
labs <- colnames(all_windows_pca)

# initialize variables used in running and storing model outputs
keep_terms    <- c("Dev", "region", "Dev:region")
keep_cols_in  <- c("Term", "Sum Sq", "Df", "F value", "Pr(>F)")
keep_cols_out <- c("Term", "SumSq", "Df", "F", "P")
pc_cols       <- paste0("PC_", seq_len(k_kept))

# Loop through real data --------------------------------------------------
all_pcs_df <- data.frame()
all_coefs  <- list()
idx_real <- 1
for (w in 1:max_windows) {
  # Assemble per-window PC frame
  window_data <- data.frame( Sample_ID = sample_names, window = w,
    start_bp = snp_positions[(w - 1) * window_size + 1], 
    end_bp = snp_positions[w * window_size],
    chr = chroms[(w - 1) * window_size + 1],
    stringsAsFactors = FALSE)
  window_data[pcs] <- NA_real_
  
  # Fill each PC column for this window
  for (pc in pcs) {
    pc_indices <- grep(paste0("^", pc, "_"), labs)
    if (length(pc_indices) == 0) next #check that there are PCs to move forward
    pc_values <- as.numeric(all_windows_pca[w, pc_indices])
    names(pc_values) <- sub("^PC_\\d+_", "", pc_labels[pc_indices])
    window_data[[pc]] <- pc_values[ match(window_data$Sample_ID, names(pc_values)) ]
  }
  
  # Append PCs for this window
  all_pcs_df <- rbind(all_pcs_df, window_data)
  
  # Merge metadata with window data to run model
  window_data <- merge(window_data, meta_data, by = "Sample_ID", all.x = TRUE)
  
  # for each pc
  for (k in 1:k_kept) {
    pc_col <- paste0("PC_", k)
    if (all(is.na(window_data[[pc_col]]))) next
    
    # Fit ANOVA and save relevant coefficients
    model <- lm(as.formula(paste(pc_col, "~ Dev * region")), data = window_data)
    a_tbl <- as.data.frame(Anova(model, type = "II"))
    a_tbl$Term <- rownames(a_tbl)
    
    rows <- a_tbl$Term %in% keep_terms
    if (length(keep_terms) == 0) next
    
    out <- a_tbl[rows, keep_cols_in, drop = FALSE]
    names(out) <- keep_cols_out
    
    # Add window data to model coefficients dataframe
    out$PC <- pc_col
    out$window <- w
    out$chr <- unique(window_data$chr)[1]
    out$start_bp <- unique(window_data$start_bp)[1]
    out$end_bp <- unique(window_data$end_bp)[1]
    
    # Add coefficent data to overall data frame and move to next index
    all_coefs[[idx_real]] <- out
    idx_real <- idx_real + 1
  }
}

#reformat output and write to csv
all_coefs <- dplyr::bind_rows(all_coefs) # Reformat list into data frame
real_df <- all_coefs %>% mutate(source = "Real")
write.csv(real_df, file ="/group/rbaygrp/katie/local_pca/uban_urchins/real_df_8_29.csv")

# Run simulations ---------------------------------------------------------
all_sim_coefs <- list()
idx_sim <- 1
for (sim in 1:num_sims) {
  cat("Starting simulation ", sim, " out of ", num_sims, "\n")
  
  # Randomize the metadata and merge it with the PC data
  set.seed(sim) #gives new seed for each loop while allowing reproducibility
  meta_data_randomized <- meta_data %>%
    mutate(
      Dev    = factor(sample(as.character(Dev)),    levels = levels(meta_data$Dev)),
      region = factor(sample(as.character(region)), levels = levels(meta_data$region))
    )
  
  sim_df <- merge(all_pcs_df, meta_data_randomized, by = "Sample_ID", all.x = TRUE)
  
  sim_by_window <- split(sim_df, sim_df$window, drop = TRUE)
  
  for (i in names(sim_by_window)) {
    window_sim_df <- sim_by_window[[i]]
    if (nrow(window_sim_df) == 0) next
    
    # Collect window data (take from first row; constant within window)
    chr      <- window_sim_df$chr[1]
    start_bp <- window_sim_df$start_bp[1]
    end_bp   <- window_sim_df$end_bp[1]
    win_id   <- as.integer(i)
    
    # Loop each pc
    for (pc_col in pc_cols) {
      # skip if PCs are NA
      colv <- window_sim_df[[pc_col]]
      if (is.null(colv) || all(is.na(colv))) next
      
      # fit model using type II Anova
      model <- lm(as.formula(paste(pc_col, "~ Dev * region")), data = window_sim_df)
      a_tbl <- as.data.frame(car::Anova(model, type = "II"))
      a_tbl$Term <- rownames(a_tbl)
      
      # keep only terms of interest
      rows <- a_tbl$Term %in% keep_terms
      if (!any(rows)) next
      
      # create tibble with coefficients
      window_pc_coefs <- a_tbl[rows, keep_cols_in, drop = FALSE]
      names(window_pc_coefs) <- keep_cols_out
      
      # add window data
      window_pc_coefs$simulation <- sim
      window_pc_coefs$PC         <- pc_col
      window_pc_coefs$window     <- win_id
      window_pc_coefs$chr        <- chr
      window_pc_coefs$start_bp   <- start_bp
      window_pc_coefs$end_bp     <- end_bp
      
      # add coefficients for each sim+window+pc combo to overall df
      all_sim_coefs[[idx_sim]] <- window_pc_coefs
      idx_sim <- idx_sim + 1
    }
  }
}

# Reformat list into data frame
all_sim_coefs <- dplyr::bind_rows(all_sim_coefs)
sim_df <- all_sim_coefs %>% mutate(source = "Simulated")
write.csv(sim_df, file ="/group/rbaygrp/katie/local_pca/uban_urchins/simulated_df_8_29.csv")

# Combine results and output one final df
all_results <- bind_rows(real_df, sim_df)
write.csv(all_results, file = "/group/rbaygrp/katie/local_pca/uban_urchins/all_results_8_29.csv")

# Load in results and vcf for downstream analysis & plotting  -----------------------------------------------------
real_results <- read.csv("/group/rbaygrp/katie/local_pca/uban_urchins/real_df_8_29.csv")
sim_results <- read.csv("/group/rbaygrp/katie/local_pca/uban_urchins/simulated_df_8_29.csv")

all_results <- bind_rows(real_results, sim_results) %>%
  mutate(
    Term = case_when(
      grepl("^Dev:region", Term) ~ "Interaction",
      grepl("^Dev", Term)        ~ "Urban",
      grepl("^region", Term)     ~ "Region",
      TRUE ~ Term
    ),
    Term = factor(Term, levels = c("Urban", "Region", "Interaction"))
  )
all_results$simulation <- NULL

# Keep only large chroms (filter out scaffolds with less than 400 SNPs )--------
chroms_to_keep <- vcf_fix %>%
  count(CHROM, name = "snp_count") %>%
  filter(snp_count >= 400) %>%
  pull(CHROM)
all_results_chr <- all_results[all_results$chr %in% chroms_to_keep,]
sim_results_chr <- all_results_chr[all_results_chr$source=="Simulated",]
real_results_chr <- all_results_chr[all_results_chr$source=="Real",]

# Look at quantiles of real vs simulated data ----------------------------------
simulated_p <- quantile(sim_results_chr$P, probs=c(0.05,0.01,0.001))
real_p <- quantile(real_results_chr$P, probs=c(0.05,0.01,0.001))
simulated_min_p <- min(sim_results_chr$P, na.rm = T)
real_min_p <- min(real_results_chr$P, na.rm = T)

# Calculate P-values for each term ---------------------------------------------
# For each term and simulation, grab the 99th percentile p-value
quanitle_per_sim <- sim_results_chr %>%
  group_by(source, Term) %>%
  summarize(onepercent = quantile(P, prob = 0.01, na.rm = TRUE))

# Grab the minimum 1% across sims, per term
dev_p <- min(quanitle_per_sim$onepercent[quanitle_per_sim$Term=="Urban"])
reg_p <- min(quanitle_per_sim$onepercent[quanitle_per_sim$Term=="Region"])
int_p <- min(quanitle_per_sim$onepercent[quanitle_per_sim$Term=="Interaction"])

# Make dataframes of only significant P-values
dev_sig_ps <- real_results_chr %>% filter(Term == "Urban", real_results_chr$P<dev_p)
reg_sig_ps <- real_results_chr %>% filter(Term == "Region", real_results_chr$P<reg_p)
int_sig_ps <- real_results_chr %>% filter(Term == "Interaction", real_results_chr$P<int_p)

# Set F-value threshold ---------------------------------------------------
f_quantile_per_sim <- sim_results_chr %>%
  group_by(source, Term) %>%
  summarize(ninetyninth = quantile(`F`, prob = 0.99, na.rm = TRUE), .groups = "drop")

dev_f <- max(f_quantile_per_sim$ninetyninth[f_quantile_per_sim$Term == "Urban"],        na.rm = TRUE)
reg_f <- max(f_quantile_per_sim$ninetyninth[f_quantile_per_sim$Term == "Region"],     na.rm = TRUE)
int_f <- max(f_quantile_per_sim$ninetyninth[f_quantile_per_sim$Term == "Interaction"], na.rm = TRUE)

# Save hits for each term that meet p and f value thresholds
dev_hits <- dev_sig_ps[dev_sig_ps$F>dev_f,]
reg_hits <- reg_sig_ps[reg_sig_ps$F>reg_f,]
int_hits <- int_sig_ps[int_sig_ps$F>int_f,]
all_hits <- bind_rows(dev_hits, reg_hits, int_hits)

dev_non_hits <- real_results_chr[real_results_chr$F<dev_f & real_results_chr$Term=="Urban",]
reg_non_hits <- real_results_chr[real_results_chr$F<reg_f & real_results_chr$Term=="Region",]
int_non_hits <- real_results_chr[real_results_chr$F<int_f & real_results_chr$Term=="Interaction",]
all_non_hits <- bind_rows(dev_non_hits, reg_non_hits, int_non_hits)

dev_mean_f_non_hit <- mean(dev_non_hits$F)
reg_mean_f_non_hit <- mean(reg_non_hits$F)
int_mean_f_non_hit <- mean(int_non_hits$F)

# Reformat and output data to be used to make manhattan plots ------------------
all_results_chr$midpoint = all_results_chr$start_bp+all_results_chr$end_bp/2
real_results_chr <- all_results_chr[all_results_chr$source=="Real",]
urban_results_chr <- real_results_chr[real_results_chr$Term=="Dev",]
region_results_chr <- real_results_chr[real_results_chr$Term=="region",]
int_results_chr <- real_results_chr[real_results_chr$Term=="Dev:region",]

write.csv(urban_results_chr, file = "/group/rbaygrp/katie/local_pca/uban_urchins/urban_pvals.csv")
write.csv(region_results_chr, file = "/group/rbaygrp/katie/local_pca/uban_urchins/region_pvals.csv")
write.csv(int_results_chr, file = "/group/rbaygrp/katie/local_pca/uban_urchins/interaction_pvals.csv")


# Density plots of P-values  ---------------------------------------------------
ggplot(all_results_chr, aes(x = P, color = source, fill = source)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~ Term, scales = "free_y") +
  scale_x_continuous(limits = c(0, 0.05)) +
  labs(title = "ANOVA p-value distributions (Real vs Simulated)",
       x = "p-value", y = "Density") +
  scale_color_manual(values = c("Real" = "blue", "Simulated" = "grey40")) +
  scale_fill_manual(values  = c("Real" = "blue", "Simulated" = "grey70")) +
  theme_minimal()

# Density plot of F-statistics ------------------------------
term_cols <- c("Urban" = "#665D4B", "Region" = "deeppink2", "Interaction" = "limegreen")

ggplot(all_hits, aes(x = `F`, color = Term, fill = Term)) +
  geom_density(alpha = 0.2, adjust = 1) +
  labs(
    title = "F-statistic Distribution of Significant Windows",
    x = "F statistic", y = "Density",
    color = "Term", fill = "Term"   # legend titles
  ) +
  scale_color_manual(values = term_cols) +
  scale_fill_manual(values  = term_cols) +
  theme_box()


ggplot(all_results_chr, aes(x = `F`, color = Term, fill = Term)) +
  geom_density(alpha = 0.3) +
  scale_color_manual(values = term_cols, name = "Term") +
  scale_fill_manual(values  = term_cols, name = "Term") +
  labs(
    title = "F-value distribution of all Windows",
    x = "F statistic", y = "Density"
  ) +
  theme_box()

ggplot(all_non_hits, aes(x = `F`, color = Term, fill = Term)) +
  geom_density(alpha = 0.3) +
  scale_color_manual(values = term_cols, name = "Term") +
  scale_fill_manual(values  = term_cols, name = "Term") +
  labs(
    title = "F-value distribution of non-significant windows",
    x = "F statistic", y = "Density"
  ) +
  theme_box()
