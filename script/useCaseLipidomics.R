library(mrmhub)
library(rgoslin)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(ggnewscale)
library(ggpmisc)
library(ggbeeswarm)
library(mirai)
library(ggvenn)
library(dplyr)
library(rmzqc)
library(tidyr)
library(patchwork)
library(here)

here("script")

# MzQC - Lipidomics Use Case
# Dataset of MRMhub provided by SLING available at zenodo: https://zenodo.org/records/15370294 (Dataset: MRMhub-Dataset3.zip; Input files for this workflow: mrmhub-workflows)
# This Use Case goes through the Workflow from MRMhub for Dataset3: https://slinghub.github.io/mrmhub-workflows/Dataset3.html#peak-picking-qc
# Additionally exports QCs to mzQC files to prove the usability of the file format for a Lipidomics Use Case

data_path <- "../data/Dataset3_MRMhub-Integrator_20251010.csv"
mexp <- MRMhubExperiment()
mexp <- import_data_mrmhub(mexp, data_path, import_metadata = TRUE)       
file_path <- "../data/Dataset3_Metadata_20251010.xlsm"
mexp <- import_metadata_msorganiser(mexp, file_path, ignore_warnings = TRUE)

# Analytical design and timeline
plot_runsequence(
  mexp, 
  show_batches = TRUE, 
  qc_types = c("SPL", "BQC", "TQC", "PBLK", "UBLK", "RQC", "SBLK", "LTR", "NIST"),
  batch_zebra_stripe = TRUE, base_font_size = 6,
  batch_fill_color = "#fffbdb", 
  segment_linewidth = 0.25,
  show_timestamp = FALSE) +
  theme(plot.title = element_text(size = 5))

# Overview Chromatographic Separation
plot_abundanceprofile(
  data = mexp,
  log_scale = FALSE,
  variable = "rt", 
  density_strip = TRUE,
  qc_types = "SPL",
  analysis_range = c(1,4000),
  show_sum = FALSE,
  #x_lim = c(6.5, 7.5),
  x_label = NA,
  feature_map = "lipidomics")

# Overall trends and check for possible outliers
fig2a <- plot_rla_boxplot(
  data = mexp,
  rla_type_batch = c("within"),
  variable = "intensity",
  qc_types = c("BQC", "TQC", "SPL", "NIST", "LTR", "SBLK"),  
  #plot_range = c(1, 4900),
  rla_limit_to_range = FALSE,
  filter_data = FALSE,
  min_feature_intensity = 1000,
  include_feature_filter = "ISTD",
  #y_lim = c(-4,4),
  show_timestamp = FALSE,
  outlier_method = "fold",
  outlier_k = c(-0.4, 0.3),
  outlier_exclude = FALSE, 
  x_gridlines = FALSE,
  show_plot = TRUE,
  batch_zebra_stripe = TRUE,
  linewidth = 0.1
)
# Get outlier in ISTD total signal excluding 3 batches with overall lower ISTD
istd_outlier <- fig2a$outliers |> 
  filter(!(batch_id %in% c("P-01", "P-02", "P-43") & val_res_median > -3))

print(istd_outlier)

# Plot peak areas of all ISTDs in different QC types against the run order
plot_runscatter(mexp, variable = "intensity", 
                include_qualifier = FALSE,
                qc_types = c("SPL", "BQC", "TQC", "PBLK", "RQC", "SBLK"),
                include_feature_filter = "IS",
                #y_min = 0.00, y_max = 0.15,
                #plot_range = c(0, 910),
                point_size = .2,
                point_border_width = 0.1,
                point_transparency = .7,
                base_font_size = 4,
                cols_page = 3,
                rows_page = 11,
                show_progress = FALSE,
                cap_outliers = TRUE)

# Check for potential systematic outliers in endogenous analytes
fig2b <- plot_rla_boxplot(
  data = mexp,
  rla_type_batch = c("within"),
  variable = "intensity",
  qc_types = c("BQC", "TQC", "SPL", "NIST", "LTR"),  
  filter_data = FALSE,
  min_feature_intensity = 1000,
  exclude_feature_filter = "ISTD",
  plot_range = c(1, 4800),
  #y_lim = c(-4,4),
  show_timestamp = FALSE,
  outlier_method = "fold",
  outlier_k = c(-2,1) ,
  outlier_exclude = FALSE, 
  x_gridlines = FALSE,
  batch_zebra_stripe = TRUE,
  show_plot = TRUE,
  linewidth = 0.1
) 
# Get outlier in ISTD total signal excluding 3 batches with conistently lower ISTD
analyte_outlier <- fig2b$outliers |> 
  filter(!(batch_id %in% c("P-01", "P-02", "P-43") & val_res_median > -3))
print(analyte_outlier)

# Outliers in analyte but not ISTD signals
setdiff(
  analyte_outlier$analysis_id,
  istd_outlier$analysis_id 
)


# Combined ISTD and analyte outlier
outlier_combined <- istd_outlier |> 
  bind_rows(analyte_outlier) |>
  bind_rows(tibble(analysis_id = "NIST-04", qc_type = "NIST")) |> 
  select(-val_res_median) |>
  distinct()

outlier_combined |> dplyr::count(qc_type)

mexp <- exclude_analyses(mexp, analyses = outlier_combined$analysis_id, clear_existing = TRUE)

# Plot the ISTDs again to check if the outliers were removed
plot_runscatter(mexp, variable = "intensity", 
                include_qualifier = FALSE,
                qc_types = c("SPL", "BQC", "TQC", "PBLK", "RQC", "SBLK"),
                include_feature_filter = "IS",
                #y_min = 0.00, y_max = 0.15,
                #plot_range = c(0, 910),
                point_size = .2,
                point_border_width = 0.1,
                point_transparency = .7,
                base_font_size = 4,
                cols_page = 3,
                rows_page = 11,
                show_progress = FALSE,
                cap_outliers = TRUE)



# Peak picking QC
plot_rt_vs_chain(
  mexp, 
  qc_types = "SPL", 
  x_var = "total_c", outlier_residual_min = 0.3,
  base_font_size <- 6, 
  cols_page = 4,
  point_size = 1) +
  theme( 
    legend.position = "right",
    legend.direction = "vertical", 
    legend.text = element_text(size = base_font_size*0.8),
    legend.title = element_text(size = base_font_size*0.8),
    legend.key.size = unit(base_font_size *0.8, "pt"),
    legend.position.inside = c(0.1, 0.7)
  )  


# Feature correlation analysis
# this below is to exclude a sample that has a very low intensity for all features, see next steps for details
plot_feature_correlations(
  mexp, 
  variable = "intensity" , 
  qc_types = c("SPL", "BQC", "TQC"), 
  point_size = 0.5, 
  cor_min = 0.98,
  point_stroke = 0.1,
  sort_by_corr = TRUE, 
  return_plot = TRUE,
  show_progress = FALSE,
  log_scale = TRUE, 
  cols_page = 5, 
  rows_page = 6, 
  font_base_size = 5)

# Exclude PC 0-36:1 due to the truncation of its peak
mexp <- exclude_features(mexp, features = "PC(O-36:1)", clear_existing = TRUE)



# Summing up LysoPL and DG isomers
mexp <- data_sum_features(mexp, feature_classes = c("LPC", "LPE", "LPG", "LPI", "LPS", "DG"))


# PCA to check for potential technical outliers
plot_pca(
  data = mexp,
  variable = "intensity",
  filter_data = FALSE,
  pca_dims = c(1,2),
  labels_threshold_mad = 5,
  labels_column = "analysis_order",
  qc_types = c("BQC", "TQC", "LTR", "NIST", "SPL"),
  ellipse_variable = "qc_type",
  log_transform = TRUE,shared_labeltext_hide = "_MS-5",
  point_size = 0.7, point_alpha = 0.7, font_base_size = 8, ellipse_alpha = 0.3,
  include_istd = FALSE,show_labels = TRUE,label_font_size = 1.5)   + theme(
    plot.title = element_blank(),
    aspect.ratio = 1,
    legend.position = "inside",
    legend.direction = "vertical", 
    legend.text = element_text(size = 9*0.7),
    legend.title = element_text(size = 9*0.7),
    legend.key.size = unit(7 *0.7, "pt"),
    legend.position.inside = c(0.83, 0.83)) 


# Loading
plot_pca_loading(
  data = mexp,
  variable = "feature_intensity",
  include_istd = FALSE,
  pca_dims = c(1,2,3,4),
  top_n = 70,
  font_base_size = 7,
  #qc_types = c("SPL", "BQC", "TQC", "LTR"),
  log_transform = TRUE,
  #point_size = 1, point_alpha = 0.7, font_base_size = 8, ellipse_alpha = 0.3,
  #include_istd = FALSE,
  #show_labels = TRUE,label_font_size = 2,
  #shared_labeltext_hide = NA
)   + theme(
  plot.title = element_blank(),
  legend.position = "inside",
  legend.direction = "vertical", 
  legend.text = element_text(size = 8*0.7),
  legend.title = element_text(size = 8*0.7),
  legend.key.size = unit(6 *0.7, "pt"),
  legend.position.inside = c(0.96, 0.06)) 



# Matrix efffects
plot_qc_matrixeffects(
  mexp, 
  variable = "intensity",
  batchwise_normalization = TRUE,
  only_istd = FALSE,include_qualifier = FALSE,
  include_feature_filter = "ISTD",
  exclude_feature_filter = "95|CL|25|d17\\:0|C1P|H2O",
  y_lim = c(50, 150),
  point_alpha = 0.05, 
  box_alpha = 0.3, 
  point_size = 0.2, 
  box_linewidth = 0.2,
  font_base_size = 7,
  min_median_value = 1000) 



# Isotope correction
mexp <- correct_interferences(mexp)
plot_qc_interferences(mexp, y_lim = c(-10, 110), qc_types = c("LTR", "NIST", "SPL", "TQC", "BQC"))



# Normalization
mexp <- normalize_by_istd(mexp, ignore_missing_annotation = FALSE) 


# Quantification
mexp <- quantify_by_istd(mexp) 


# PCA plot for normalized and quantified data
plot_pca(
  data = mexp,
  variable = "conc",
  filter_data = FALSE,
  pca_dims = c(1,2),labels_column = "analysis_order",
  labels_threshold_mad = 4,
  qc_types = c("BQC", "TQC", "LTR", "NIST", "SPL"),
  ellipse_variable = "qc_type",
  log_transform = TRUE,
  point_size = 0.7, point_alpha = 0.7, font_base_size = 8, ellipse_alpha = 0.3,
  include_istd = FALSE,show_labels = TRUE,label_font_size = 1.5)


# Remove batches with identified technical issues
ids_reruns <- mexp@annot_analyses |> filter(batch_id %in% c("P-43", "P-42")) |> pull(analysis_id)

# ToDo: Currently exclude_analyses usage requires re-processing of the data. Needs to be updated.
# mexp_final <- exclude_analyses(mexp, analyses = ids_reruns, clear_existing = FALSE)

# Workaround is to filter the dataset directly
mexp@dataset <- mexp@dataset |> filter(!analysis_id %in% ids_reruns)
mexp@annot_analyses <- mexp@annot_analyses |> filter(!analysis_id %in% ids_reruns)
mexp@annot_responsecurves <- mexp@annot_responsecurves |> filter(!analysis_id %in% ids_reruns)


# Inspection of individual features concentrations across the run order
plot_runscatter(mexp, 
                variable = "conc", filter_data = FALSE,
                qc_types = c("SPL", "BQC", "TQC", "LTR"), 
                include_feature_filter = "PC 34\\:2|PC 36\\:4|PE 38\\:4|PI 38\\:4|
    PC 38\\:6|CE 181\\:1|Cer d18\\:1\\/16\\:0|TG 52\\:3|SM 34\\:1|
    TG 56\\:3|LPC 18\\:1",
                #y_min = 0.00,
                #y_max = 0.15,
                #plot_range = c(0, 910), 
                show_reference_lines = TRUE,ref_qc_types = "SPL",
                reference_fill_color = "#111111",
                reference_k_sd = 3,
                point_size = 0.5,
                point_border_width = 0.1,
                base_font_size = 6,
                cols_page = 2,
                rows_page = 6,
                cap_outliers = TRUE,
                reference_sd_shade = FALSE,
                #batch_zebra_stripe = TRUE, 
                output_pdf = FALSE,
                path = "rt_all.pdf")


# Drift and Batch correction
mexp <- correct_drift_gaussiankernel(
  mexp, 
  variable = "conc", 
  ref_qc_types = "SPL", 
  batch_wise = TRUE, 
  kernel_size = 10, 
  outlier_filter = TRUE, 
  outlier_ksd = 5, 
  recalc_trend_after = TRUE, 
  show_progress = FALSE) 

mexp <- correct_batch_centering(
  mexp, 
  ref_qc_types = "SPL", 
  variable = "conc")


# Runscatter Plot
plot_runscatter(mexp, variable = "conc", 
                filter_data = FALSE,
                #include_feature_filter = "ISTD", 
                qc_types = c("SPL", "BQC", "TQC", "LTR"), 
                include_feature_filter = "PC 34\\:2|PC 36\\:4|PE 38\\:4|PI 38\\:4|
    PC 38\\:6|CE 181\\:1|Cer d18\\:1\\/16\\:0|TG 52\\:3|SM 34\\:1|
    TG 56\\:3|LPC 18\\:1",
                #y_min = 0.00,
                #y_max = 0.15,
                #plot_range = c(0, 910), 
                show_reference_lines = TRUE,ref_qc_types = "SPL",
                reference_fill_color = "#111111",
                reference_k_sd = 3,
                point_size = 0.5,
                point_border_width = 0.1,
                base_font_size = 6,
                cols_page = 2,
                rows_page = 6,
                cap_outliers = TRUE,
                reference_sd_shade = FALSE,
                #batch_zebra_stripe = TRUE, 
                output_pdf = FALSE,
                path = "rt_all.pdf")




# Runscatter plots of final concentrations of all features
plot_runscatter(mexp, variable = "conc_raw", 
                #include_feature_filter = "ISTD", 
                qc_types = c("SPL", "BQC", "TQC", "LTR"), 
                #include_feature_filter = example_species,
                #y_min = 0.00,
                #y_max = 0.9,
                #plot_range = c(0, 910), 
                show_reference_lines = TRUE,
                ref_qc_types = "SPL",
                reference_fill_color = "#111111",
                reference_k_sd = 3,
                show_trend = TRUE,
                point_size = 1,
                base_font_size = 6,
                cols_page = 2,
                rows_page = 3,
                cap_outliers = FALSE,
                reference_sd_shade = FALSE,
                show_progress = FALSE,
                #batch_zebra_stripe = TRUE, 
                output_pdf = TRUE,
                path = "../output/Dataset3_runscatter_rawConc_all.pdf")

plot_runscatter(mexp, variable = "conc", 
                #include_feature_filter = "ISTD", 
                qc_types = c("SPL", "BQC", "TQC", "LTR"), 
                #include_feature_filter = example_species,
                #y_min = 0.00,
                #y_max = 0.9,
                #plot_range = c(0, 910), 
                show_reference_lines = TRUE,
                ref_qc_types = "SPL",
                reference_fill_color = "#111111",
                reference_k_sd = 3,
                show_trend = TRUE,
                point_size = 1,
                base_font_size = 6,
                cols_page = 2,
                rows_page = 3,
                cap_outliers = FALSE,
                show_progress = FALSE,
                reference_sd_shade = FALSE,
                #batch_zebra_stripe = TRUE, 
                output_pdf = TRUE,
                path = "../output/Dataset3_runscatter_FinalConc_all.pdf")


# QC of normalization and drift/batch correction
mexp <- calc_qc_metrics(mexp, use_robust_cv = FALSE, use_batch_medians = TRUE)

plot_normalization_qc(
  plot_type = "diff",
  data = mexp,
  before_norm_var = "intensity",
  after_norm_var = "conc",
  y_lim = c(-15,15), 
  x_lim = c(0,75),
  qc_types = c("TQC", "BQC", "SPL", "NIST"),
  cols_page = 5,
  font_base_size = 5,
  point_size = 0.5,
  facet_by_class = TRUE,
  include_qualifier = FALSE,
) 


# Process vs instrument variability
mrmhub::plot_qcmetrics_comparison(
  mexp,
  plot_type = "diff",
  y_shared = TRUE,
  "conc_cv_tqc",
  "conc_cv_bqc", 
  log_scale = FALSE,
  equality_line = TRUE,
  facet_by_class = TRUE,
  point_size = 2,
  font_base_size = 5,
  x_lim = c(0, 25),
  y_lim = c(-15, 15))


# Effect of feature intensity on technical variability
mrmhub::plot_qcmetrics_comparison(
  mexp,
  plot_type = "scatter",
  y_shared = FALS,
  "intensity_median_bqc",
  "intensity_cv_bqc", 
  log_scale = FALSE,
  equality_line = FALSE,
  facet_by_class = FALSE,
  point_size = 0.5,
  font_base_size = 5,
  x_lim = c(10, Inf),
  y_lim = c(0, 100))+ 
  ggplot2::geom_smooth(method = "loess", se = FALSE, span = 0.75)+
  geom_hline(yintercept = 20, linetype = "dashed", color = "grey70") +
  scale_x_log10(expand = ggplot2::expansion(mult = c(0, 0.00)),
                breaks = c(1E2, 1E3, 1E4, 1E5, 1E6, 1E7, 1E8)) +
  scale_y_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.00)),breaks = c(0, 20, 40, 60, 80))





# Response curves
sel_species <- c("PC 26:0 (ISTD)", "SM 34:1", "PC 34:2", "CE 18:2", "TG 50:1 [NL-16:0]", "CE 20:4", "LPC 18:1")
plot_responsecurves(
  data = mexp,
  variable = "intensity", 
  max_regression_value = 100, 
  filter_data = FALSE,
  font_base_size = 6, 
  line_width = 0.5, point_size = 1.2,
  include_feature_filter = sel_species,
  output_pdf = FALSE, path = "response-curves-dataset3.pdf",
  show_progress = FALSE,
  cols_page = 3, 
  rows_page = 3,
  return_plots = TRUE)[[1]]+
  theme(panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.direction = "horizontal", 
        legend.text = element_text(size = 6*0.6),
        legend.title = element_blank(),
        legend.key.size = unit(6 *0.6, "pt"),
        strip.text = element_text(size = 6),
        legend.position.inside = c(0.9, 0.01))




# Generate pdf files of response curves
plot_responsecurves(
  data = mexp,
  variable = "intensity", max_regression_value = 100, 
  filter_data = FALSE,font_base_size = 6, line_width = 0.5, point_size = 1.2,
  output_pdf = TRUE, path = "../output/dataset3-response-curves.pdf",show_progress = FALSE,
  cols_page = 6, rows_page = 5,
  return_plots = FALSE)




# Feature filter
mexp <- calc_qc_metrics(mexp, use_robust_cv = FALSE, use_batch_medians = TRUE, include_response_stats = TRUE)

mexp  <- filter_features_qc(
  data =  mexp, 
  clear_existing = TRUE,
  use_batch_medians = TRUE,
  include_qualifier = FALSE,
  include_istd = FALSE,
  response.curves.selection = 1,
  response.curves.summary = "mean",
  min.rsquare.response = 0.8,
  min.slope.response = 0.5,
  max.yintercept.response = 0.5,
  min.signalblank.median.spl.pblk = 10,
  min.intensity.median.spl = 100,
  max.cv.conc.bqc = 25,
  #max.dratio.sd.conc.bqc = 0.75,
  max.prop.missing.conc.spl = 100,
  features.to.keep = c("CE 20:4", "CE 22:5", "CE 22:6", "CE 16:0", "CE 18:0")
)
table(mexp@metrics_qc$feature_class)

# Results of feature filter
 plot_qc_summary_byclass(mexp) + 
  theme(
    strip.text = ggplot2::element_text(size = 5),
    #aspect.ratio = 0.9,
    legend.position = "inside",
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8, face = "plain"),
    axis.text.y.right =  element_text(size = 7, face = "plain"),
    legend.direction = "vertical", 
    legend.text = element_text(size = 8*0.7),
    legend.title = element_blank(),
    legend.key.size = unit(8 *0.7, "pt"),
    legend.position.inside = c(0.77, 0.27)
  )  



plot_qc_summary_overall(mexp) 


# Lipidome profile
plot_abundanceprofile(
  data = mexp,
  log_scale = TRUE,
  filter_data = TRUE,
  variable = "conc", 
  qc_types = "SPL",
  #x_lim = c(-6, 2),
  x_label = NA,
  feature_map = "lipidomics")


# Export final dataset
save_dataset_csv(
  data = mexp, 
  path = "../output/sperfect_UNFILTERED_RAW-feature_conc_uM_20250518a.csv",
  variable = "conc", 
  qc_types = "SPL", 
  include_qualifier = FALSE,
  filter_data = FALSE)


# End of MRMhub workflow


# ---------------------------------------------------------------------------------------------------------------------------------------------

# MzQC Export


# Compute example QC metrics
# TO-DO: Discuss and update the metrics selection and its terms

# Prepare input
feature_qc_metrics <- mexp@metrics_qc %>%
  select(
    feature_id,
    
    # pooled QC CVs
    conc_cv_bqc,
    conc_cv_tqc,
    
    # Signal to blank
    sb_ratio_pblk,
    
    # Linearity
    r2_rqc_B,
    slopenorm_rqc_B,
    y0norm_rqc_B
  )

# Define input file metadata
# TO-DO: Discuss which term to include in the psi-ms-CV https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo
inputFiles = list(
  rmzqc::MzQCinputFile$new(
    name = "Dataset3_MRMhub-Integrator_20251010.csv",
    location = paste0("file://", data_path),
    fileFormat = rmzqc::MzQCcvParameter$new(
      # TO-DO: discuss accession number
      accession = "MS:1001040",
      name = "intermediate analysis format"
    )
  )
)

metadata <- rmzqc::MzQCmetadata$new(
  label = "MRMhub Dataset3",
  
  inputFiles = inputFiles,
  
  analysisSoftware = list(
    rmzqc::MzQCanalysisSoftware$new(
      name = "MRMhub",
      version = as.character(
        packageVersion("mrmhub")
      ),
      accession = "MS:XXXXXXX", # to be defined
      uri = "https://github.com/SLINGhub/MRMhub/"
    )
  )
)


quality_metrics <- list(
  # BQC CV
  rmzqc::MzQCqualityMetric$new(
    accession = "MS:XXXXXXX", # to be defined
    name = "pooled batch quality control concentration coefficient of variation",
    value = list(
      feature_id = feature_qc_metrics$feature_id,
      value = feature_qc_metrics$conc_cv_bqc
    )
  ),
  # TQC CV
  rmzqc::MzQCqualityMetric$new(
    accession = "MS:XXXXXXX", # to be defined
    name = "pooled technical quality control concentration coefficient of variation",
    value = list(
      feature_id = feature_qc_metrics$feature_id,
      value = feature_qc_metrics$conc_cv_tqc
    )
  ),
  
  # Signal-to-blank
  # TO-DO: Add term to in the psi-ms-CV
  rmzqc::MzQCqualityMetric$new(
    accession = "MS:XXXXXXX", # to be defined
    name = "median-based signal-to-blank ratio",
    value = list(
      feature_id = feature_qc_metrics$feature_id,
      value = feature_qc_metrics$sb_ratio_pblk
    )
  ),
  
  # Linearity
  # TO-DO: Add term to in the psi-ms-CV
  rmzqc::MzQCqualityMetric$new(
    accession = "MS:XXXXXXX", # to be defined
    name = "response curve linearity coefficient of determination (R²)",
    value = list(
      feature_id = feature_qc_metrics$feature_id,
      value = feature_qc_metrics$r2_rqc_B
    )
  ),
  
  rmzqc::MzQCqualityMetric$new(
    accession = "MS:XXXXXXX", # to be defined
    name = "response curve linearity normalized slope",
    value = list(
      feature_id = feature_qc_metrics$feature_id,
      value = feature_qc_metrics$slopenorm_rqc_B
    )
  ),
  
  rmzqc::MzQCqualityMetric$new(
    accession = "MS:XXXXXXX", # to be defined
    name = "response curve linearity normalized intercept",
    value = list(
      feature_id = feature_qc_metrics$feature_id,
      value = feature_qc_metrics$y0norm_rqc_B
    )
  )
)

set_quality <- rmzqc::MzQCsetQuality$new(
  metadata = metadata,
  qualityMetrics = quality_metrics
)

controlled_vocabularies <- list(
  rmzqc::MzQCcontrolledVocabulary$new(
    name = "Proteomics Standards Initiative Mass Spectrometry Controlled Vocabulary",
    uri = "https://github.com/HUPO-PSI/psi-ms-CV",
    version = "4.1.257"
  )
)

mzqc_obj <- rmzqc::MzQCmzQC$new(
  version = "1.0.0",
  creationDate = rmzqc::MzQCDateTime$new(),
  contactName = "Franziska Nicolaus",
  contactAddress = "f.nicolaus@fz-juelich.de",
  description = "MRMhub lipidomics workflow mzQC export containing pooled QC concentration coefficients of variation, signal-to-blank ratios, and response curve linearity metrics.",
  runQualities = list(),
  setQualities = list(set_quality),
  controlledVocabularies = controlled_vocabularies
)

output_file <- "../output/Dataset3_testQCmetrics.mzQC"

rmzqc::writeMZQC(
  filepath = output_file,
  mzqc = mzqc_obj
)
cat("Yay! The mzQC file has been successfully written to:", output_file, "\n")



# Further Analysis for the use case


# Histograms

# BQC CV
p1 <- ggplot(feature_qc_metrics, aes(x = conc_cv_bqc)) +
  geom_histogram(
    bins = 30,
    fill = "#4C78A8",
    color = "white"
  ) +
  scale_x_continuous(
    breaks = seq(0, 170, by = 20),
    limits = c(0, 170)
  ) +
  labs(
    title = "BQC Coefficient of Variation",
    x = "Concentration CV (%)",
    y = "Number of features"
  ) +
  theme_qc

# TQC CV
p2 <- ggplot(feature_qc_metrics, aes(x = conc_cv_tqc)) +
  geom_histogram(
    bins = 30,
    fill = "#59A14F",
    color = "white"
  ) +
  scale_x_continuous(
    breaks = seq(0, 170, by = 20),
    limits = c(0, 170)
  ) +
  labs(
    title = "TQC Coefficient of Variation",
    x = "Concentration CV (%)",
    y = "Number of features"
  ) +
  theme_qc

# Signal-to-Blank
p3 <- ggplot(feature_qc_metrics, aes(x = sb_ratio_pblk)) +
  geom_histogram(
    bins = 30,
    fill = "#F28E2B",
    color = "white"
  ) +
  scale_x_log10() +
  labs(
    title = "Signal-to-Blank",
    x = "Median Signal-to-Blank Ratio (log10)",
    y = "Number of features"
  ) +
  theme_qc

# R²
p4 <- ggplot(feature_qc_metrics, aes(x = r2_rqc_B)) +
  geom_histogram(
    bins = 20,
    fill = "#E15759",
    color = "white"
  ) +
  labs(
    title = expression(R^2),
    x = expression(R^2),
    y = "Number of features"
  ) +
  theme_qc

# Arrange in a grid
(p1 | p2) / (p3 | p4)



# Additional calculations

# Mean / median of CVs
qc %>%
  summarise(
    BQC_median = median(conc_cv_bqc, na.rm = TRUE),
    BQC_mean = mean(conc_cv_bqc, na.rm = TRUE),
    TQC_median = median(conc_cv_tqc, na.rm = TRUE),
    TQC_mean = mean(conc_cv_tqc, na.rm = TRUE)
  )

# Bar CV + N of Lipids
hist_data_bqc <- hist(feature_qc_metrics$conc_cv_bqc, plot = FALSE, breaks = 30)
i <- which.max(hist_data_bqc$counts)

cat(
  "Highest bar BQC CV:",
  hist_data_bqc$breaks[i], "-",
  hist_data_bqc$breaks[i + 1], "% CV\n",
  "Lipids:", hist_data_bqc$counts[i], "\n"
)

hist_data_tqc <- hist(feature_qc_metrics$conc_cv_tqc, plot = FALSE, breaks = 30)
i <- which.max(hist_data_tqc$counts)

cat(
  "Highest bar TQC CV:",
  hist_data_tqc$breaks[i], "-",
  hist_data_tqc$breaks[i + 1], "% CV\n",
  "Lipids:", hist_data_tqc$counts[i], "\n"
)

