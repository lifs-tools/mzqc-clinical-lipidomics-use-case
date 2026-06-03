# target lipidomics workflow: Lipidomics Data Processing
# MRMhub developed by SLING

# Link to tutorial: https://slinghub.github.io/MRMhub/articles/T01_targetlipidomics_workflow.html

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


myexp <- mrmhub::MRMhubExperiment(title = "sPerfect")

# Data obtained from MRM_version_20260315 https://github.com/SLINGhub/MRMhub/releases
data_path <- "../data/long.csv"
myexp <- import_data_mrmhub(data = myexp, path = data_path, import_metadata = TRUE)
print(myexp)

# A glimpse into imported data
# View(myexp@dataset_orig)
head(myexp@dataset_orig)

#  Analytical design and timeline
plot_runsequence(
  myexp, 
  qc_types = NA, 
  show_batches = TRUE, 
  batch_zebra_stripe = TRUE, 
  batch_fill_color = "#fffbdb", 
  segment_linewidth = 0.5,
  show_timestamp = FALSE)

# Preview chromatic separation
plot_abundanceprofile(
  data = myexp,
  log_scale = FALSE,
  variable = "rt", 
  density_strip = TRUE,
  qc_types = "SPL",
  show_sum = FALSE,
  x_label = NA,
  feature_map = "lipidomics")

# Peak picking QC
plot_rt_vs_chain(
  myexp, 
  qc_types = "SPL", 
  x_var = "total_c", outlier_residual_min = 0.3,
  base_font_size = 8, 
  point_size = 1)  

# Signal trends of Internal Standards
plot_runscatter(
  data = myexp,
  variable = "intensity",
  qc_types = c("BQC", "TQC", "SPL", "PBLK", "SBLK"),
  #analysis_range = NA, #get_batch_boundaries(myexp, c(1,6)), 
  include_feature_filter = "ISTD", 
  exclude_feature_filter = "Hex|282",
  cap_outliers = TRUE,
  log_scale = FALSE, 
  show_batches = TRUE,base_font_size = 5,
  output_pdf = TRUE, # change here if saving as PDF favoured or not
  path = "../output/runscatter_istd.pdf",
  cols_page = 4, rows_page = 3
)


# Adding detailed metadata
file_path <- "../data/metadata-template.xlsx"
myexp <- import_metadata_msorganiser(myexp, path = file_path, ignore_warnings = TRUE)
myexp <- set_analysis_order(myexp, order_by = "timestamp")
myexp <- set_intensity_var(myexp, variable_name = "area")


# Overall trends and possible outlier
mrmhub::plot_rla_boxplot(
  data = myexp,
  rla_type_batch = c("within"),
  variable = "intensity",
  qc_types = c("BQC", "SPL", "RQC", "TQC", "PBLK"), 
  filter_data = FALSE, 
  #analysis_range = get_batch_boundaries(myexp, batch_indices = c(5,6)), 
  #y_lim = c(-3,3),
  show_timestamp = FALSE,
  outlier_exclude= FALSE, x_gridlines = FALSE,
  batch_zebra_stripe = FALSE,
  linewidth = 0.1
)


# PCA plot of all QC types
plot_pca(
  data = myexp, 
  variable = "feature_intensity", 
  filter_data = FALSE,
  pca_dims = c(1,2),
  labels_threshold_mad = 3, 
  qc_types = c("SPL", "BQC", "TQC"), #, "PBLK", "RQC"
  log_transform = TRUE,  
  point_size = 2, point_alpha = 0.7, font_base_size = 8, ellipse_alpha = 0.3, 
  include_istd = FALSE)


# Exclude technical outliers
# Exclude the sample from the processing
myexp <- exclude_analyses(myexp, analyses = c("Longit_batch6_51"), clear_existing  = TRUE)

# Replot the PCA
plot_pca(
  data = myexp,
  variable = "intensity",
  filter_data = FALSE,
  pca_dim = c(1,3),
  labels_threshold_mad = 3,
  qc_types = c("SPL", "BQC", "TQC"),
  log_transform = TRUE,
  point_size = 2, point_alpha = 0.7, font_base_size = 8, ellipse_alpha = 0.3,
  include_istd = FALSE,
  shared_labeltext_hide = NA)


# Response curves
# Exclude very low abundant features
myexp <- mrmhub::filter_features_qc(myexp, 
                                    include_qualifier = FALSE,
                                    include_istd = TRUE,
                                    min.intensity.median.spl = 200)
plot_responsecurves(
  data = myexp,
  variable = "intensity",
  filter_data = TRUE,
  #include_feature_filter = "^PC 3[0-5]", # here we use regular expressions
  output_pdf = TRUE, path = "../output/response-curves.pdf",
  cols_page = 5, rows_page = 4,
)


# Isotope interference correction
myexp <- mrmhub::correct_interferences(myexp)
plot_qc_interferences(myexp, qc_types = c("BQC", "SPL", "TQC", "LTR"))


# Normalization and quantification based on ISTDs
myexp <- mrmhub::normalize_by_istd(myexp)
myexp <- mrmhub::quantify_by_istd(myexp)


# Examine the effects of class-wide ISTD normalization
myexp <- mrmhub::filter_features_qc(myexp, 
                                    include_qualifier = FALSE,
                                    include_istd = TRUE,
                                    min.intensity.median.spl = 1000)
mrmhub::plot_normalization_qc(
  data = myexp, 
  before_norm_var = "intensity", 
  after_norm_var = "norm_intensity",
  plot_type = "diff",
  point_size = 2,
  facet_by_class = TRUE,
  qc_types = c("TQC", "BQC", "SPL"),
  y_lim = c(-5, 15))


# Drift correction
myexp <- mrmhub::correct_drift_gaussiankernel(
  data = myexp,
  ignore_istd = TRUE,
  variable = "conc",
  ref_qc_types = c("SPL"),
  batch_wise = TRUE,
  replace_previous = TRUE,
  recalc_trend_after = TRUE,
  kernel_size = 10,
  #conditional_correction = 
  outlier_filter = FALSE,
  outlier_ksd = 5,
  location_smooth = TRUE,
  scale_smooth = FALSE, 
  show_progress = FALSE  # set to FALSE when rendering
)


# Define a wrapper function
my_trend_plot <- function(variable, feature){
  plot_runscatter(
    data = myexp,
    variable = variable,
    qc_types = c("BQC", "TQC", "SPL"),
    include_feature_filter = feature,
    exclude_feature_filter = "ISTD",
    cap_outliers = TRUE,
    log_scale = FALSE,
    show_trend = TRUE,
    output_pdf = FALSE,
    path = "../output/runscatter_PC408_beforecorr.pdf",
    cols_page = 1, rows_page = 1, 
  )
}

my_trend_plot("conc_before", "PC 40:8")
my_trend_plot("conc", "PC 40:8")


# Batch effect correction
myexp <- mrmhub::correct_batch_centering(
  myexp, 
  variable = "conc",
  ref_qc_types = "SPL",
  replace_previous = TRUE,
  correct_location = TRUE, 
  correct_scale = TRUE, 
  log_transform_internal = TRUE)

my_trend_plot("conc", "PC 40:8")


# Saving runscatter plots of all features as PDF
plot_runscatter(
  data = myexp,
  variable = "conc",
  qc_types = c("BQC", "TQC", "SPL"),
  include_feature_filter =  NA,
  exclude_feature_filter = "ISTD",
  cap_outliers = TRUE,
  log_scale = FALSE,
  show_trend = TRUE,
  output_pdf = TRUE,
  path = "../output/runscatter_after-drift-batch-correction.pdf",
  cols_page = 2, 
  rows_page = 2,
  show_progress = TRUE
)


# QC-based feature filtering
myexp <- filter_features_qc(
  data = myexp, 
  clear_existing = TRUE,
  use_batch_medians = TRUE,
  include_qualifier = FALSE,
  include_istd = FALSE,
  response.curves.selection = c(1,2),
  response.curves.summary = "mean",
  min.rsquare.response = 0.8,
  min.slope.response = 0.75,
  max.yintercept.response = 0.5,
  min.signalblank.median.spl.pblk = 10,
  min.intensity.median.spl = 100,
  max.cv.conc.bqc = 25,
  features.to.keep = c("CE 20:4", "CE 22:5", "CE 22:6", "CE 16:0", "CE 18:0")
)


# Summary of the filtering
mrmhub::plot_qc_summary_byclass(myexp) 
mrmhub::plot_qc_summary_overall(myexp)


# Lipidome Profile
plot_abundanceprofile(
  data = myexp,
  log_scale = TRUE,
  filter_data = TRUE,
  variable = "conc", 
  qc_types = "SPL",
  #x_lim = c(-6, 2),
  x_label = NA,
  feature_map = "lipidomics")


# Report with data, metadata and processing details
report_path <- "../output/mrmhub_report.xlsx"
mrmhub::save_report_xlsx(myexp,path = report_path)

csv_path <- "../output/processed_dataset.csv"
mrmhub::save_dataset_csv(
  data = myexp,
  path = csv_path,
  variable = "conc",
  qc_types = "SPL",
  include_qualifier = FALSE,
  filter_data = TRUE
)


# Sharing the MRMhubExperiment dataset
path = tempfile(fileext = ".rds")
saveRDS(myexp, file = path, compress = TRUE)
my_saved_exp <- readRDS(file = path)
print(myexp)




# MzQC export example
# CV terms from https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo, to be discussed

# Extract processed dataset
df <- myexp@dataset

# Define input file metadata
# TO-DO: check for fitting CV term
inputFiles = list(
  rmzqc::MzQCinputFile$new(
    name = "long.csv",
    location = paste0("file://", data_path),
    fileFormat = rmzqc::MzQCcvParameter$new(
      accession = "MS:1001040",
      name = "intermediate analysis format"
    )
  )
)

# Compute example QC metrics
# TO-DO: Discuss the metrics selection
qc_summary <- df %>%
  filter(qc_type %in% c("BQC", "TQC", "SPL")) %>%
  group_by(analysis_id) %>%
  summarise(
    cv = sd(feature_conc) /
    mean(feature_conc) * 100,
    median_rt = median(feature_rt),
    median_fwhm = median(feature_fwhm),
    .groups = "drop"
  )


# Convert each run into mzQC RunQuality
run_qualities <- lapply(1:nrow(qc_summary), function(i){
  
  row <- qc_summary[i, ]
  
  # Assign QC metrics to CV term, TO-DO: check for existing CV
  metrics <- list(
    rmzqc::MzQCqualityMetric$new(
      accession = "MS:1001883",
      name = "coefficient of variation",
      value = row$cv
    ),
    rmzqc::MzQCqualityMetric$new(
      accession = "MS:4000016",
      name = "retention time metric",
      value = row$median_rt
    ),
    rmzqc::MzQCqualityMetric$new(
      accession = "MS:1000086",
      name = "full width at half-maximum",
      value = row$median_fwhm
    )
  )
  
  # Build run quality object
  rmzqc::MzQCrunQuality$new(
    metadata = rmzqc::MzQCmetadata$new(
      label = as.character(row$analysis_id),
      inputFiles = inputFiles,
      analysisSoftware = list(
        rmzqc::MzQCanalysisSoftware$new(
          name = "MRMhub",
          version = as.character(
            packageVersion("mrmhub")
          ),
          accession = "MS:1000531",
          uri = "https://github.com/SLINGhub/MRMhub/"
        )
      )
    ),
    qualityMetrics = metrics
  )
})


# Create mzQC object
mzqc_obj <- rmzqc::MzQCmzQC$new(
  version = "1.0.0",
  creationDate = rmzqc::MzQCDateTime$new(Sys.time()),
  description = "mzQC export from MRMhub lipidomics workflow",
  runQualities = run_qualities
)


# Save file
output_path <- "../output/qc_results.mzQC"
writeMZQC(filepath = output_path, mzqc_obj = mzqc_obj)
cat("Yay! MzQC file successfully written to: ", output_path)
