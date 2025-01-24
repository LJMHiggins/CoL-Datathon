library(TCGAbiolinks)


# Download TNBC BRCA TCGA subset info -------------------------------------

# Article:
# Proteogenomic analysis of triple-negative breast cancer
# identifies subtype-specific therapeutic vulnerabilities and
# epigenetic immune suppression

# This article inferred TNBC using TNBCtype https://cbc.app.vumc.org/tnbc/
# by evaluating distribution of ER, PGR and HER2 using RNA, protein
# (RPPA/mass spec) and DNA copy number annotated with clinical
# assessment (IHC and FISH) provided by TCGA

tnbc <- readxl::read_xlsx(
  "Table_S2_TNBCsubtype clinical information and signatures.xlsx", # downloaded from https://github.com/TransBioInfoLab/TNBC_analysis/blob/cbfe410e4df1fd35cdd97bda548f176362844c08/data/Table_S2_TNBCsubtype%20clinical%20information%20and%20signatures.xlsx
  sheet = "A-TCGA_TNBC_subtype"
)

tnbc_patient <- tnbc$patient


# Download deconvoluted cell types ----------------------------------------

# The Immune Landscape of Cancer Immunity, 2019, https://www.sciencedirect.com/science/article/pii/S1074761318301213?via%3Dihub#app2

# This article deconvoluted cell types abundance in TCGA sampels
tcga_deconvoluted <- readxl::read_xlsx("1-s2.0-S1074761318301213-mmc2.xlsx")

# filter out TNBC
tcga_deconvoluted <- tcga_deconvoluted[tcga_deconvoluted$`TCGA Participant Barcode` %in% tnbc_patient, ]
data.table::fwrite(tcga_deconvoluted, "TCGA_TNBC_tcga_deconvoluted.csv")

# Retrieve TCGA Clinical data ---------------------------------------------

# Download data
query_clin <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Clinical",
  data.format = "BCR XML",
  barcode = tnbc_patient
)

GDCdownload(query_clin)

# Prepare patient data
patient_data <- GDCprepare_clinic(
  query = query_clin,
  clinical.info = "patient"
)
data.table::fwrite(patient_data, "TCGA_TNBC_patient.csv")


# Prepare treatment data
treatment_data <- GDCprepare_clinic(
  query_clin,
  clinical.info = "drug"
)
data.table::fwrite(treatment_data, "TCGA_TNBC_treatment.csv")


# Prepare follow up data
fu_data <- GDCprepare_clinic(
  query_clin,
  clinical.info = "follow_up"
)
data.table::fwrite(fu_data, "TCGA_TNBC_fu_data.csv")


# Retrieve TCGA RNAseq data -----------------------------------------------

# Query RNAseq data
query_rna <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  platform = "Illumina",
  barcode = tnbc_patient,
  experimental.strategy = "RNA-Seq"
)

# Download RNAseq data
GDCdownload(query_rna)

# Prepare RNA data
BRCA_exp <- GDCprepare(
  query = query_rna
)
save(BRCA_exp, file = "TCGA_TNBC_RNAseq_TPMs.RData")

tpm <- BRCA_exp@assays@data$tpm_unstrand
rownames(tpm) <- rownames(BRCA_exp)
colnames(tpm) <- colnames(BRCA_exp)

tpm_hugo <- misomix::ens2hugo(t(tpm))

tpm_hugo <- as.data.frame(tpm_hugo)
data.table::fwrite(t(tpm_hugo[, -which(duplicated(colnames(tpm_hugo)))]), "TCGA_TNBC_RNAseq_TPMs.tsv", sep = "\t", row.names = TRUE)

colData_df <- as.data.frame(BRCA_exp@colData)
list_columns <- sapply(colData_df, is.list)
for (col in names(list_columns[list_columns])) {
  colData_df[[col]] <- sapply(colData_df[[col]], function(x) paste(x, collapse = ", "))
}
data.table::fwrite(colData_df, "TCGA_TNBC_RNAseq_Metadata.tsv", sep = "\t")
