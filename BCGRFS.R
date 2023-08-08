# 1.TCGA (Validation set) ####
library(TCGAbiolinks)
library(SummarizedExperiment)

setwd("C:\\Users\\hourg\\Desktop\\BLCA\\TCGA\\TCGA")

# 1-1.download TCGA data ####
BLCA = GDCquery(
  project = "TCGA-BLCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = "Primary Tumor"
)
GDCdownload(BLCA, method = "api")

# 1-2.integrate data ####
expr = GDCprepare(
  query = BLCA,
  save = TRUE,
  save.filename = "BLCA.RData",
  summarizedExperiment = TRUE
)

rowdata = rowData(expr)
names(rowdata)
table(rowdata$gene_type)

# 1-3.get gene expression ####
mrna = expr[rowdata$gene_type == "protein_coding", ]

# 1-4.RNAseq counts ####
TPM = as.data.frame(assay(mrna, i = "tpm_unstrand"))
FPKM = as.data.frame(assay(mrna, i = "fpkm_unstrand"))
Counts = as.data.frame(assay(mrna, i = "unstranded"))
# Counts i= "unstranded"
# tpm i="tpm_unstrand"
# fpkm i=" fpkm_unstrand"

save(mrna, TPM, FPKM, Counts, file = "expr_data.RData")

# 1-5. addgene symbol ####
load("expr_data.RData")

symbol_mrna = rowData(mrna)$gene_name
head(symbol_mrna)

expr_TPM = cbind(data.frame(symbol_mrna), as.data.frame(TPM))
# expr_counts=cbind(data.frame(symbol_mrna), as.data.frame(Counts))

# remove duliplicate genes
suppressPackageStartupMessages(library(tidyverse))
library(dplyr)

data = expr_TPM #expr_counts #expr_read
expr_Counts = data %>%
  as_tibble() %>%
  mutate(meanrow = rowMeans(.[, -1]), .before = 2) %>%
  arrange(desc(meanrow)) %>%
  distinct(symbol_mrna, .keep_all = T) %>%
  select(-meanrow) %>%
  column_to_rownames(var = "symbol_mrna") %>%
  as.data.frame()

save(expr_read, expr_Counts, clin, file = "exp_symbol.RData")
save(expr_read, expr_Counts, clin, file = "exp_symbol_TPM.RData")

# 1-6.import metadata ####
clin = GDCquery_clinic(project = "TCGA-BLCA",
                       type = "Clinical",
                       save.csv = TRUE)

# 1-7.merge exp & metadata ####
library(dplyr)
library(openxlsx)
load("exp_symbol_TPM.RData")

surv_RFS = read.xlsx("TCGA_RFS.xlsx")

surv$Treatment = NULL

percentile_threshold1 = 0.1  
threshold1 = quantile(expr_read, probs = percentile_threshold1, na.rm =
                        TRUE)
proportions1 = rowMeans(expr_read > threshold1)
percentage_threshold1 = 0.1 
g = rownames(expr_read)[proportions1 > percentage_threshold1]
dt = expr_read[g, ]

percentile_threshold2 = 0.9
threshold2 = quantile(dt, probs = percentile_threshold2, na.rm = TRUE)
proportions2 = rowMeans(dt < threshold2)
percentage_threshold2 = 0.9
g2 = rownames(dt)[proportions2 < percentage_threshold2]
expr_read = dt[g2, ]


data = expr_Counts #expr_read
exp = data.frame(t(data))
id = rownames(exp)
barcode = paste(surv[[1]])

result = data.frame()
for (i in barcode) {
  is_present = grepl(i, id)
  if (any(is_present)) {
    result = rbind(result, data.frame(barcode = i, id = id[is_present]))
  }
}

exp$id = rownames(exp)
dt = merge(result, exp, by = "id")

#tpm
dt = merge(dt, surv, by = "barcode")
dt$barcode = NULL
dt$id = gsub("-", "_", dt$id)

#count
dt$id = NULL
dt$barcode = gsub("-", "_", dt$barcode)

save(dt, file = "TCGA_tpm.RData")
save(dt, file = "TCGA_count.RData")

# 1-8. ####

# 2.GEO ####
library(GEOquery)
library(openxlsx)
setwd("C:\\Users\\hourg\\Desktop\\BLCA\\data\\exp_meta")

# 2-1.get GEO data ####
GSE = "GSE31684" # GSE19423
GPL = "GPL570"

gset = getGEO(GSE, GSEMatrix = TRUE, getGPL = TRUE)
if (length(gset) > 1){
  idx = grep(GPL, attr(gset, "names"))
}else{ idx = 1}

gset = gset[[idx]]
data = exprs(gset)

metadata <- gset@phenoData@data

write.xlsx(metadata, file = "metadata_GSE31684.xlsx")

# 2-2.covert probe ID ####
setwd("C:\\Users\\hourg\\Desktop\\BLCA\\data\\exp_meta")

GSE19423 = read.table("GSE19423_raw_data.txt",  header = TRUE, sep = "\t")
GSE31684 = read.table("GSE31684_raw_data.txt",  header = TRUE, sep = "\t")
GSE13507 = read.table("GSE13507_raw_data.txt",  header = TRUE, sep = "\t")

library(illuminaHumanv2.db) # hgu133plus2.db (GPL570)  illuminaHumanv2.db (GPL6102)
probe2symbol_df = toTable(get("illuminaHumanv2SYMBOL")) # hgu133plus2SYMBOL (GPL570)  illuminaHumanv2SYMBOL (GPL6102)
length(unique(probe2symbol_df$probe_id))
length(unique(probe2symbol_df$symbol))
data = as.data.frame(GSE31684)

data = data[-1, ]
rownames(data) = data$PROBE_ID
data$PROBE_ID = NULL

# GSE13507
GSE13507 = GSE13507[-1, ]
rownames(GSE13507) = GSE13507[, 1]
GSE13507[, 1] = NULL
#colnames(GSE13507)[1]="PROBE_ID"

data = as.data.frame(GSE13507)
for (col in names(data)) {
  data[[col]] <- as.numeric(data[[col]])
}

library(dplyr)
library(tibble)
data = data %>%
  rownames_to_column("probe_id") %>%
  inner_join(probe2symbol_df, by = "probe_id") %>%
  dplyr::select(-probe_id) %>%
  dplyr::select(symbol, everything()) %>%
  mutate(rowMean = rowMeans(.[, -1])) %>%序
  arrange(desc(rowMean)) %>%
  distinct(symbol, .keep_all = T) %>%
  dplyr::select(-rowMean) %>%
  column_to_rownames("symbol")


percentile_threshold1 = 0.1  # 
threshold1 = quantile(data, probs = percentile_threshold1, na.rm = TRUE)
proportions1 = rowMeans(data > threshold1)
percentage_threshold1 = 0.1  # 
g = rownames(data)[proportions1 > percentage_threshold1]
dt = data[g, ]

percentile_threshold2 = 0.9
threshold2 = quantile(dt, probs = percentile_threshold2, na.rm = TRUE)
proportions2 = rowMeans(dt < threshold2)
percentage_threshold2 = 0.9
g2 = rownames(dt)[proportions2 < percentage_threshold2]
dt = data[g2, ]

GSE31684 = data.frame(t(dt))
GSE13507 = data.frame(t(data))

# 2-3.merge gene experssion id and metadata barcode ####
#used when it is inconsistent with the barcode and ID, such as GSE19423
meta = read.xlsx("metadata_GSE19423.xlsx", sheet = 2)
barcode = rownames(GSE19423)
id = paste(meta[[1]])

result = data.frame()
for (i in id) {
  is_present = grepl(i, barcode)
  if (any(is_present)) {
    result = rbind(result, data.frame(id = i, barcode = barcode[is_present]))
  }
}

GSE19423$barcode = rownames(GSE19423)
dt = merge(result, GSE19423, by = "barcode")

dt$barcode = NULL

rownames(dt) = dt$id
dt$id = NULL

GSE19423 = dt

#used when it is inconsistent with the barcode and ID, such as GSE31684
meta = read.xlsx("metadata_GSE31684.xlsx",
                 sheet = 2,
                 rowNames = T)
dt = merge(meta, GSE31684, by = 0)
colnames(dt)[1] = "id"

GSE31684 = dt


save(GSE31684, meta, file = "GSE31684.RData")

# 2-4.Uni-cox regression of OS/RFS (GSE31684) ####
library(survival)
library(data.table)

load("data/GSE31684.RData") #GSE31684
dt = GSE31684

uni_cox = function(dt) {
  coxY = c("RFS_Status", "RFS_Time") #RFS_Status: 0(no recurrence) 1(recurrence)
  colnames = colnames(dt)
  allGenes = setdiff(colnames, c(coxY, "id", "OS_Status", "OS_Time"))
  outcome = NULL

  cat("unicox......................................")

  for (gene in allGenes) {
    #debug:
    #gene=allGenes[9887]
    #..............

    tryCatch({
      fm = formula(paste("Surv(RFS_Time, RFS_Status)~", gene))
      fm

      cox = coxph(fm, data = dt)
      coxSummary = summary(cox)

      tmp = data.frame(
        cbind(
          gene,
          HR = coxSummary$conf.int[, "exp(coef)"],
          HR.95L = coxSummary$conf.int[, "lower .95"],
          HR.95H = coxSummary$conf.int[, "upper .95"],
          P.value = coxSummary$coefficients[, "Pr(>|z|)"]
        )
      )
      tmp

      outcome = rbind(outcome, tmp)
      coxY = c(coxY, gene)

      # fwrite(outcome, file="Uni.csv", append=T)
      # write.csv(outcome,"Uni.csv")

    },

    error = function(e) {
      message('An Error Occurred')
    },

    warning = function(w) {
      message('A Warning Occurred')
    })


  }

  cat("filter sig genes......................................")

  outcome$P.adjusted = p.adjust(outcome$P.value, "fdr")
  uni = outcome
  uni$HR = as.numeric(uni$HR)
  uni$HR = sprintf("%.3f", uni$HR)

  # uni_adjust=filter(uni, P.adjusted<0.05)
  write.xlsx(uni, 'uni_sig.xlsx')

  # coxYa=c("id", "OS_Status", "OS_Time")
  # coxYb=c(coxYa,uni_adjust$gene)
  # sigGenes=dt[,coxYb]

  return(uni)
}
uni = uni_cox(dt)


# 2-5.DEGs of grading (GSE19423) ####
library(DESeq2)

load("data/GSE19423.RData")
countData = data.frame(t(GSE19423))
countData = round(countData)

colData = meta

countData = countData[rowMeans(countData) > 1, ]
dds = DESeqDataSetFromMatrix(countData = countData,
                             colData = colData,
                             design =  ~ Grade)
head(dds)

# rm read<10
keep = rowSums(counts(dds)) >= 10
dds = dds[keep, ]

dds$Grade = relevel(dds$Grade, ref = "low")
dds1 = DESeq(
  dds,
  fitType = 'mean',
  minReplicatesForReplace = 7,
  parallel = FALSE
)

res = results(dds1)
summary(res)

res1 = data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
res1 = res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)),]

res1_up = res1[which(res1$log2FoldChange >= 1 &
                       res1$pvalue < 0.05), ]
res1_down = res1[which(res1$log2FoldChange <= -1 &
                         res1$pvalue < 0.05), ]
res1_total = rbind(res1_up, res1_down)

save(res1_up, res1_down, res1_total, file = "DEGs_grading.RData")

# 2-6.Training set (GSE154261) ####
GSE154261 = read.table("GSE154261_raw_data.txt", header = TRUE, sep = "\t")
GSE154261 = filter(GSE154261, gene_biotype == "protein_coding")
GSE154261$gene_name[which(GSE154261$gene_name == "HSPA14")[2]] = "MSANTD7"

rownames(GSE154261) = GSE154261$gene_name
GSE154261 = GSE154261[-(1:7)]
data = data.frame(lapply(GSE154261, as.numeric))
rownames(data) = rownames(GSE154261)

# raw counts to TPM
library(edgeR)
dge = DGEList(counts = data)
dge = calcNormFactors(dge, method = "TMM")
dge = estimateCommonDisp(dge)
expr = data.frame(cpm(dge, normalized.lib.sizes = TRUE) * 1e6)

save(expr, file = "GSE154261.RData")

# 3.filter genes ####
# 3-1.remove genes with expression median absolute deviation (MAD) ≤0.01 and at the bottom 25% genes ####
load("C:\\Users\\hourg\\Desktop\\BLCA\\data\\GSE154261.RData")
dt = data.frame(t(expr))

my_mad = function(x) {
  mad(x, na.rm = TRUE)
} 
m.mad = apply(dt, 2, my_mad)

dt = dt[, which(m.mad > max(quantile(m.mad, probs = seq(0, 1, 0.25))[2], 0.01))]

# 3-2.rm genes with undetectable expression in >75% of samples ####
detection_rate = apply(dt, 1, function(x)
  sum(x != 0) / length(x))
undetectable_genes = detection_rate < 0.75
dt = dt[!undetectable_genes,]

save(dt, file = "GSE154261_del.RData")

# 4.gene union ####
library(ggVennDiagram)
library(VennDiagram)
library(openxlsx)
setwd("C:\\Users\\hourg\\Desktop\\BLCA\\result")

# Mutation, CNV, Structure variant
Mutation = read.xlsx("bioPortal.xlsx", sheet = 1)
threshold1 = quantile(Mutation$Freq, 0.95)
filtered_Mutation = subset(Mutation, Freq >= threshold1)

CNV = read.xlsx("bioPortal.xlsx", sheet = 2)
threshold2 = quantile(CNV$Freq, 0.95)
filtered_CNV = subset(CNV, Freq >= threshold2)

SV = read.xlsx("bioPortal.xlsx", sheet = 3)
threshold3 = quantile(SV$Freq, 0.95)
filtered_SV = subset(SV, Freq >= threshold3)

venn_bio = list(
  "Mutation" = filtered_Mutation[, 1],
  "Copy number variant" = filtered_CNV[, 1],
  "Structure variation" = filtered_SV[, 1]
)
union_bio = data.frame(union(union(venn_bio[["Mutation"]], venn_bio[["Copy number variant"]]), venn_bio[["Structure variation"]])) #聯集

ggVennDiagram(
  venn_bio,
  edge_lty = "solid",
  edge_size = 3,
  set_size = 8,
  label = "both",
  label_color = "black",
  label_size = 5,
  label_alpha = 0
) + 
  scale_fill_gradient(low = "white", high = "#b9292b", name = "gene count") +
  scale_x_continuous(expand = expansion(mult = .4))

# GSE19423, GSE31684
OS = read.xlsx("uni_DEGs_genes.xlsx", sheet = 1)
RFS = read.xlsx("uni_DEGs_genes.xlsx", sheet = 2)
Grading = read.xlsx("uni_DEGs_genes.xlsx", sheet = 5)

venn = list(OS = OS[, 1],
            RFS = RFS[, 1],
            Grading = Grading[, 1])
union = data.frame(union(union(venn[["OS"]], venn[["RFS"]]), venn[["Grading"]])) #union

ggVennDiagram(
  venn,
  edge_lty = "solid",
  edge_size = 3,
  set_size = 8,
  label = "both",
  label_color = "black",
  label_size = 5,
  label_alpha = 0
) + 
  scale_fill_gradient(low = "white", high = "#b9292b", name = "gene count") +
  scale_x_continuous(expand = expansion(mult = .4))

# Total
union_total = data.frame(union(union[, 1], union_bio[, 1])) 

venn_total = list("Mutation and CNV and SV" = union_bio[, 1],
                  "Genes of OS and RFS and Grading" = union[, 1])
ggVennDiagram(
  venn_total,
  edge_lty = "solid",
  edge_size = 3,
  set_size = 5,
  label = "both",
  label_color = "black",
  label_size = 5,
  label_alpha = 0
) + 
  scale_fill_gradient(low = "white", high = "#b9292b", name = "gene count") +
  scale_x_continuous(expand = expansion(mult = .4))

save(union_total, file = "union.RData")

# GSE154261
load("GSE154261_del.RData")
load("C:\\Users\\hourg\\Desktop\\BLCA\\result\\union.RData")

expr = data.frame(t(dt))
colnames(union_total)[1] = "gene"
expr$gene = rownames(expr)
dt = merge(expr, union_total, by = "gene")
rownames(dt) = dt$gene
dt = data.frame(t(dt[-1]))

save(dt, file = "GSE154261_union.RData")

# 5.classify high/low cohorts ####
library(Rtsne)
library(ggplot2)

load("data\\GSE154261.RData") #expr (19937*73)
load("TCGA\\TCGA_tpm.RData")  #dt (34*19942)
expr = data.frame(t(dt))

gene = list(
  Ta = c("FGFR3", "KDM6A", "PIK3CA", "STAG2", "CREBBP", "KMT2C"),
  T1 = c(
    "KDM6A",
    "FGFR3",
    "KMT2D",
    "ARID1A",
    "ERCC2",
    "TP53",
    "ELF3",
    "CREBBP"
  )
)
union = data.frame(union(gene[["Ta"]], gene[["T1"]])) 

genes = union[, 1]
selected_data = expr[genes, ]

transposed_data = data.frame(t(selected_data))
for (col_name in genes) {
  transposed_data[[col_name]] <-
    as.numeric(transposed_data[[col_name]])
}

# tSNE
set.seed(222)
n = 7 #Training
n = 10 #Validation
tsne_result = Rtsne(
  transposed_data,
  dims = 2,
  perplexity = n,
  verbose = TRUE
)
tsne_coords = tsne_result$Y

hc.norm = hclust(dist(tsne_coords))
tmp = factor(cutree(hc.norm, 2))

tsne_coords = data.frame(tsne_coords)
tsne_coords$cluster = tmp
tsne_coords$cluster = factor(tsne_coords$cluster,
                             levels = c(1, 2),
                             labels = c("High", "Low"))
rownames(tsne_coords) = rownames(transposed_data)

tsne_plot = ggplot(tsne_coords, aes(x = X1, y = X2, color = cluster)) +
  geom_point(size = 7) +
  theme_bw() +
  labs(title = "tSNE") + xlab("") + ylab("") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 30),
    axis.text.x = element_text(face = "bold", size = 20),
    axis.text.y = element_text(face = "bold", size = 20),
    axis.title.x = element_text(hjust = 0.5, face = "bold", size = 25),
    axis.title.y = element_text(
      hjust = 0.5,
      angle = 90,
      face = "bold",
      size = 25
    ),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 18)
  ) +
  scale_color_discrete(name = "Immune score")
tsne_plot

##tsne_coords
#saveRDS(tsne_coords,"data/GSE154261_high_low_classification.RData")

ggsave(
  "result/Immune score/tSNE_Training set.png",
  tsne_plot,
  width = 12,
  height = 12,
  dpi = 600
)
ggsave(
  "result/Immune score/tSNE_Validation set.png",
  tsne_plot,
  width = 12,
  height = 12,
  dpi = 600
)

save(tsne_result, tsne_coords, file = "result/Immune score/tSNE_Training set.RData")
save(tsne_result, tsne_coords, file = "result/Immune score/tSNE_Validation set.RData")

# TCGA
tsne_coords$id = dt$id
dt = merge(dt, tsne_coords, by = "id")
dt = dt[-(42:43)] #34*19943

save(dt, file = "TCGA\\TCGA_tpm.RData")

# 6.compute ESTIMATE ####
library(estimate)
library(edgeR)

load("data\\GSE154261.RData") #expr


write.table(
  expr,
  file = "data\\Estimate\\geneExpSet.txt",
  sep = "\t",
  quote = F,
  row.names = T
)

# to gct file
filterCommonGenes(input.f = "data\\Estimate\\geneExpSet.txt",
                  output.f = "data\\Estimate\\geneExpSet.gct",
                  id = "GeneSymbol")

estimateScore(input.ds = "data\\Estimate\\geneExpSet.gct",
              output.ds = "data\\Estimate\\estimate_score.gct",
              platform = "illumina")

scores = read.table("data\\Estimate\\estimate_score.gct",
                    skip = 2,
                    header = T)
rownames(scores) = scores[, 1]
scores = data.frame(t(scores[, 3:ncol(scores)]))

TumourPurity = cos(0.6049872018 + 0.0001467884 * scores$ESTIMATEScore)
scores$TumourPurity = TumourPurity

scores

saveRDS(scores, "data\\Estimate\\scores.RData") #me

write.table(scores,
            file = "data\\Estimate\\estimate_score.tsv",
            sep = "\t",
            quote = F)


library(datat.table)

if (!file.exists("data\\Estimate\\estimated_purity_plots")) {
  dir.create(
    "data\\Estimate\\estimated_purity_plots",
    showWarnings = FALSE,
    recursive = TRUE
  )
}


# 7.WGCNA ####
library(WGCNA)
library(preprocessCore)
library(doParallel)
library(fastcluster)
library(dynamicTreeCut)
library(openxlsx)

setwd("C:\\Users\\hourg\\Desktop\\BLCA\\result\\WGCNA\\不分組")
setwd("C:\\Users\\hourg\\Desktop\\BLCA\\result\\WGCNA\\High")
setwd("C:\\Users\\hourg\\Desktop\\BLCA\\result\\WGCNA\\Low")

load("result\\GSE154261_union.RData") #dt 73*2360
meta = read.xlsx("data\\exp_meta\\metadata_GSE154261.xlsx", sheet = 2)
meta
# id RFS_Status RFS_Time
# 1  JMT1_84          1    72.00
# 2  JMT1_32          1    16.57
# 3  JMT1_26          1     3.00
# 4  JMT1_81          0    60.00
# 5   JMT1_1          0    56.20

ESTIMATE = read.xlsx("data\\exp_meta\\metadata_GSE154261.xlsx", sheet =
                       3)

options(stringsAsFactors = FALSE) 

# 7-1.metadata ####
metadata = merge(ESTIMATE, meta, by = "id")
rownames(metadata) = metadata$id

metadata = metadata[1:6]
metadata$cluster = recode(metadata$cluster, "Low" = 0, "High" = 1)
colnames(metadata)[3] = "Stromal Score"
colnames(metadata)[4] = "Immune Score"
colnames(metadata)[5] = "ESTIMATE Score"
colnames(metadata)[6] = "Tumour Purity"


meta = metadata

meta_low = filter(metadata, cluster == 1)
meta_high = filter(metadata, cluster == 0)

meta = meta_high
meta = meta_low
meta$cluster = NULL

dim(meta)
names(meta)

# 7-2.genotype ####
exp = dt


gp_high = meta[, 1]
gp_low = meta[, 1]

exp = dt[gp_high, ]
exp = dt[gp_low, ]

dim(exp)
head(exp)

# check missing and outliers
gsg = goodSamplesGenes(exp, verbose = 3)
gsg$allOK

sampleTree = hclust(dist(exp), method = "average")
sizeGrWindow(20, 12) 
par(cex = 0.8)
par(mar = c(6, 10, 6, 0)) 
plot(
  sampleTree,
  main = "Sample clustering to detect outliers",
  sub = "",
  xlab = "",
  cex.lab = 1.5,
  cex.axis = 1,
  cex.main = 3
)

# rm outliers
height = 1.07e+10  # high 1.1e+10 , low 9.2e+9 , 不分組 1.07e+10
abline(h = height, col = "red")
clust = cutreeStatic(sampleTree, cutHeight = height, minSize = 10)
table(clust)
keepSamples = (clust == 1)  
exp = exp[keepSamples,]
nGenes = ncol(exp)
nSamples = nrow(exp)


Samples = rownames(exp) 
dataRows = match(Samples, meta$id)
Metadata = meta[dataRows,-1] 
rownames(Metadata) = meta[dataRows, 1]
collectGarbage()

sampleTree2 = hclust(dist(exp), method = "average")
clinicalColors = numbers2colors(Metadata, signed = FALSE)
plotDendroAndColors(
  sampleTree2,
  clinicalColors,
  groupLabels = names(Metadata),
  main = "Sample dendrogram and heatmap",
  cex.axis = 0.5,
  cex.dendroLabels = 0.1
)

save(exp, Metadata, file = "BLCA_01_dataInput.RData")

# 7-3.WGCNA ####
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
lnames = load(file = "BLCA_01_dataInput.RData")
lnames

# 7-3-1 ####
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(exp, powerVector = powers, verbose = 5)

sizeGrWindow(12, 3) 

png(
  "02_power.png",
  width = 3000,
  height = 1500,
  units = "px",
  bg = "white",
  res = 216
)
par(mfrow = c(1, 2))
cex1 = 0.9

# Scale independence
plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit,signed R^2",
  type = "n",
  main = paste("Scale independence")
)
text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers,
  cex = cex1,
  col = "red"
)
abline(h = 0.90, col = "red") 

# Mean connectivity
plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = paste("Mean connectivity")
)
text(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = powers,
  cex = cex1,
  col = "red"
)
dev.off()


sft$powerEstimate

# 7-3-2. ####
cor = WGCNA::cor # 
net = blockwiseModules(
  exp,
  power = sft$powerEstimate,
  TOMType = "unsigned",
  minModuleSize = 50,
  #High 50，Low 80
  reassignThreshold = 0,
  mergeCutHeight = 0.2,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "BCTOM",
  verbose = 3
)
table(net$colors) 
cor = stats::cor 

sizeGrWindow(12, 9)
png(
  "03_cluster dendogram.png",
  width = 1500,
  height = 1000,
  units = "px",
  bg = "white",
  res = 216
)
mergedColors = labels2colors(net$colors)
plotDendroAndColors(
  net$dendrograms[[1]],
  mergedColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
dev.off()

table(mergedColors)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs

geneTree = net$dendrograms[[1]]

save(MEs, moduleLabels, moduleColors, geneTree,
     file = "BLCA_02_networkConstruction.RData")

# 7-4 ####
library(WGCNA)
options(stringsAsFactors = FALSE)

lnames = load(file = "BLCA_01_dataInput.RData")
lnames = load(file = "BLCA_02_networkConstruction.RData")


nGenes = ncol(exp)
nSamples = nrow(exp)

MEs0 = moduleEigengenes(exp, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, Metadata, use = "p")

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
names(MEs)

sizeGrWindow(25, 15)

textMatrix = paste(signif(moduleTraitCor, 2),
                   "\n(",
                   signif(moduleTraitPvalue, 1),
                   ")",
                   sep = "")

dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(0, 10, 2, 0)) 


labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(Metadata),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  setStdMargins = T,
  textMatrix = textMatrix,
  textAdj = c(0.5, 0.5),
  cex.text = 0.8,
  cex.lab.x = 0.8,
  cex.lab.y = 0.8,
  xLabelsAngle = 90,
  yColorWidth = 0.03,
  zlim = c(-1, 1),
  yColorOffset = 0.01,
  main = paste("Module-trait relationships")
)

# 7-5. ####
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(as.numeric(str_replace_all(
  colnames(MEs), "ME", ""
))))
MEs_col = orderMEs(MEs_col)

plotEigengeneNetworks(
  MEs_col,
  "Eigengene adjacency heatmap",
  marDendro = c(3, 3, 2, 8),
  #下、左、上、右
  marHeatmap = c(3, 4, 2, 2),
  #下、左、上、右
  plotDendrograms = T,
  xLabelsAngle = 90
)

# 7-6.hub genes ####
datKME = signedKME(exp, MEs, outputColumnName = "kME_MM.")

hubgenes = function(datKME) {
  FilterGenes = abs(module) > 0.8
  hubgene = dimnames(data.frame(exp))[[2]][FilterGenes]
}

library(data.table)
total = NULL

module = datKME$kME_MM.red
total$red = hubgenes(datKME)

module = datKME$kME_MM.yellow
total$yellow = hubgenes(datKME)

module = datKME$kME_MM.brown
total$brown = hubgenes(datKME)

module = datKME$kME_MM.turquoise
total$turquoise = hubgenes(datKME)

total_hubgenes = data.frame(unlist(total, use.names = FALSE))
names(total_hubgenes)[1] = "gene"

save(MEs,
     moduleLabels,
     moduleColors,
     geneTree,
     total_hubgenes,
     file = "BLCA_02_networkConstruction.RData")

# 7-7.####
module_genes = data.frame(Gene = colnames(exp), Module = moduleColors)
gene_list = split(module_genes$Gene, module_genes$Module)
module_brown_genes = gene_list[[3]]
module_red_genes = gene_list[[7]]
module_yellow_genes = gene_list[[8]]

table(module_brown_genes %in% module_red_genes)
table(module_yellow_genes %in% module_brown_genes)

# 7-8.Enrichment ####
setwd("C:\\Users\\hourg\\Desktop\\BLCA\\result\\WGCNA\\High")
lnames = load(file = "BLCA_01_dataInput.RData")
lnames = load(file = "BLCA_02_networkConstruction.RData")


moduletotalgene = function(modColors, modNames) {
  moduleGenes = modColors == modNames
  module_total_gene = as.data.frame(dimnames(data.frame(exp))[[2]][moduleGenes])
  module_total_gene = data.frame(gene = module_total_gene[, 1])
  return(module_total_gene)
}

genes = NULL
modNames = substring(names(MEs), 3)
gene_list = lapply(modNames, moduletotalgene, modColors = moduleColors)

gene1 = data.frame(do.call(cbind, gene_list[[1]]))
gene1$module = "green"

gene2 = data.frame(do.call(cbind, gene_list[[2]]))
gene2$module = "brown"

gene3 = data.frame(do.call(cbind, gene_list[[3]]))
gene3$module = "red"

gene4 = data.frame(do.call(cbind, gene_list[[4]]))
gene4$module = "black"

gene5 = data.frame(do.call(cbind, gene_list[[5]]))
gene5$module = "blue"

gene6 = data.frame(do.call(cbind, gene_list[[6]]))
gene6$module = "turquoise"

gene7 = data.frame(do.call(cbind, gene_list[[7]]))
gene7$module = "yellow"

gene = do.call(rbind, list(gene1, gene2, gene3, gene4, gene5, gene6, gene7))

# Annotation
a = gene$gene
b = gene$module
Enrichment = userListEnrichment(
  a,
  b,
  useBrainLists = F,
  useBloodAtlases = F,
  omitCategories = F,
  outputCorrectedPvalues = TRUE,
  outputGenes = F,
  minGenesInCategory = 1,
  useStemCellLists = F,
  useBrainRegionMarkers = F,
  useImmunePathwayLists = F,
  usePalazzoloWang = T,
  nameOut = "Enrichment.csv"
)


# 8.Cox regression ####
library(openxlsx)
# setwd("C:\\Users\\hourg\\Desktop\\BLCA\\result\\Cox")

# Training set
load("result\\WGCNA\\High\\BLCA_01_dataInput.RData")
load("result\\WGCNA\\High\\BLCA_02_networkConstruction.RData")

load("result\\WGCNA\\Low\\BLCA_01_dataInput.RData")
load("result\\WGCNA\\Low\\BLCA_02_networkConstruction.RData")

load("result\\WGCNA\\不分組\\BLCA_01_dataInput.RData")
load("result\\WGCNA\\不分組\\BLCA_02_networkConstruction.RData")

metadata = read.xlsx("data\\exp_meta\\metadata_GSE154261.xlsx",
                     sheet = 2,
                     rowNames = T)

# 8-1-1.Training set ####

total_hubgenes = total_hubgenes[, 1]  
train_dt = exp[, total_hubgenes]
train_dt = merge(metadata, train_dt, by = 0)
colnames(train_dt)[1] = "id"

dim(train_dt) #69*27

# High/Low immune
# genes=total_hubgenes[,1]
genes = total_hubgenes
train_dt = exp[, genes]

train_dt = merge(metadata, train_dt, by = 0)
colnames(train_dt)[1] = "id"

# 8-1-2.Validation set ####
# load("C:\\Users\\hourg\\Desktop\\BLCA\\TCGA\\TCGA\\TCGA_tpm.RData")
load("TCGA\\TCGA_tpm.RData") #dt
table(dt$cluster) #23 high 11 low
test_dt = dt[-(19940:19942)] #234*19939

# High/Low
# dt=dt[-(19940:19941)]
test_dt = filter(dt, cluster == "High")
test_dt = filter(dt, cluster == "Low")
test_dt$cluster = NULL

save(train_dt, test_dt, file = "Cox_data_high.RData")
save(train_dt, test_dt, file = "Cox_data_low.RData")
save(train_dt, test_dt, file = "Cox_data_nongrouped.RData")

# 8-2.Lasso & Multi cox ####
library(glmnet)
library(survival)
library(pROC)

# setwd("C:\\Users\\hourg\\Desktop\\BLCA\\result\\Cox")
load("result\\Cox\\high\\Cox_data_high.RData") #test_dt (23*19939) train_dt (35*37)
load("result\\Cox\\low\\Cox_data_low.RData") #test_dt (11*19939) train_dt (33*37)
load("result\\Cox\\Cox_data_nongrouped.RData") #test_dt (103*19939) train_dt(69*37)

train = train_dt[4:ncol(train_dt)]
train = log2(train_dt[4:ncol(train_dt)] + 1)
# boxplot(t(train))
# boxplot(train)
train = cbind(train_dt[1:3], train)

test = test_dt[2:19937]
test = log2(test[1:ncol(test)] + 1)
# boxplot(t(test))
# boxplot(test[1:10])
test = cbind(test_dt[, 1], test_dt[19938:19939], test)

# load("result\\Cox\\treated.RData") 

multi = function(train, test, pFilter = 0.05) {
  x = train[-(1:3)]
  y = Surv(train$RFS_Time, train$RFS_Status)
  cat("lasso regression...................")
  lasso_fit = glmnet(
    data.matrix(x),
    as.matrix(y),
    family = "cox",
    alpha = 1,
    type.measure = 'deviance',
    maxit = 100000
  )
  cvfit = cv.glmnet(
    data.matrix(x),
    as.matrix(y),
    family = "cox",
    alpha = 1,
    type.measure = 'deviance'
  )

  png(
    "lasso.png",
    width = 3000,
    height = 1500,
    units = "px",
    bg = "white",
    res = 216
  )
  par(mfrow = c(1, 2))
  plot(lasso_fit,
       xvar = "lambda",
       label = TRUE,
       cex.lab = 1.2)
  plot(cvfit, cex.lab = 1.2)
  abline(v = log(c(cvfit$lambda.min, cvfit$lambda.1se)), lty = "dashed")
  dev.off()

  coef = coef(lasso_fit, s = cvfit$lambda.min)
  index = which(coef != 0)
  actCoef = coef[index]

  coxY = c("RFS_Time", "RFS_Status")
  lasso_gene = row.names(coef)[index]
  coxY_lasso = c(coxY, lasso_gene)
  train = train[, coxY_lasso]

  # multi
  coxY = c("RFS_Time", "RFS_Status", "id")
  colnames = colnames(train) #train
  genes = setdiff(colnames, coxY)

  fm = as.formula(paste0(
    "Surv(RFS_Time, RFS_Status)~",
    paste(genes, sep = "", collapse = "+")
  ))

  Cox = coxph(fm, data = train) #train
  Cox2 = step(Cox, direction = "both")

  Cox2_fm = Cox2$formula
  Cox2_fm
  Cox2_pre = coxph(Cox2_fm, data = train) #train

  coxSummary = summary(Cox2_pre)
  tmp = data.frame(
    cbind(
      coef = coxSummary$coefficients[, "coef"],
      HR = coxSummary$conf.int[, "exp(coef)"],
      HR.95L = coxSummary$conf.int[, "lower .95"],
      HR.95H = coxSummary$conf.int[, "upper .95"],
      P.value = coxSummary$coefficients[, "Pr(>|z|)"]
    )
  )
  tmp

  tmp$HR = sprintf("%.3f", tmp$HR)
  tmp$HR.95L = sprintf("%.3f", tmp$HR.95L)
  tmp$HR.95H = sprintf("%.3f", tmp$HR.95H)
  tmp$"95CI%" = paste(tmp[, 3], tmp[, 4], sep = "-")
  tmp$P.value = sprintf("%.3f", tmp$P.value)
  tmp

  write.csv(tmp, "g_pvalue.csv")

  # ROC
  #Cox2_fm=as.formula(Surv(RFS_Time, RFS_Status)~CSGALNACT1+FAM135B+MAP1A+PDE2A+PRRX1+SYT11+TRPA1)
  #Cox2_pre=coxph(Cox2_fm, data=train)

  train_predict = predict(Cox2_pre, type = 'risk', newdata = train)
  train_roc = roc(
    as.numeric(train$RFS_Status),
    train_predict,
    levels = c(0, 1),
    direction = "<"
  )
  train_AUC = as.numeric(auc(train_roc))

  test_predict = predict(Cox2_pre, type = 'risk', newdata = test) #test
  test_roc = roc(
    as.numeric(test$RFS_Status),
    test_predict,
    levels = c(0, 1),
    direction = "<"
  )
  test_AUC = as.numeric(auc(test_roc))

  formula_to_string = function(Cox2_fm) {
    if (is.character(Cox2_fm)) {
      return(Cox2_fm)
    }

    parts = as.character(Cox2_fm)
    complete_formula_str = paste(parts, collapse = " ")

    return(complete_formula_str)
  }
  formula = formula_to_string(Cox2_fm)

  result = data.frame(formula, train_AUC, test_AUC)
  return(result)

}
#取得train test sets的最佳組合基因AUC
multicox = multi(train, test, pFilter = 0.05)

save(train, test, Cox2, tmp, result, file = "result\\Cox\\high\\Cox_high_result.RData")
save(train, test, Cox2, tmp, result, file = "result\\Cox\\low\\Cox_low_result.RData")
save(train, test, Cox2, tmp, result, file = "result\\Cox\\不分組\\Cox_nongrouped_result.RData")

# 9.figures ####
# 9-1.ROC curve ####
training_plot = plot(
  train_roc,
  print.auc = TRUE,
  auc.polygon = TRUE,
  grid = c(0.1, 0.2),
  grid.col = c("green", "red"),
  max.auc.polygon = TRUE,
  print.thres = TRUE,
  main = "ROC curve for the Training set (Genes and Variates)"
)
test_plot = plot(
  test_roc,
  print.auc = TRUE,
  auc.polygon = TRUE,
  grid = c(0.1, 0.2),
  grid.col = c("green", "red"),
  max.auc.polygon = TRUE,
  print.thres = TRUE,
  main = "ROC curve for the Validation (Genes and Variates)"
)

png(
  "roc_.png",
  width = 1500,
  height = 1500,
  units = "px",
  bg = "white",
  res = 216
)
plot(
  smooth(training_plot),
  col = "orange",
  legacy.axes = T,
  main = "ROC curve of low immune score"
)
plot((test_plot), add = TRUE, col = "cadetblue2")
legend(
  "bottomright",
  c("Training set (AUC) 0.944",
    "Validation set (AUC) 0.056"),
  cex = 1,
  text.col = c("orange", "cadetblue2"),
  col = c("orange", "cadetblue2")
)
dev.off()

# 9-2.KM plot ####
library(survival)
library(survminer)
library(cutpointr)
library(cowplot)
library(dplyr)

setwd("C:\\Users\\hourg\\Desktop\\BLCA\\result\\Cox")
load("C:\\Users\\hourg\\Desktop\\BLCA\\result\\Cox\\high\\Cox_high_result.RData")
load("C:\\Users\\hourg\\Desktop\\BLCA\\result\\Cox\\low\\Cox_low_result.RData")
load(
  "C:\\Users\\hourg\\Desktop\\BLCA\\result\\Cox\\不分組\\Cox_nongrouped_result.RData"
)

fm = as.formula(Surv(RFS_Time, RFS_Status) ~ ACKR1 + CACNA2D1 + PTPRT +
                  SELE)
fm = as.formula(Surv(RFS_Time, RFS_Status) ~ CSGALNACT1 + FAM135B + MAP1A +
                  PDE2A + PRRX1 + SYT11 + TRPA1)
fm = as.formula(
  Surv(RFS_Time, RFS_Status) ~ ATP8B2 + CD48 + CD84 + CELF2 + DOCK10 + FCRL3 +
    HCK + IL2RA + PYHIN1 + RCSD1 + SLAMF1 + SLAMF6 + SP140
)

# Training
pre = coxph(fm, data = train)
train_predict = predict(pre, type = "risk", newdata = train)
train_predict = as.numeric(train_predict)
summary(train_predict)

train$category = cut(
  train_predict,
  breaks = c(-Inf, 363.549, Inf),
  labels = c("low", "high")
)
table(train$category)

fit_training = survfit(Surv(RFS_Time, RFS_Status) ~ category, data = train)
summary(fit_training)

#P value & HR
data.survdiff = survdiff(Surv(RFS_Time, RFS_Status) ~ category, data = train)
p.val_train = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR_train = (data.survdiff$obs[2] / data.survdiff$exp[2]) / (data.survdiff$obs[1] /
                                                              data.survdiff$exp[1])
train_index = data.frame(p.val_train, HR_train)

training = ggsurvplot(
  fit_training,
  pval = F,
  font.title = 25,
  font.tickslab = 10,
  font.x = 15,
  font.y = 15,
  legend = c(0.86, 0.2),
  legend.title = "Risk score",
  legend.labs = c("Low", "High"),
  title = "Training set",
  ggtheme = theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
)

# Validation
pre1 = coxph(fm, data = test)
test_predict = predict(pre1, type = 'risk', newdata = test)
test_predict = as.numeric(test_predict)
summary(test_predict)

test$category = cut(
  test_predict,
  breaks = c(-Inf, 363.549, Inf),
  labels = c("high", "low")
)
table(test$category)

fit_test = survfit(Surv(RFS_Time, RFS_Status) ~ category, data = test)

#P value & HR
data.survdiff1 = survdiff(Surv(RFS_Time, RFS_Status) ~ category, data =
                            test)
p.val_test = 1 - pchisq(data.survdiff1$chisq, length(data.survdiff1$n) - 1)
HR_test = (data.survdiff1$obs[2] / data.survdiff1$exp[2]) / (data.survdiff1$obs[1] /
                                                               data.survdiff1$exp[1])
test_index = data.frame(p.val_test, HR_test)

testing = ggsurvplot(
  fit_test,
  pval = F,
  font.title = 25,
  font.tickslab = 10,
  font.x = 15,
  font.y = 15,
  legend = c(0.86, 0.2),
  legend.title = "Risk score",
  legend.labs = c("High", "Low"),
  title = "Validation set",
  ggtheme = theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
)

KM = plot_grid(training$plot, testing$plot, ncol = 2)
ggsave("KM.png",
       KM,
       width = 12,
       height = 6,
       dpi = 300)

# 9-3.Nomogram ####
library(survival)
library(rms)
setwd("C:\\Users\\hourg\\Desktop\\BLCA\\result\\Cox")
load("C:\\Users\\hourg\\Desktop\\BLCA\\result\\Cox\\high\\Cox_high_result.RData")
load("C:\\Users\\hourg\\Desktop\\BLCA\\result\\Cox\\low\\Cox_low_result.RData")
load(
  "C:\\Users\\hourg\\Desktop\\BLCA\\result\\Cox\\不分組\\Cox_nongrouped_result.RData"
)

dd <- datadist(train)
options(datadist = 'dd')

# Nomogram
# Training set
fm = as.formula(Surv(RFS_Time, RFS_Status) ~ ACKR1 + CACNA2D1 + PTPRT +
                  SELE)
fm = as.formula(Surv(RFS_Time, RFS_Status) ~ CSGALNACT1 + FAM135B + MAP1A +
                  PDE2A + PRRX1 + SYT11 + TRPA1)

coxfit = cph(
  fm,
  data = train,
  x = T,
  y = T,
  surv = T
)

surv = Survival(coxfit)
surv1 = function(x)
  surv(12, x) # 1yr
surv3 = function(x)
  surv(36, x) # 3yr
surv5 = function(x)
  surv(60, x) # 5yr

nom = nomogram(
  coxfit,
  fun = list(surv1, surv3, surv5),
  lp = F,
  funlabel = c(
    "1-year survival Probability",
    "3-year survival Probability",
    "5-year survival Probability"
  ),
  maxscale = 100,
  fun.at = c(0.99, 0.8, 0.6, 0.4, 0.2, 0.05)
)

png(
  "Nomogram_Training.png",
  width = 3000,
  height = 1500,
  units = "px",
  bg = "white",
  res = 216
)
plot(
  nom,
  label.every = 2,
  cex.axis = 0.5,
  xfrac = 0.3,
  col.grid = gray(c(0.8, 0.95))
)
dev.off()

# validate Nomogram
coxfit = cph(
  fm,
  data = train,
  x = T,
  y = T,
  surv = T
)
rcorrcens(Surv(RFS_Time, RFS_Status) ~ predict(coxfit), data = train)
cindex_train = 1 - 0.069

coxfit1 = cph(
  fm,
  data = test,
  x = T,
  y = T,
  surv = T
)
rcorrcens(Surv(RFS_Time, RFS_Status) ~ predict(coxfit1), data = test)
cindex_valid = 1 - 0.268

Nomogram_cindex = data.frame(cindex_train, cindex_valid)

# Validation set
dd1 <- datadist(test)
options(datadist = 'dd1')

fm = as.formula(Surv(RFS_Time, RFS_Status) ~ ACKR1 + CACNA2D1 + PTPRT +
                  SELE)
fm = as.formula(Surv(RFS_Time, RFS_Status) ~ CSGALNACT1 + FAM135B + MAP1A +
                  PDE2A + PRRX1 + SYT11 + TRPA1)
coxfit1 = cph(
  fm,
  data = test,
  x = T,
  y = T,
  surv = T
)

surv_T = Survival(coxfit1)
surv1_1 = function(x)
  surv(12, x) # 1yr OS
surv3_1 = function(x)
  surv(36, x) # 3yr OS
surv5_1 = function(x)
  surv(60, x) # 5yr OS

nom1 = nomogram(
  coxfit1,
  fun = list(surv1_1, surv3_1, surv5_1),
  lp = F,
  funlabel = c(
    "1-year survival Probability",
    "3-year survival Probability",
    "5-year survival Probability"
  ),
  maxscale = 100,
  fun.at = c(0.99, 0.8, 0.6, 0.4, 0.2, 0.05)
)

png(
  "Nomogram_Validation.png",
  width = 3000,
  height = 1500,
  units = "px",
  bg = "white",
  res = 216
)
plot(
  nom1,
  label.every = 2,
  cex.axis = 0.7,
  xfrac = 0.3,
  col.grid = gray(c(0.8, 0.95))
)
dev.off()

# 9-4.Calibration ####
library(survival)
library(rms)
# setwd("C:\\Users\\hourg\\Desktop\\BLCA\\result\\Cox")
# load("C:\\Users\\hourg\\Desktop\\BLCA\\result\\Cox\\high\\Cox_high_result.RData")
# load("C:\\Users\\hourg\\Desktop\\BLCA\\result\\Cox\\low\\Cox_low_result.RData")
# load("C:\\Users\\hourg\\Desktop\\BLCA\\result\\Cox\\不分組\\Cox_nongrouped_result.RData")

load("result\\Cox\\high\\Cox_high_result.RData")
load("result\\Cox\\low\\Cox_low_result.RData")
load("result\\Cox\\不分組\\Cox_nongrouped_result.RData")



dd <- datadist(train)
options(datadist = 'dd')

dv <- datadist(test)
options(datadist = 'dv')

fm = as.formula(Surv(RFS_Time, RFS_Status) ~ ACKR1 + CACNA2D1 + PTPRT +
                  SELE)
fm = as.formula(Surv(RFS_Time, RFS_Status) ~ CSGALNACT1 + FAM135B + MAP1A +
                  PDE2A + PRRX1 + SYT11 + TRPA1)

# 1 Year
model1 = cph(
  fm,
  data = train,
  x = TRUE,
  y = TRUE,
  surv = TRUE,
  time.inc = 12
)
cal1 = calibrate(
  model1,
  u = 5,
  cmethod = 'KM',
  m = 100,
  B = 1000
)

model1_v = cph(
  fm,
  data = test,
  x = TRUE,
  y = TRUE,
  surv = TRUE,
  time.inc = 12
)
cal1_v = calibrate(
  model1_v,
  u = 5,
  cmethod = 'KM',
  m = 20,
  B = 1000
)

png(
  "1-Year Calibration.png",
  width = 1500,
  height = 1500,
  units = "px",
  bg = "white",
  res = 216
)
par(mar = c(6, 4, 1, 1), cex = 0.9)
plot(
  cal1,
  lty = 1,
  lwd = 2.5,
  xlim = c(0.2, 1),
  ylim = c(0.2, 1),
  errbar.col = "#00AFBB",
  col = "#00AFBB",
  xlab = "Nomogram-Predicted Probability of 1-Year RFS",
  ylab = "Observed RFS (%)",
  mgp = c(2, 1, 0)
)
plot(
  cal1_v,
  add = T,
  lty = 1,
  lwd = 2.5,
  conf.int = T,
  subtitles = F,
  cex.subtitles = 0.8,
  errbar.col = "#FC4E07",
  col = "#FC4E07"
)
legend(
  "bottomright",
  legend = c("Training set", "Validation"),
  col = c("#00AFBB", "#FC4E07"),
  lwd = 2
)
abline(0, 1, lty = 3, lwd = 1, col = "grey")
dev.off()

# 3 Year
model3 = cph(
  fm,
  data = train,
  x = TRUE,
  y = TRUE,
  surv = TRUE,
  time.inc = 36
)
cal3 = calibrate(
  model3,
  u = 5,
  cmethod = 'KM',
  m = 100,
  B = 1000
)

model3_v = cph(
  fm,
  data = test,
  x = TRUE,
  y = TRUE,
  surv = TRUE,
  time.inc = 36
)
cal3_v = calibrate(
  model3_v,
  u = 5,
  cmethod = 'KM',
  m = 20,
  B = 1000
)

png(
  "3-Year Calibration.png",
  width = 1500,
  height = 1500,
  units = "px",
  bg = "white",
  res = 216
)
par(mar = c(6, 4, 1, 1), cex = 0.9)
plot(
  cal3,
  lty = 1,
  lwd = 2.5,
  xlim = c(0.2, 1),
  ylim = c(0.2, 1),
  errbar.col = "#00AFBB",
  col = "#00AFBB",
  xlab = "Nomogram-Predicted Probability of 3-Year RFS",
  ylab = "Observed RFS (%)",
  mgp = c(2, 1, 0)
)
plot(
  cal3_v,
  add = T,
  lty = 1,
  lwd = 2.5,
  conf.int = T,
  subtitles = F,
  cex.subtitles = 0.8,
  errbar.col = "#FC4E07",
  col = "#FC4E07"
)
legend(
  "bottomright",
  legend = c("Training set", "Validation"),
  col = c("#00AFBB", "#FC4E07"),
  lwd = 2
)
abline(0, 1, lty = 3, lwd = 1, col = "grey")
dev.off()

# 5 Year
model5 = cph(
  fm,
  data = train,
  x = TRUE,
  y = TRUE,
  surv = TRUE,
  time.inc = 60
)
cal5 = calibrate(
  model5,
  u = 5,
  cmethod = 'KM',
  m = 100,
  B = 1000
)

model5_v = cph(
  fm,
  data = test,
  x = TRUE,
  y = TRUE,
  surv = TRUE,
  time.inc = 60
)
cal5_v = calibrate(
  model5_v,
  u = 5,
  cmethod = 'KM',
  m = 20,
  B = 1000
)

png(
  "5-Year Calibration.png",
  width = 1500,
  height = 1500,
  units = "px",
  bg = "white",
  res = 216
)
par(mar = c(6, 4, 1, 1), cex = 0.9)
plot(
  cal5,
  lty = 1,
  lwd = 2.5,
  xlim = c(0.2, 1),
  ylim = c(0.2, 1),
  errbar.col = "#00AFBB",
  col = "#00AFBB",
  xlab = "Nomogram-Predicted Probability of 5-Year RFS",
  ylab = "Observed RFS (%)",
  mgp = c(2, 1, 0)
)
plot(
  cal5_v,
  add = T,
  lty = 1,
  lwd = 2.5,
  conf.int = T,
  subtitles = F,
  cex.subtitles = 0.8,
  errbar.col = "#FC4E07",
  col = "#FC4E07"
)
legend(
  "bottomright",
  legend = c("Training set", "Validation"),
  col = c("#00AFBB", "#FC4E07"),
  lwd = 2
)
abline(0, 1, lty = 3, lwd = 1, col = "grey")
dev.off()


#10-correlation matrix and ggpairs----------------------------------
load("result\\GSE154261_union.RData") #dt 73*2360
geneLow = c('CSGALNACT1',
            'FAM135B',
            'MAP1A',
            'PDE2A',
            'PRRX1',
            'SYT11',
            'TRPA1')

geneHigh = c('ACKR1',
             'CACNA2D1',
             'PTPRT',
             'SELE')

dt2 = dt[, c(geneLow, geneHigh)]
dim(dt2)

#Import ESTIMATE results
scores = readRDS("data\\ESTIMATE\\scores.RData")
sum(row.names(scores) != row.names(dt2))

dtscores = cbind(dt2, scores)

saveRDS(dtscores, "data\\ESTIMATE\\dtscores.RData")


#correlation matrix plot
dtscores = readRDS(file = "data\\ESTIMATE\\dtscores.RData")


library(corrplot)
correlation_matrix <- cor(dtscores)
corrplot(correlation_matrix, method = "circle")
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
p.mat <- cor.mtest(dtscores) #p values

corrplot(
  correlation_matrix,
  type = "upper",
  col = brewer.pal(n = 8, name = "RdYlBu"),
  p.mat = p.mat,
  sig.level = 0.05
)

#11-heatmap of significant genes, Estimates, Recurrence######################
##GSE154261 and ESTIMATE

dtscores = readRDS("data\\Estimate\\dtscores.RData")
dim(dtscores) #73*15
head(dtscores)

#survival data
meta = read.xlsx("data\\exp_meta\\metadata_GSE154261.xlsx", sheet = 2)
meta

#high low immune classified
immunC = readRDS("data/GSE154261_high_low_classification.RData")


#match case id
sort(row.names(dtscores))
sort(meta$id)


dtscmt = merge(dtscores, meta, by.x = 0, by.y = "id")
dtscmt = merge(dtscmt, immunC, by.x = "Row.names", by.y = 0)
dtscmt
dim(dtscmt) #73 *21
dtscmt$RFS_Status = as.factor(dtscmt$RFS_Status)


dtscmt$cluster = car::recode(dtscmt$cluster, "'High'='Low';'Low'='High'")
table(dtscmt$cluster)
# High  Low
# 36   37
saveRDS(dtscmt, file = "data/dtscmt_GSE154261_exper_Estimate_survival_data.RData")


#grouped boxplot-----------
dtscmt = readRDS(file = "data/dtscmt_GSE154261_exper_Estimate_survival_data.RData")
library(ggplot2)

# group histogram--------------
ggplot(dtscmt, aes(x = CSGALNACT1, fill = cluster)) +
  geom_histogram(position = "identity",
                 alpha = 0.7,
                 bins = 20) +
  facet_grid(cluster ~ .) +
  labs(title = "Histogram by Group", x = "Value", y = "Frequency") +
  theme_minimal()

# Create a grouped boxplot-----------------
##cluster vs immune/stromal index

a <-
  ggplot(dtscmt, aes(x = cluster, y = TumourPurity, fill = cluster)) +
  geom_boxplot() +
  labs(y = "Tumor purity", x = "") +
  theme_minimal()

# Add p-value comparisons to the plot
a <- a + stat_compare_means(
  comparisons = list(c("High", "Low")),
  method = "t.test",
  label.y = max(dtscmt$TumourPurity)
)


b <-
  ggplot(dtscmt, aes(x = cluster, y = StromalScore, fill = cluster)) +
  geom_boxplot() +
  labs(y = "Stromal Score", x = "") +
  theme_minimal()

# Add p-value comparisons to the plot
b <- b + stat_compare_means(
  comparisons = list(c("High", "Low")),
  method = "t.test",
  #"t.test"
  label.y = max(dtscmt$StromalScore)
)



c <-
  ggplot(dtscmt, aes(x = cluster, y = ImmuneScore, fill = cluster)) +
  geom_boxplot() +
  labs(y = "Immune Score", x = "") +
  theme_minimal()

# Add p-value comparisons to the plot
c <- c + stat_compare_means(
  comparisons = list(c("High", "Low")),
  method = "t.test",
  label.y = max(dtscmt$ImmuneScore)
)


d <-
  ggplot(dtscmt, aes(x = cluster, y = ESTIMATEScore, fill = cluster)) +
  geom_boxplot() +
  labs(y = "ESTIMATE Score", x = "") +
  theme_minimal()

# Add p-value comparisons to the plot
d <- d + stat_compare_means(
  comparisons = list(c("High", "Low")),
  method = "t.test",
  label.y = max(dtscmt$ESTIMATEScore)
)

library(ggpubr)
ggarrange(a, b, c, d + rremove("x.text"),
          # labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)

fN = file.path(
  "result",
  paste0(
    today(),
    "box plot of immune groups with stromal_immune_tumorpurity.png"
  )
)
fN

ggsave(fN)


# Create a grouped scatter plot using ggplot2-------------

ggplot(dtscmt, aes(x = ImmuneScore , y = StromalScore, color = cluster)) +
  geom_point() +
  labs(
    title = "Grouped Scatter Plot",
    x = "X-axis",
    y = "Y-axis",
    color = "Group"
  ) +
  theme_minimal()


ggplot(dtscmt, aes(x = ImmuneScore , y = ESTIMATEScore, color = cluster)) +
  geom_point() +
  labs(
    title = "Grouped Scatter Plot",
    x = "X-axis",
    y = "Y-axis",
    color = "Group"
  ) +
  theme_minimal()

ggplot(dtscmt, aes(x = ImmuneScore , y = TumourPurity, color = cluster)) +
  geom_point() +
  labs(
    title = "Grouped Scatter Plot",
    x = "X-axis",
    y = "Y-axis",
    color = "Group"
  ) +
  theme_minimal()

ggplot(dtscmt, aes(x = StromalScore , y = TumourPurity, color = cluster)) +
  geom_point() +
  labs(
    title = "Grouped Scatter Plot",
    x = "X-axis",
    y = "Y-axis",
    color = "Group"
  ) +
  theme_minimal()


ggplot(dtscmt, aes(x = StromalScore , y = ESTIMATEScore, color = cluster)) +
  geom_point() +
  labs(
    title = "Grouped Scatter Plot",
    x = "X-axis",
    y = "Y-axis",
    color = "Group"
  ) +
  theme_minimal()


ggplot(dtscmt, aes(x = TumourPurity , y = ESTIMATEScore, color = cluster)) +
  geom_point() +
  labs(
    title = "Grouped Scatter Plot",
    x = "X-axis",
    y = "Y-axis",
    color = "Group"
  ) +
  theme_minimal()




#ggpairs by group--------------------
library(GGally)
library(RColorBrewer)
library(lubridate)

indicators = c("StromalScore",
               "ImmuneScore",
               "ESTIMATEScore",
               "TumourPurity")
genes = c(
  "CSGALNACT1",
  "FAM135B",
  "MAP1A",
  "PDE2A" ,
  "PRRX1",
  "SYT11",
  "TRPA1",
  "ACKR1",
  "CACNA2D1",
  "PTPRT",
  "SELE"
)
gp = c("RFS_Status", "cluster")


for (i in seq_along(genes)) {
  ##group by two variables
  p <- ggpairs(dtscmt[, c(genes[i], indicators, gp)],
               mapping = aes(color = RFS_Status , shape = cluster))
  ##group by one variables
  p <- ggpairs(dtscmt[, c(genes[i], indicators, gp)],
               mapping = aes(color = RFS_Status))

  fN = paste0(today(),
              " ggpiar of ",
              genes[i],
              " by group_",
              paste(gp, collapse = "."),
              ".pdf")
  fN = file.path("result", fN)
  fN

  ggsave(fN, p, width = 12, height = 12)

}


#correlation matrix group by cluster----------------------
library(corrplot)

genesH = c("ACKR1", "CACNA2D1", "PTPRT", "SELE")

genesL = c("CSGALNACT1",
           "FAM135B",
           "MAP1A",
           "PDE2A" ,
           "PRRX1",
           "SYT11",
           "TRPA1")

dtscoresH = dtscmt[dtscmt$cluster == "High", c(genesH, indicators)]
dtscoresL = dtscmt[dtscmt$cluster == "Low", c(genesL, indicators)]

dtscores = dtscoresH
dtscores = dtscoresL

correlation_matrix <- cor(dtscores)
corrplot(correlation_matrix, method = "circle")
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
p.mat <- cor.mtest(dtscores) #p values

corrplot(
  correlation_matrix,
  type = "upper",
  col = brewer.pal(n = 8, name = "RdYlBu"),
  p.mat = p.mat,
  sig.level = 0.05
)

# Load the required library
library(survival)

# Example usage of coxph
result <- coxph(Surv(RFS_Time, RFS_Status) ~ cluster,
                id = Row.names,
                data = dtscmt)


# 12- Cox regression of high and low immune -------------------------------
library(survival)
library(ggplot2)

dt = readRDS(file = "data/dtscmt_GSE154261_exper_Estimate_survival_data.RData")
dim(dt) #73*21
table(dt$cluster, dt$RFS_Status)
# 0  1
# High 23 14
# Low  19 17

require("survival")
library("survminer")

#covert 0/1 to 1/2 and numeric format
dt$RFS_Status = as.numeric(car::recode(dt$RFS_Status, "1=2;0=1"))
table(dt$cluster, dt$RFS_Status)
str(dt)

fit <- survfit(Surv(RFS_Time, RFS_Status) ~ cluster , data = dt)
# Basic survival curves
ggsurvplot(fit, data = dt)

fit <-
  survfit(Surv(RFS_Time, RFS_Status) ~ CSGALNACT1 + strata(cluster) , data = dt)
# Basic survival curves
ggsurvplot(fit, data = dt)



aa_fit <-
  survival::aareg(Surv(RFS_Time, RFS_Status) ~ cluster , data = dt)
ggplot2::autoplot(aa_fit)


# cox_model <- coxph(Surv(RFS_Time, RFS_Status) ~ cluster, data = dt,id=Row.names)
summary(cox_model)

predicted_surv <- survfit(cox_model$formula, data = dt)


# Plot the survival curves by group
ggplot(predicted_surv, aes(x = RFS_Time, y = RFS_Status, color = cluster)) +
  geom_step() +
  labs(x = "Time", y = "Survival Probability", color = "Group") +
  theme_minimal()

