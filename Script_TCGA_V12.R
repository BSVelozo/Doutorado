# --------------------------------------------------------------------------
# Script Principal para Análise TCGA-PRAD com TPM Oficial e Enriquecimento
# --------------------------------------------------------------------------

# Configuração inicial --------------------------------------------------------
work_dir <- "/home/bsvelozo/TCGA/TCGA_PRAD"
download_dir <- "/home/bsvelozo/TCGA/TCGA_downloads"

dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(download_dir, showWarnings = FALSE, recursive = TRUE)
setwd(work_dir)

options(timeout = 600)

# Carregar bibliotecas --------------------------------------------------------
required_packages <- c("dplyr", "magrittr", "ggplot2", "SummarizedExperiment", "DESeq2")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    if (pkg %in% c("SummarizedExperiment", "DESeq2")) {
      if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

if (!require("TCGAbiolinks", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("TCGAbiolinks")
}
library(TCGAbiolinks)

# 1. Query e download ---------------------------------------------------------
query <- GDCquery(
  project = "TCGA-PRAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  experimental.strategy = "RNA-Seq"
)

gdc_client_available <- system("which gdc-client", ignore.stderr = TRUE) == 0

tryCatch({
  if (gdc_client_available) {
    GDCdownload(query, directory = download_dir, method = "client", files.per.chunk = 3)
  } else {
    stop("Forçando uso do método API")
  }
}, error = function(e) {
  GDCdownload(query, directory = download_dir, method = "api", files.per.chunk = 2)
})

# 2. Preparar dados -----------------------------------------------------------
data <- GDCprepare(query, directory = download_dir, save = TRUE, 
                   save.filename = "TCGA-PRAD_data.rda")

# 3. Processar metadados e matrizes -------------------------------------------
cat("\n=== PROCESSANDO METADADOS E MATRIZES ===\n")
sample_info <- as.data.frame(colData(data))

# Criar dataframe básico
valid_samples <- data.frame(
  sample_id = rownames(sample_info),
  condition = ifelse(
    sample_info$definition == "Primary solid Tumor", "tumor",
    ifelse(sample_info$definition == "Solid Tissue Normal", "normal",
           ifelse(grepl("benign|hyperplasia|bph", sample_info$definition, 
                        ignore.case = TRUE), "bph", NA))
  ),
  row.names = rownames(sample_info),
  stringsAsFactors = FALSE
)

# Filtrar apenas amostras tumor e normal válidas para análise principal
valid_samples <- valid_samples[valid_samples$condition %in% c("tumor", "normal"), ]
count_matrix <- assay(data, "unstranded")
count_matrix <- count_matrix[, rownames(valid_samples)]

# Extrair matriz de TPM oficial (calculada pelo TCGA com comprimento dos genes)
if ("tpm_unstranded" %in% names(assays(data))) {
  tpm_matrix <- assay(data, "tpm_unstranded")
  tpm_matrix <- tpm_matrix[, rownames(valid_samples)]
  cat("Matriz de TPM oficial carregada (tpm_unstranded).\n")
} else {
  stop("Matriz de TPM não encontrada nos dados baixados. Verifique o objeto 'data'.")
}

# 4. Salvar dados intermediários ----------------------------------------------
write.table(count_matrix, "TCGA-PRAD_count_matrix.txt", sep = "\t", quote = FALSE)
write.csv(valid_samples, "TCGA-PRAD_sample_info.csv", row.names = TRUE)

# 5. Análise DESeq2 -----------------------------------------------------------
cat("\n=== INICIANDO ANÁLISE DESEQ2 ===\n")
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = valid_samples,
  design = ~ condition
)

# Pré-filtragem de genes com baixa expressão
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("Genes após pré-filtragem:", nrow(dds), "\n")

# Executar análise diferencial
dds <- DESeq(dds)

# Obter resultados com nível de significância de 5%
res <- results(dds, alpha = 0.05, contrast = c("condition", "tumor", "normal"))

# Filtrar genes significativos:
# - FDR (padj) <= 0.05
# - |log2FoldChange| >= 1
# - Sem NAs no FDR
resSig <- res[which(res$padj <= 0.05 & abs(res$log2FoldChange) >= 1 & !is.na(res$padj)), ]

# Ordenar pelo FDR ajustado (menor para maior)
resSig <- resSig[order(resSig$padj), ]

# 6. Salvar resultados da análise --------------------------------------------
write.csv(as.data.frame(res), "TCGA-PRAD_DESeq2_results_all.csv")
write.csv(as.data.frame(resSig), "TCGA-PRAD_DESeq2_results_significant.csv")

# 7. Calcular médias de TPM por condição (usando matriz oficial) -------------
cat("\n=== CALCULANDO MÉDIAS DE TPM (OFICIAL) ===\n")

# Separar amostras por condição
normal_samples <- rownames(valid_samples[valid_samples$condition == "normal", ])
tumor_samples <- rownames(valid_samples[valid_samples$condition == "tumor", ])

# Calcular médias de TPM por grupo
if (length(normal_samples) > 0) {
  tpm_normal <- rowMeans(tpm_matrix[, normal_samples, drop = FALSE], na.rm = TRUE)
} else {
  tpm_normal <- rep(NA, nrow(tpm_matrix))
  names(tpm_normal) <- rownames(tpm_matrix)
}

if (length(tumor_samples) > 0) {
  tpm_tumor <- rowMeans(tpm_matrix[, tumor_samples, drop = FALSE], na.rm = TRUE)
} else {
  tpm_tumor <- rep(NA, nrow(tpm_matrix))
  names(tpm_tumor) <- rownames(tpm_matrix)
}

# 8. Preparar tabelas de upregulated e downregulated com TPM oficial ---------

# Separar upregulated e downregulated
upregulated <- resSig[resSig$log2FoldChange >= 1, ]
downregulated <- resSig[resSig$log2FoldChange <= -1, ]

# Ordenar
upregulated <- upregulated[order(upregulated$log2FoldChange, decreasing = TRUE), ]
downregulated <- downregulated[order(downregulated$log2FoldChange), ]

# Função auxiliar para remover versão do Ensembl ID
strip_ensembl_version <- function(ids) {
  gsub("\\..*$", "", ids)
}

# Converter para dataframes e adicionar TPM (sem versão no rowname)
if (nrow(upregulated) > 0) {
  upregulated_df <- as.data.frame(upregulated)
  # Remover versão dos rownames para facilitar merge posterior
  rownames(upregulated_df) <- strip_ensembl_version(rownames(upregulated_df))
  upregulated_df$TPM_mean_normal <- tpm_normal[rownames(upregulated_df)]
  upregulated_df$TPM_mean_tumor <- tpm_tumor[rownames(upregulated_df)]
  write.csv(upregulated_df, "TCGA_DESeq2_upregulated.csv", row.names = TRUE)
}

if (nrow(downregulated) > 0) {
  downregulated_df <- as.data.frame(downregulated)
  rownames(downregulated_df) <- strip_ensembl_version(rownames(downregulated_df))
  downregulated_df$TPM_mean_normal <- tpm_normal[rownames(downregulated_df)]
  downregulated_df$TPM_mean_tumor <- tpm_tumor[rownames(downregulated_df)]
  write.csv(downregulated_df, "TCGA_DESeq2_downregulated.csv", row.names = TRUE)
}

# 9. Análise de Enriquecimento para Proteínas de Membrana e Interatoras -----
cat("\n=== ANÁLISE DE ENRIQUECIMENTO ===\n")

# Criar diretório para resultados de enriquecimento
enrich_dir <- file.path(work_dir, "Gene_Enrichment_Results")
dir.create(enrich_dir, showWarnings = FALSE, recursive = TRUE)

# Função para ler listas de genes (assumindo um ID por linha, sem cabeçalho)
read_gene_list <- function(file_path) {
  if (file.exists(file_path)) {
    genes <- readLines(file_path)
    genes <- strip_ensembl_version(genes)  # remover versão, se houver
    genes <- genes[genes != ""]
    return(genes)
  } else {
    warning("Arquivo não encontrado: ", file_path, ". Pulando.")
    return(character(0))
  }
}

# Carregar listas de proteínas
membrane_genes <- read_gene_list("ENSG_Membrane_Proteins.txt")
interactor_genes <- read_gene_list("ENSG_Interactors_Proteins.txt")

cat("Proteínas de membrana carregadas:", length(membrane_genes), "\n")
cat("Proteínas interatoras carregadas:", length(interactor_genes), "\n")

# Função para classificar enriquecimento com base no ratio (tumor/normal ou normal/tumor)
classify_enrichment <- function(ratio) {
  ifelse(ratio >= 10, "Enriquecido10X",
         ifelse(ratio >= 5, "Enriquecido5X",
                ifelse(ratio >= 2, "Enriquecido2X", NA)))
}

# Inicializar listas para armazenar contagens e dados detalhados
enrichment_summary <- list()
membrane_details <- data.frame()
interactor_details <- data.frame()

# ---- Análise para genes upregulados (cancerosos, tumor > normal) ----
if (nrow(upregulated_df) > 0) {
  # Calcular ratio tumor/normal
  up_ratio <- upregulated_df$TPM_mean_tumor / upregulated_df$TPM_mean_normal
  up_ratio[is.infinite(up_ratio)] <- 100  # se normal == 0, considerar >10X
  up_category <- classify_enrichment(up_ratio)
  
  # Adicionar ao dataframe
  upregulated_df$enrichment_category <- up_category
  upregulated_df$enriched_in <- "tumor"
  upregulated_df$tissue <- "prostate"
  
  # Contagens totais para cancerosos
  cancer_counts <- table(factor(up_category, levels = c("Enriquecido10X", "Enriquecido5X", "Enriquecido2X")))
  enrichment_summary$cancer_total <- cancer_counts
  enrichment_summary$cancer_total_genes <- nrow(upregulated_df[!is.na(up_category), ])
  
  # Subconjunto para proteínas de membrana
  up_membrane <- upregulated_df[rownames(upregulated_df) %in% membrane_genes, ]
  if (nrow(up_membrane) > 0) {
    mem_cancer_counts <- table(factor(up_membrane$enrichment_category, 
                                      levels = c("Enriquecido10X", "Enriquecido5X", "Enriquecido2X")))
    enrichment_summary$cancer_membrane <- mem_cancer_counts
    enrichment_summary$cancer_membrane_total <- nrow(up_membrane)
    membrane_details <- rbind(membrane_details, 
                              data.frame(ID_Ensembl = rownames(up_membrane),
                                         Tecido_enriquecido = "prostate",
                                         Categoria = up_membrane$enrichment_category,
                                         TPM = up_membrane$TPM_mean_tumor,
                                         Condicao = "canceroso",
                                         stringsAsFactors = FALSE))
  }
  
  # Subconjunto para proteínas interatoras
  up_interactor <- upregulated_df[rownames(upregulated_df) %in% interactor_genes, ]
  if (nrow(up_interactor) > 0) {
    int_cancer_counts <- table(factor(up_interactor$enrichment_category,
                                      levels = c("Enriquecido10X", "Enriquecido5X", "Enriquecido2X")))
    enrichment_summary$cancer_interactor <- int_cancer_counts
    enrichment_summary$cancer_interactor_total <- nrow(up_interactor)
    interactor_details <- rbind(interactor_details,
                                data.frame(ID_Ensembl = rownames(up_interactor),
                                           Tecido_enriquecido = "prostate",
                                           Categoria = up_interactor$enrichment_category,
                                           TPM = up_interactor$TPM_mean_tumor,
                                           Condicao = "canceroso",
                                           stringsAsFactors = FALSE))
  }
}

# ---- Análise para genes downregulados (saudáveis, normal > tumor) ----
if (nrow(downregulated_df) > 0) {
  # Calcular ratio normal/tumor
  down_ratio <- downregulated_df$TPM_mean_normal / downregulated_df$TPM_mean_tumor
  down_ratio[is.infinite(down_ratio)] <- 100
  down_category <- classify_enrichment(down_ratio)
  
  downregulated_df$enrichment_category <- down_category
  downregulated_df$enriched_in <- "normal"
  downregulated_df$tissue <- "prostate"
  
  # Contagens totais para saudáveis
  healthy_counts <- table(factor(down_category, levels = c("Enriquecido10X", "Enriquecido5X", "Enriquecido2X")))
  enrichment_summary$healthy_total <- healthy_counts
  enrichment_summary$healthy_total_genes <- nrow(downregulated_df[!is.na(down_category), ])
  
  # Subconjunto para proteínas de membrana
  down_membrane <- downregulated_df[rownames(downregulated_df) %in% membrane_genes, ]
  if (nrow(down_membrane) > 0) {
    mem_healthy_counts <- table(factor(down_membrane$enrichment_category,
                                       levels = c("Enriquecido10X", "Enriquecido5X", "Enriquecido2X")))
    enrichment_summary$healthy_membrane <- mem_healthy_counts
    enrichment_summary$healthy_membrane_total <- nrow(down_membrane)
    membrane_details <- rbind(membrane_details,
                              data.frame(ID_Ensembl = rownames(down_membrane),
                                         Tecido_enriquecido = "prostate",
                                         Categoria = down_membrane$enrichment_category,
                                         TPM = down_membrane$TPM_mean_normal,
                                         Condicao = "saudavel",
                                         stringsAsFactors = FALSE))
  }
  
  # Subconjunto para proteínas interatoras
  down_interactor <- downregulated_df[rownames(downregulated_df) %in% interactor_genes, ]
  if (nrow(down_interactor) > 0) {
    int_healthy_counts <- table(factor(down_interactor$enrichment_category,
                                       levels = c("Enriquecido10X", "Enriquecido5X", "Enriquecido2X")))
    enrichment_summary$healthy_interactor <- int_healthy_counts
    enrichment_summary$healthy_interactor_total <- nrow(down_interactor)
    interactor_details <- rbind(interactor_details,
                                data.frame(ID_Ensembl = rownames(down_interactor),
                                           Tecido_enriquecido = "prostate",
                                           Categoria = down_interactor$enrichment_category,
                                           TPM = down_interactor$TPM_mean_normal,
                                           Condicao = "saudavel",
                                           stringsAsFactors = FALSE))
  }
}

# 10. Escrever arquivo de sumário --------------------------------------------
cat("\n=== GERANDO ARQUIVOS DE ENRIQUECIMENTO ===\n")

summary_file <- file.path(enrich_dir, "Protein_Enrichment_Numbers.txt")
sink(summary_file)

cat("================================================================================\n")
cat("NÚMEROS DE ENRIQUECIMENTO DE PROTEÍNAS\n")
cat("================================================================================\n\n")

# Função para extrair contagem com segurança
get_count <- function(counts, level) {
  if (is.null(counts)) return(0)
  as.numeric(counts[level])
}

# Seção Saudável
cat("=== CLASSIFICAÇÃO DE ENRIQUECIMENTO SAUDÁVEL ===\n")
healthy_tot <- enrichment_summary$healthy_total_genes %||% 0
cat("Enriquecido10X:", get_count(enrichment_summary$healthy_total, "Enriquecido10X"), "genes\n")
cat("Enriquecido5X:", get_count(enrichment_summary$healthy_total, "Enriquecido5X"), "genes\n")
cat("Enriquecido2X:", get_count(enrichment_summary$healthy_total, "Enriquecido2X"), "genes\n")
cat("Total de genes enriquecidos (saudáveis):", healthy_tot, "\n\n")

cat("=== CLASSIFICAÇÃO DE ENRIQUECIMENTO CANCEROSO ===\n")
cancer_tot <- enrichment_summary$cancer_total_genes %||% 0
cat("Enriquecido10X:", get_count(enrichment_summary$cancer_total, "Enriquecido10X"), "genes\n")
cat("Enriquecido5X:", get_count(enrichment_summary$cancer_total, "Enriquecido5X"), "genes\n")
cat("Enriquecido2X:", get_count(enrichment_summary$cancer_total, "Enriquecido2X"), "genes\n")
cat("Total de genes enriquecidos (cancerosos):", cancer_tot, "\n\n")

cat("=== FILTROS APLICADOS ===\n\n")

cat("PROTEÍNAS DE MEMBRANA:\n")
cat("  Saudável:\n")
cat("    Enriquecido10X:", get_count(enrichment_summary$healthy_membrane, "Enriquecido10X"), "genes\n")
cat("    Enriquecido5X:", get_count(enrichment_summary$healthy_membrane, "Enriquecido5X"), "genes\n")
cat("    Enriquecido2X:", get_count(enrichment_summary$healthy_membrane, "Enriquecido2X"), "genes\n")
cat("    Total:", enrichment_summary$healthy_membrane_total %||% 0, "\n")
cat("  Canceroso:\n")
cat("    Enriquecido10X:", get_count(enrichment_summary$cancer_membrane, "Enriquecido10X"), "genes\n")
cat("    Enriquecido5X:", get_count(enrichment_summary$cancer_membrane, "Enriquecido5X"), "genes\n")
cat("    Enriquecido2X:", get_count(enrichment_summary$cancer_membrane, "Enriquecido2X"), "genes\n")
cat("    Total:", enrichment_summary$cancer_membrane_total %||% 0, "\n\n")

cat("PROTEÍNAS INTERACTORAS:\n")
cat("  Saudável:\n")
cat("    Enriquecido10X:", get_count(enrichment_summary$healthy_interactor, "Enriquecido10X"), "genes\n")
cat("    Enriquecido5X:", get_count(enrichment_summary$healthy_interactor, "Enriquecido5X"), "genes\n")
cat("    Enriquecido2X:", get_count(enrichment_summary$healthy_interactor, "Enriquecido2X"), "genes\n")
cat("    Total:", enrichment_summary$healthy_interactor_total %||% 0, "\n")
cat("  Canceroso:\n")
cat("    Enriquecido10X:", get_count(enrichment_summary$cancer_interactor, "Enriquecido10X"), "genes\n")
cat("    Enriquecido5X:", get_count(enrichment_summary$cancer_interactor, "Enriquecido5X"), "genes\n")
cat("    Enriquecido2X:", get_count(enrichment_summary$cancer_interactor, "Enriquecido2X"), "genes\n")
cat("    Total:", enrichment_summary$cancer_interactor_total %||% 0, "\n")

sink()

# 11. Escrever arquivos CSV detalhados ---------------------------------------
if (nrow(membrane_details) > 0) {
  write.csv(membrane_details, 
            file.path(enrich_dir, "all_enriched_membrane_proteins_TCGA.csv"),
            row.names = FALSE)
  cat("Arquivo de proteínas de membrana criado:", nrow(membrane_details), "linhas.\n")
} else {
  cat("Nenhuma proteína de membrana enriquecida encontrada.\n")
}

if (nrow(interactor_details) > 0) {
  write.csv(interactor_details,
            file.path(enrich_dir, "all_enriched_interactors_proteins_TCGA.csv"),
            row.names = FALSE)
  cat("Arquivo de proteínas interatoras criado:", nrow(interactor_details), "linhas.\n")
} else {
  cat("Nenhuma proteína interatora enriquecida encontrada.\n")
}

# 12. Visualização -----------------------------------------------------------
cat("\n=== GERANDO VISUALIZAÇÕES ===\n")
vsd <- vst(dds, blind = FALSE)

# Estatísticas das amostras no PCA
cat("=== ESTATÍSTICAS DAS AMOSTRAS NO PCA ===\n")
cat("Total de amostras processadas no PCA:", ncol(vsd), "\n")
cat("Distribuição por condição:\n")
amostras_pca <- table(colData(vsd)$condition)
print(amostras_pca)
cat("Normal:", sum(colData(vsd)$condition == "normal"), "amostras\n")
cat("Tumor:", sum(colData(vsd)$condition == "tumor"), "amostras\n")
cat("Proporção Tumor/Normal:", round(sum(colData(vsd)$condition == "tumor") / 
                                       sum(colData(vsd)$condition == "normal"), 2), "\n\n")

estatisticas_pca <- data.frame(
  Grupo = names(amostras_pca),
  Numero_Amostras = as.numeric(amostras_pca),
  Proporcao = round(as.numeric(amostras_pca) / ncol(vsd) * 100, 1)
)
write.csv(estatisticas_pca, "TCGA-PRAD_PCA_sample_stats.csv", row.names = FALSE)

# PCA plot
pca_title <- paste0("TCGA-PRAD - Tumor (", sum(colData(vsd)$condition == "tumor"), 
                    ") vs Normal (", sum(colData(vsd)$condition == "normal"), ")")
pca_plot <- plotPCA(vsd, intgroup = "condition") + 
  ggtitle(pca_title) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("TCGA-PRAD_PCA.png", plot = pca_plot, width = 8, height = 6)

# Volcano plot
volcano_data <- as.data.frame(res)
volcano_data$direction <- "ns"
volcano_data$direction[volcano_data$padj <= 0.05 & volcano_data$log2FoldChange >= 1 & !is.na(volcano_data$padj)] <- "up"
volcano_data$direction[volcano_data$padj <= 0.05 & volcano_data$log2FoldChange <= -1 & !is.na(volcano_data$padj)] <- "down"
volcano_data <- volcano_data[is.finite(-log10(volcano_data$padj)), ]

volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = subset(volcano_data, direction == "ns"), 
             color = "grey70", alpha = 0.6, size = 1.5) +
  geom_point(data = subset(volcano_data, direction == "down"), 
             color = "#1F78B4", alpha = 0.8, size = 2.5) +
  geom_point(data = subset(volcano_data, direction == "up"), 
             color = "#E31A1C", alpha = 0.8, size = 2.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.5) +
  labs(
    title = "FDR e p-value <= 0.05, |LFC| >= 1",
    x = bquote(~Log[2]~ "Fold Change"),
    y = bquote(~-Log[10]~ "FDR")
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
ggsave("TCGA-PRAD_volcano.png", plot = volcano_plot, width = 10, height = 8, dpi = 300)

# 13. Mensagem final ---------------------------------------------------------
cat("\n=== ANÁLISE CONCLUÍDA ===\n")
cat("Diretório de trabalho:", work_dir, "\n")
cat("Arquivos gerados:\n")
cat("- TCGA-PRAD_count_matrix.txt\n")
cat("- TCGA-PRAD_sample_info.csv\n")
cat("- TCGA-PRAD_DESeq2_results_all.csv\n")
cat("- TCGA-PRAD_DESeq2_results_significant.csv\n")
cat("- TCGA_DESeq2_upregulated.csv (com TPM médio oficial)\n")
cat("- TCGA_DESeq2_downregulated.csv (com TPM médio oficial)\n")
cat("- Gene_Enrichment_Results/Protein_Enrichment_Numbers.txt\n")
cat("- Gene_Enrichment_Results/all_enriched_membrane_proteins_TCGA.csv\n")
cat("- Gene_Enrichment_Results/all_enriched_interactors_proteins_TCGA.csv\n")
cat("- TCGA-PRAD_PCA_sample_stats.csv\n")
cat("- TCGA-PRAD_PCA.png\n")
cat("- TCGA-PRAD_volcano.png\n")

cat("\nEstatísticas finais:\n")
cat("  - Total de genes analisados:", nrow(res), "\n")
cat("  - Genes significativos (FDR ≤ 0.05, |LFC| ≥ 1):", nrow(resSig), "\n")
cat("  - Genes upregulated:", nrow(upregulated_df), "\n")
cat("  - Genes downregulated:", nrow(downregulated_df), "\n")
cat("  - Amostras totais:", ncol(vsd), "\n")
cat("  - Amostras tumor:", sum(colData(vsd)$condition == "tumor"), "\n")
cat("  - Amostras normal:", sum(colData(vsd)$condition == "normal"), "\n")

# Salvar sessão completa
save.image("TCGA-PRAD_analysis_complete.RData")
cat("\nSessão salva em: TCGA-PRAD_analysis_complete.RData\n")