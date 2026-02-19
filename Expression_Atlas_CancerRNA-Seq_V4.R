# Script: Expression_Atlas_RNAseq_Only.R
# Versão: 2.3 - Com correções para classificação de prostate/prostate_gland
# Descrição: Analisa experimentos RNA-seq do Expression Atlas, classificando amostras de próstata em Cancer/Saudavel/HPB
# Data: 2024
# Autor: Adaptado para projeto de análise de expressão gênica

# =============================================================================
# SEÇÃO 1: CARREGAMENTO DE PACOTES E BIBLIOTECAS
# =============================================================================

library(readr)      # Leitura e escrita de arquivos CSV
library(tibble)     # Estrutura de dados tibble
library(purrr)      # Programação funcional
library(tidyr)      # Manipulação de dados (pivotagem)
library(dplyr)      # Manipulação de dados
library(SummarizedExperiment) # Para dados de expressão
library(qs)         # Leitura de arquivos .qs
library(S4Vectors)  # Estruturas de dados S4

# =============================================================================
# SEÇÃO 2: CONFIGURAÇÕES E PARÂMETROS
# =============================================================================

# Diretórios de entrada e saída
experiment_dir <- "/home/bsvelozo/Expression_Atlas/Experiment_Files/Projeto_Cancer_Prostata"
gene_expression_dir <- "/home/bsvelozo/Expression_Atlas/Expression_Files_Cancer"

# Arquivos de configuração
genes_file <- "All_IDs_ENSG.txt"           # Lista de genes de interesse
log_file <- "gene_analysis_log.txt"         # Arquivo de log
not_analyzed_file <- "NotAnalyzed.txt"      # Lista de experimentos não analisados
classification_file <- "experimentos_categorizados.csv"  # Classificação

# Criar diretório de saída se não existir
dir.create(gene_expression_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# SEÇÃO 3: FUNÇÕES AUXILIARES E DE LOG
# =============================================================================

#' Função de logging para registrar mensagens
#' @param msg Mensagem a ser registrada
log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  full_msg <- paste0("[", timestamp, "] ", msg)
  cat(full_msg, "\n")
  write(full_msg, file = log_file, append = TRUE)
}

#' Carrega a lista de genes de interesse do arquivo
#' @return Vetor com IDs dos genes
load_target_genes <- function() {
  if (!file.exists(genes_file)) {
    log_message(paste("❌ Arquivo de genes não encontrado:", genes_file))
    return(character(0))
  }
  
  genes <- readLines(genes_file, warn = FALSE) %>% 
    trimws() %>% 
    .[. != ""]
  log_message(paste("✅ Carregados", length(genes), "genes de interesse"))
  return(genes)
}

#' Carrega a classificação de amostras (Cancer/Saudavel/HPB/Indeciso/Descartado)
#' @return DataFrame com classificação ou NULL se arquivo não encontrado
load_sample_classification <- function() {
  if (!file.exists(classification_file)) {
    log_message(paste("❌ Arquivo de classificação não encontrado:", classification_file))
    return(NULL)
  }
  
  classification_df <- read.csv(classification_file, stringsAsFactors = FALSE)
  
  # Verificar se as colunas esperadas existem
  required_cols <- c("Experimento", "Corrida", "Tipo_Amostra")
  if (!all(required_cols %in% colnames(classification_df))) {
    log_message(paste("❌ Arquivo de classificação não tem as colunas esperadas. Encontradas:", 
                     paste(colnames(classification_df), collapse = ", ")))
    return(NULL)
  }
  
  log_message(paste("✅ Carregada classificação de", nrow(classification_df), "amostras"))
  
  # Estatísticas da classificação com os novos valores
  stats <- table(classification_df$Tipo_Amostra)
  for (tipo in names(stats)) {
    log_message(paste("   ", tipo, ":", stats[tipo], "amostras"))
  }
  
  # Contar experimentos únicos
  unique_experiments <- unique(classification_df$Experimento)
  log_message(paste("   Total de experimentos na classificação:", length(unique_experiments)))
  
  return(classification_df)
}

# =============================================================================
# SEÇÃO 4: PROCESSAMENTO DE EXPERIMENTOS RNA-SEQ (CORRIGIDO)
# =============================================================================

#' Processa um experimento RNA-seq, aplicando classificação de amostras
#' @param file Caminho do arquivo do experimento
#' @param target_genes Vetor de genes de interesse
#' @param classification_df DataFrame com classificação de amostras
#' @return Lista com dados processados ou NULL em caso de erro
process_rnaseq_experiment <- function(file, target_genes, classification_df) {
  acc <- gsub("\\.(rds|qs)$", "", basename(file))
  
  tryCatch({
    # -------------------------------------------------------------------------
    # SUBSESSÃO 4.1: CARREGAMENTO E VERIFICAÇÃO DO ARQUIVO
    # -------------------------------------------------------------------------
    
    # Carregar arquivo
    if (grepl("\\.rds$", file)) {
      exp_data <- readRDS(file)
    } else {
      exp_data <- qs::qread(file)
    }
    
    # Verificar se é RNA-seq
    is_rnaseq <- FALSE
    rse <- NULL
    
    if (inherits(exp_data, "SimpleList") && "rnaseq" %in% names(exp_data)) {
      rse <- exp_data$rnaseq
      is_rnaseq <- TRUE
    } else if (inherits(exp_data, "RangedSummarizedExperiment")) {
      rse <- exp_data
      is_rnaseq <- TRUE
    } else if (inherits(exp_data, "SummarizedExperiment")) {
      rse <- exp_data
      is_rnaseq <- TRUE
    }
    
    if (!is_rnaseq || is.null(rse)) {
      return(list(data = NULL, has_organism_part = FALSE, reason = "Não é RNA-seq"))
    }
    
    # Verificar se tem counts
    if (!"counts" %in% assayNames(rse)) {
      return(list(data = NULL, has_organism_part = FALSE, reason = "Sem assay 'counts'"))
    }
    
    # -------------------------------------------------------------------------
    # SUBSESSÃO 4.2: EXTRAÇÃO DE DADOS DE EXPRESSÃO
    # -------------------------------------------------------------------------
    
    # Obter matriz de expressão
    expr_matrix <- assay(rse, "counts")
    gene_ids <- rownames(expr_matrix)
    
    # Obter metadados
    sample_metadata <- as.data.frame(colData(rse))
    
    # Verificar se tem organism_part
    if (!"organism_part" %in% colnames(sample_metadata)) {
      return(list(data = NULL, has_organism_part = FALSE, reason = "Sem coluna 'organism_part'"))
    }
    
    log_message(paste("📊 RNA-seq", acc, "-", nrow(expr_matrix), "genes x", ncol(expr_matrix), "amostras"))
    
    # Filtrar genes de interesse
    common_genes <- intersect(gene_ids, target_genes)
    
    if (length(common_genes) == 0) {
      return(list(data = NULL, has_organism_part = TRUE, reason = "Nenhum gene de interesse"))
    }
    
    log_message(paste("🎯 Genes encontrados em", acc, ":", length(common_genes)))
    
    # Criar dataframe de expressão no formato longo
    expr_df <- as.data.frame(expr_matrix[common_genes, , drop = FALSE]) %>%
      tibble::rownames_to_column("GeneID") %>%
      pivot_longer(cols = -GeneID, names_to = "Sample", values_to = "Expression") %>%
      mutate(Experiment = acc, Sample_ID = paste(acc, Sample, sep = "_"))
    
    # -------------------------------------------------------------------------
    # SUBSESSÃO 4.3: CLASSIFICAÇÃO DE AMOSTRAS
    # -------------------------------------------------------------------------
    
    # Extrair metadados de tecido
    sample_meta <- sample_metadata %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Sample") %>%
      select(Sample, Organism_Part = organism_part)
    
    log_message(paste("   🔍 Procurando correspondência para", nrow(sample_meta), "amostras na classificação"))
    
    # Estratégias de correspondência de nomes de amostras
    sample_meta <- sample_meta %>%
      mutate(
        Sample_Clean = gsub("[^a-zA-Z0-9_.-]", "", Sample),      # Remove caracteres especiais
        Sample_Simple = gsub("^[A-Za-z]+\\.\\.", "", Sample_Clean) # Remove prefixos
      )
    
    # Filtrar classificação para este experimento
    classification_subset <- classification_df %>% 
      filter(Experimento == acc)
    
    # Se não houver classificação, usar organism_part padrão
    if (nrow(classification_subset) == 0) {
      log_message(paste("   ⚠️ Nenhuma classificação encontrada para o experimento", acc))
      expr_df <- left_join(expr_df, sample_meta %>% select(Sample, Tissue = Organism_Part), by = "Sample")
      expr_df$Tissue <- as.character(expr_df$Tissue)
      return(list(data = expr_df, has_organism_part = TRUE, reason = "Sem classificação específica"))
    }
    
    # -------------------------------------------------------------------------
    # SUBSESSÃO 4.4: CORRESPONDÊNCIA DE AMOSTRAS
    # -------------------------------------------------------------------------
    
    # Tentar diferentes estratégias de junção
    joined_data <- NULL
    
    # Estratégia 1: junção direta
    temp_meta <- sample_meta %>%
      left_join(classification_subset, by = c("Sample" = "Corrida"))
    
    matched_count <- sum(!is.na(temp_meta$Tipo_Amostra))
    log_message(paste("   📋 Estratégia 1 (junção direta):", matched_count, "amostras correspondentes"))
    
    if (matched_count > 0) {
      joined_data <- temp_meta
    } else {
      # Estratégia 2: usar Sample_Clean
      temp_meta <- sample_meta %>%
        left_join(classification_subset, by = c("Sample_Clean" = "Corrida"))
      
      matched_count <- sum(!is.na(temp_meta$Tipo_Amostra))
      log_message(paste("   📋 Estratégia 2 (Sample_Clean):", matched_count, "amostras correspondentes"))
      
      if (matched_count > 0) {
        joined_data <- temp_meta
      } else {
        # Estratégia 3: usar Sample_Simple
        temp_meta <- sample_meta %>%
          left_join(classification_subset, by = c("Sample_Simple" = "Corrida"))
        
        matched_count <- sum(!is.na(temp_meta$Tipo_Amostra))
        log_message(paste("   📋 Estratégia 3 (Sample_Simple):", matched_count, "amostras correspondentes"))
        
        if (matched_count > 0) {
          joined_data <- temp_meta
        }
      }
    }
    
    # Se nenhuma estratégia funcionou, usar organism_part padrão
    if (is.null(joined_data)) {
      log_message(paste("   ⚠️ Nenhuma correspondência encontrada para amostras de", acc))
      expr_df <- left_join(expr_df, sample_meta %>% select(Sample, Tissue = Organism_Part), by = "Sample")
      expr_df$Tissue <- as.character(expr_df$Tissue)
      return(list(data = expr_df, has_organism_part = TRUE, reason = "Sem correspondência na classificação"))
    }
    
    # -------------------------------------------------------------------------
    # SUBSESSÃO 4.5: APLICAÇÃO DA CLASSIFICAÇÃO (CORRIGIDA)
    # -------------------------------------------------------------------------
    
    # CORREÇÃO: Unificar "prostate" e "prostate gland" nas categorias apropriadas
    joined_data <- joined_data %>%
      mutate(
        # Normalizar os nomes de organism_part (tratar "prostate gland" como "prostate")
        Organism_Part_Norm = case_when(
          tolower(Organism_Part) %in% c("prostate", "prostate gland", "prostate.gland") ~ "prostate",
          TRUE ~ tolower(Organism_Part)
        ),
        # Aplicar classificação baseada no Tipo_Amostra
        Tissue = case_when(
          Tipo_Amostra == "Cancer" & Organism_Part_Norm == "prostate" ~ "Prostate_Cancer",
          Tipo_Amostra == "Saudavel" & Organism_Part_Norm == "prostate" ~ "Prostate_Healthy", 
          Tipo_Amostra == "HPB" & Organism_Part_Norm == "prostate" ~ "Prostate_HPB",
          # Se não for próstata, manter o organism_part original
          TRUE ~ Organism_Part_Norm
        )
      ) %>%
      mutate(Tissue = as.character(Tissue))
    
    # Remover amostras indecisas, descartadas e NAs
    amostras_antes <- nrow(joined_data)
    joined_data <- joined_data %>% 
      filter(!Tipo_Amostra %in% c("Indeciso", "Descartado") & !is.na(Tissue))
    amostras_depois <- nrow(joined_data)
    
    if (amostras_antes > amostras_depois) {
      removidas <- amostras_antes - amostras_depois
      log_message(paste("   🗑️ Removidas", removidas, "amostras (Indeciso/Descartado/NA)"))
    }
    
    # Juntar com dados de expressão
    expr_df <- left_join(expr_df, joined_data %>% select(Sample, Tissue), by = "Sample")
    expr_df$Tissue <- as.character(expr_df$Tissue)
    
    # Remover amostras sem Tissue (NA)
    expr_df <- expr_df %>% filter(!is.na(Tissue))
    
    # Verificar se há dados após remoção
    if (nrow(expr_df) == 0) {
      return(list(data = NULL, has_organism_part = TRUE, reason = "Todas as amostras foram removidas"))
    }
    
    # -------------------------------------------------------------------------
    # SUBSESSÃO 4.6: RELATÓRIO FINAL DO EXPERIMENTO
    # -------------------------------------------------------------------------
    
    # Verificar quantas amostras foram classificadas
    classified_data <- expr_df %>%
      distinct(Sample, Tissue)
    
    if (nrow(classified_data) > 0) {
      tissue_counts <- table(classified_data$Tissue)
      for (tissue_name in names(tissue_counts)) {
        log_message(paste("   🎯", tissue_name, ":", tissue_counts[tissue_name], "amostras"))
      }
    }
    
    log_message(paste("   ✅ Classificação aplicada -", nrow(expr_df), "amostras após filtro"))
    
    return(list(data = expr_df, has_organism_part = TRUE, reason = "Sucesso"))
    
  }, error = function(e) {
    log_message(paste("❌ Erro em", acc, ":", e$message))
    return(list(data = NULL, has_organism_part = FALSE, reason = paste("Erro:", e$message)))
  })
}

# =============================================================================
# SEÇÃO 5: FUNÇÕES PARA SALVAMENTO E PROCESSAMENTO DE DADOS
# =============================================================================

#' Salva lista de tecidos encontrados e estatísticas
#' @param all_expression_data Lista com dados de expressão de todos os experimentos
save_tissues_list <- function(all_expression_data) {
  if (length(all_expression_data) == 0) return()
  
  combined_data <- bind_rows(all_expression_data)
  
  # Garantir que a coluna Tissue é caractere
  combined_data$Tissue <- as.character(combined_data$Tissue)
  
  all_tissues <- unique(combined_data$Tissue)
  all_tissues <- sort(all_tissues[!is.na(all_tissues)])
  
  tissues_file <- file.path(gene_expression_dir, "tissues_found.txt")
  writeLines(all_tissues, tissues_file)
  log_message(paste("📝 Lista de", length(all_tissues), "tecidos salva em", tissues_file))
  
  # Contar ocorrências de cada tecido
  tissue_counts <- combined_data %>%
    distinct(Experiment, Sample, Tissue) %>%
    group_by(Tissue) %>%
    summarise(n = n()) %>%
    arrange(desc(n))
  
  log_message("📊 Distribuição de tecidos:")
  for (i in 1:min(10, nrow(tissue_counts))) {
    log_message(paste("   ", tissue_counts$Tissue[i], ":", tissue_counts$n[i], "amostras"))
  }
}

#' Corrige nomes de tecidos substituindo caracteres problemáticos
#' @param df DataFrame com nomes de colunas a serem corrigidos
#' @return DataFrame com nomes corrigidos
correct_tissue_names <- function(df) {
  colnames(df) <- gsub("\\s+", "_", colnames(df))
  colnames(df) <- gsub("[()]", "", colnames(df))
  colnames(df) <- gsub("\\.", "_", colnames(df))
  colnames(df) <- gsub(",", "_", colnames(df))
  return(df)
}

#' Salva lista de experimentos não analisados
#' @param not_analyzed_list Lista de experimentos não analisados
save_not_analyzed <- function(not_analyzed_list) {
  if (length(not_analyzed_list) > 0) {
    not_analyzed_df <- data.frame(
      Experiment = names(not_analyzed_list),
      Reason = unlist(not_analyzed_list)
    )
    write_csv(not_analyzed_df, not_analyzed_file)
    log_message(paste("📝 Salvo", length(not_analyzed_list), "experimentos não analisados em", not_analyzed_file))
  }
}

# =============================================================================
# SEÇÃO 6: CRIAÇÃO DE ARQUIVOS CSV (CORRIGIDA)
# =============================================================================

#' Cria arquivos CSV no formato wide (0_ e 1_)
#' @param all_expression_data Lista com dados de expressão
create_gene_csv_files_wide <- function(all_expression_data) {
  if (length(all_expression_data) == 0) {
    log_message("❌ Nenhum dado de expressão disponível")
    return()
  }
  
  combined_data <- bind_rows(all_expression_data)
  genes <- unique(combined_data$GeneID)
  log_message(paste("📁 Criando arquivos CSV (formato wide) para", length(genes), "genes"))
  
  for (gene in genes) {
    tryCatch({
      log_message(paste("🔍 Processando gene:", gene))
      
      gene_data <- combined_data %>% filter(GeneID == gene)
      
      # Verificar se há dados para o gene
      if (nrow(gene_data) == 0) {
        log_message(paste("⚠️ Nenhum dado encontrado para o gene", gene))
        next
      }
      
      # -----------------------------------------------------------------------
      # SUBSESSÃO 6.1: ARQUIVO 0_ (DADOS BRUTOS)
      # -----------------------------------------------------------------------
      
      # Filtrar apenas amostras com Tissue definido (remover NAs)
      gene_data_clean <- gene_data %>% filter(!is.na(Tissue))
      
      # Criar formato wide com tecidos como colunas
      wide_data_raw <- gene_data_clean %>%
        select(Sample_ID, Tissue, Expression) %>%
        distinct(Sample_ID, Tissue, .keep_all = TRUE) %>%
        pivot_wider(
          names_from = Tissue,
          values_from = Expression,
          values_fill = list(Expression = NA)
        )
      
      # Corrigir nomes de tecidos
      wide_data_raw <- correct_tissue_names(wide_data_raw)
      
      # Remover coluna "NA_" se existir
      if ("NA_" %in% colnames(wide_data_raw)) {
        wide_data_raw <- wide_data_raw %>% select(-"NA_")
      }
      
      # Remover colunas que são apenas NA
      wide_data_raw <- wide_data_raw[, colSums(is.na(wide_data_raw)) < nrow(wide_data_raw)]
      
      # Reordenar colunas: Sample_ID primeiro, depois tecidos
      sample_col <- wide_data_raw %>% select(Sample_ID)
      tissue_cols <- wide_data_raw %>% select(-Sample_ID)
      wide_data_raw_final <- bind_cols(sample_col, tissue_cols)
      
      # Criar diretório para o gene
      gene_dir <- file.path(gene_expression_dir, gene)
      dir.create(gene_dir, showWarnings = FALSE, recursive = TRUE)
      
      # Salvar arquivo bruto (0_)
      output_file_raw <- file.path(gene_dir, paste0("0_", gene, ".csv"))
      write_csv(wide_data_raw_final, output_file_raw)
      log_message(paste("   ✅ Arquivo 0_ salvo:", output_file_raw))
      
      # -----------------------------------------------------------------------
      # SUBSESSÃO 6.2: ARQUIVO 1_ (MÉDIAS POR EXPERIMENTO)
      # -----------------------------------------------------------------------
      
      # Processar o arquivo 0_ para criar o 1_ com médias por experimento
      raw_data <- read.csv(output_file_raw, stringsAsFactors = FALSE)
      
      # Remover colunas vazias ou todas NA
      raw_data <- raw_data[, colSums(is.na(raw_data)) < nrow(raw_data)]
      
      sample_ids <- raw_data$Sample_ID
      tissue_cols <- colnames(raw_data)[-1]  # Todas as colunas exceto Sample_ID
      
      # Criar uma lista para armazenar as médias por experimento
      experiment_means <- list()
      
      # Processar cada Sample_ID para extrair o experimento
      for (i in 1:length(sample_ids)) {
        sample_id <- sample_ids[i]
        # Extrair código do experimento (parte antes do primeiro _)
        experiment_code <- strsplit(sample_id, "_")[[1]][1]
        
        # Para cada tecido, acumular os valores por experimento
        for (tissue in tissue_cols) {
          value <- raw_data[i, tissue]
          if (!is.na(value)) {
            key <- paste(experiment_code, tissue, sep = "|")
            if (!key %in% names(experiment_means)) {
              experiment_means[[key]] <- list(values = c(), experiment = experiment_code, tissue = tissue)
            }
            experiment_means[[key]]$values <- c(experiment_means[[key]]$values, value)
          }
        }
      }
      
      # Calcular médias e criar estrutura para o novo dataframe
      unique_experiments <- unique(sapply(strsplit(sample_ids, "_"), function(x) x[1]))
      mean_df <- data.frame(Sample_ID = unique_experiments)
      
      # Adicionar colunas para cada tecido
      for (tissue in tissue_cols) {
        mean_df[[tissue]] <- NA
      }
      
      # Preencher com as médias calculadas - ARREDONDANDO PARA INTEIROS
      for (key in names(experiment_means)) {
        exp_data <- experiment_means[[key]]
        experiment <- exp_data$experiment
        tissue <- exp_data$tissue
        mean_val <- round(mean(exp_data$values, na.rm = TRUE))
        
        # Encontrar a linha correspondente ao experimento e atualizar o valor
        row_idx <- which(mean_df$Sample_ID == experiment)
        if (length(row_idx) > 0) {
          mean_df[row_idx, tissue] <- mean_val
        }
      }
      
      # Salvar arquivo de médias (1_)
      output_file_mean <- file.path(gene_dir, paste0("1_", gene, ".csv"))
      write.csv(mean_df, output_file_mean, row.names = FALSE)
      log_message(paste("   ✅ Arquivo 1_ salvo:", output_file_mean))
      
      log_message(paste("✅", gene, "- 0 e 1 criados -", 
                        nrow(wide_data_raw_final), "amostras,", 
                        nrow(mean_df), "experimentos,",
                        ncol(wide_data_raw_final)-1, "tecidos"))
      
    }, error = function(e) {
      log_message(paste("❌ Erro com gene", gene, ":", e$message))
    })
  }
}

#' Cria arquivos 2_ com agrupamentos de tecidos
create_gene_csv_files_2 <- function() {
  # Encontrar todos os diretórios de genes
  gene_dirs <- list.dirs(gene_expression_dir, recursive = FALSE, full.names = TRUE)
  gene_dirs <- gene_dirs[grepl("ENSG", basename(gene_dirs))]
  
  log_message(paste("📁 Criando arquivos 2_ para", length(gene_dirs), "genes"))
  
  # ---------------------------------------------------------------------------
  # SUBSESSÃO 6.3: DEFINIÇÃO DOS AGRUPAMENTOS DE TECIDOS
  # ---------------------------------------------------------------------------
  
  # Manter todos os tecidos individuais, apenas agrupar sinônimos
  tissue_groups <- list(
    adipose = c("adipose", "adipose_tissue"),
    adrenal = c("adrenal", "adrenal_gland"),
    appendix = c("appendix", "vermiform_appendix"),
    bladder = c("bladder", "urinary_bladder"),
    blood = c("blood", "venous_blood", "leukocyte"),
    brain_general = c("brain", "brain_fragment"),
    diencephalon = c("diencephalon", "diencephalon_and_midbrain", "midbrain"),
    forebrain = c("forebrain", "forebrain_and_midbrain", "forebrain_fragment"),
    frontal_lobe = c("frontal_lobe", "prefrontal_cortex"),
    hindbrain = c("hindbrain", "hindbrain_fragment", "hindbrain_without_cerebellum"),
    intestine = c("intestine", "small_intestine", "duodenum"),
    brain_stem = c("medulla_oblongata", "spinal_cord"),
    salivary_gland = c("saliva_secreting_gland", "salivary_gland"),
    skeletal_muscle = c("skeletal_muscle", "skeletal_muscle_tissue"),
    skin = c("skin", "zone_of_skin"),
    thyroid = c("thyroid", "thyroid_gland"),
    umbilical_cord = c("umbilical_cord_blood", "umbilical_vein")
  )
  
  # NOTA: Prostate_Cancer, Prostate_Healthy, Prostate_HPB serão mantidos separadamente
  # Não vamos agrupar prostate_gland pois ele já foi classificado nas categorias acima
  
  # Tecidos para remover completamente
  tissues_to_remove <- c("animal_ovary", "reproductive_organ", "NA_")
  
  for (gene_dir in gene_dirs) {
    tryCatch({
      gene <- basename(gene_dir)
      log_message(paste("🔍 Processando arquivo 2_ para:", gene))
      
      # Caminho para o arquivo 1_
      file_1 <- file.path(gene_dir, paste0("1_", gene, ".csv"))
      
      if (!file.exists(file_1)) {
        log_message(paste("⚠️ Arquivo 1_ não encontrado para", gene))
        next
      }
      
      # Ler o arquivo 1_
      data_1 <- read.csv(file_1, stringsAsFactors = FALSE)
      
      # Criar cópia para o arquivo 2_
      data_2 <- data_1
      
      # Remover tecidos especificados
      cols_to_remove <- which(colnames(data_2) %in% tissues_to_remove)
      if (length(cols_to_remove) > 0) {
        data_2 <- data_2[, -cols_to_remove, drop = FALSE]
      }
      
      # Aplicar agrupamentos de tecidos (apenas sinônimos)
      for (group_name in names(tissue_groups)) {
        group_tissues <- tissue_groups[[group_name]]
        
        # Verificar quais tecidos do grupo estão presentes nos dados
        existing_tissues <- group_tissues[group_tissues %in% colnames(data_2)]
        
        if (length(existing_tissues) > 0) {
          # Calcular a média para cada linha (experimento) - ARREDONDANDO
          group_means <- apply(data_2[, existing_tissues, drop = FALSE], 1, function(row) {
            # Filtrar valores não-NA e > 0
            valid_values <- row[!is.na(row) & row > 0]
            if (length(valid_values) > 0) {
              return(round(mean(valid_values)))  # ARREDONDAMENTO
            } else {
              return(NA)
            }
          })
          
          # Adicionar a nova coluna com o nome do grupo
          data_2[[group_name]] <- group_means
          
          # Remover as colunas originais que foram agrupadas
          data_2 <- data_2[, !colnames(data_2) %in% existing_tissues, drop = FALSE]
        }
      }
      
      # Remover colunas que contêm apenas NA
      cols_before <- ncol(data_2)
      na_columns <- sapply(data_2, function(col) all(is.na(col)))
      # Manter a coluna Sample_ID mesmo se for só NA
      if ("Sample_ID" %in% colnames(data_2)) {
        na_columns["Sample_ID"] <- FALSE
      }
      data_2 <- data_2[, !na_columns, drop = FALSE]
      cols_after <- ncol(data_2)
      
      # Salvar arquivo 2_
      output_file_2 <- file.path(gene_dir, paste0("2_", gene, ".csv"))
      write.csv(data_2, output_file_2, row.names = FALSE)
      
      log_message(paste("   ✅ Arquivo 2_ salvo:", output_file_2))
      log_message(paste("   📊 Colunas:", paste(colnames(data_2), collapse = ", ")))
      
    }, error = function(e) {
      log_message(paste("❌ Erro ao processar arquivo 2_ para", basename(gene_dir), ":", e$message))
    })
  }
  
  log_message(paste("✅ Arquivos 2_ criados com sucesso"))
}

#' Cria arquivos 3_ com as médias globais por gene
create_gene_csv_files_3 <- function() {
  # Encontrar todos os diretórios de genes
  gene_dirs <- list.dirs(gene_expression_dir, recursive = FALSE, full.names = TRUE)
  gene_dirs <- gene_dirs[grepl("ENSG", basename(gene_dirs))]
  
  log_message(paste("📁 Criando arquivos 3_ para", length(gene_dirs), "genes"))
  
  for (gene_dir in gene_dirs) {
    tryCatch({
      gene <- basename(gene_dir)
      log_message(paste("🔍 Processando arquivo 3_ para:", gene))
      
      # Caminho para o arquivo 2_
      file_2 <- file.path(gene_dir, paste0("2_", gene, ".csv"))
      
      if (!file.exists(file_2)) {
        log_message(paste("⚠️ Arquivo 2_ não encontrado para", gene))
        next
      }
      
      # Ler o arquivo 2_
      data_2 <- read.csv(file_2, stringsAsFactors = FALSE)
      
      # Remover coluna Sample_ID
      tissue_data <- data_2[, !colnames(data_2) %in% "Sample_ID", drop = FALSE]
      
      # Calcular a média para cada tecido (coluna), considerando apenas valores >0
      mean_values <- sapply(tissue_data, function(col) {
        # Filtrar valores não-NA e >0
        valid_values <- col[!is.na(col) & col > 0]
        if (length(valid_values) > 0) {
          return(round(mean(valid_values)))  # ARREDONDAMENTO
        } else {
          return(NA)
        }
      })
      
      # Criar data frame com uma linha: Gene e os valores médios
      result_df <- data.frame(Gene = gene, as.list(mean_values), stringsAsFactors = FALSE)
      
      # Corrigir nomes de tecidos
      result_df <- correct_tissue_names(result_df)
      
      # Remover coluna "NA_" se existir
      if ("NA_" %in% colnames(result_df)) {
        result_df <- result_df %>% select(-"NA_")
      }
      
      # Salvar arquivo 3_
      output_file_3 <- file.path(gene_dir, paste0("3_", gene, ".csv"))
      write.csv(result_df, output_file_3, row.names = FALSE)
      log_message(paste("   ✅ Arquivo 3_ salvo:", output_file_3))
      
    }, error = function(e) {
      log_message(paste("❌ Erro ao processar arquivo 3_ para", basename(gene_dir), ":", e$message))
    })
  }
  
  log_message(paste("✅ Arquivos 3_ criados com sucesso"))
}

#' Cria arquivo final com todos os genes
create_all_genes_final_file <- function() {
  log_message("📁 Criando arquivo All_Genes_Final.csv")
  
  # Encontrar todos os diretórios de genes
  gene_dirs <- list.dirs(gene_expression_dir, recursive = FALSE, full.names = TRUE)
  gene_dirs <- gene_dirs[grepl("ENSG", basename(gene_dirs))]
  
  all_genes_data <- list()
  
  for (gene_dir in gene_dirs) {
    tryCatch({
      gene <- basename(gene_dir)
      
      # Caminho para o arquivo 3_
      file_3 <- file.path(gene_dir, paste0("3_", gene, ".csv"))
      
      if (!file.exists(file_3)) {
        log_message(paste("⚠️ Arquivo 3_ não encontrado para", gene))
        next
      }
      
      # Ler o arquivo 3_
      gene_data <- read.csv(file_3, stringsAsFactors = FALSE)
      
      # Remover coluna "NA_" se existir
      if ("NA_" %in% colnames(gene_data)) {
        gene_data <- gene_data %>% select(-"NA_")
      }
      
      # Adicionar à lista
      all_genes_data[[gene]] <- gene_data
      
    }, error = function(e) {
      log_message(paste("❌ Erro ao ler arquivo 3_ para", basename(gene_dir), ":", e$message))
    })
  }
  
  # Combinar todos os dados
  if (length(all_genes_data) > 0) {
    combined_data <- bind_rows(all_genes_data)
    
    # Reordenar colunas: Gene primeiro, depois as demais colunas (tecidos)
    gene_col <- combined_data %>% select(Gene)
    tissue_cols <- combined_data %>% select(-Gene)
    final_data <- bind_cols(gene_col, tissue_cols)
    
    # Salvar arquivo final
    output_file <- file.path(gene_expression_dir, "All_Genes_Final.csv")
    write.csv(final_data, output_file, row.names = FALSE)
    log_message(paste("✅ Arquivo All_Genes_Final.csv salvo com", nrow(final_data), "genes e", ncol(final_data)-1, "tecidos"))
  } else {
    log_message("❌ Nenhum dado de gene encontrado para criar All_Genes_Final.csv")
  }
}

#' Cria arquivos de sumário com estatísticas
#' @param all_expression_data Lista com dados de expressão
create_summary_file <- function(all_expression_data) {
  if (length(all_expression_data) == 0) return()
  combined_data <- bind_rows(all_expression_data)
  
  # Sumário por experimento
  experiment_summary <- combined_data %>%
    group_by(Experiment) %>%
    summarise(
      Amostras = n_distinct(Sample),
      Tecidos = n_distinct(Tissue),
      Genes = n_distinct(GeneID),
      .groups = 'drop'
    )
  
  # Sumário por tecido
  tissue_summary <- combined_data %>%
    group_by(Tissue) %>%
    summarise(
      Experimentos = n_distinct(Experiment),
      Amostras = n_distinct(Sample),
      Genes = n_distinct(GeneID),
      .groups = 'drop'
    )
  
  # Salvar sumários
  summary_file <- file.path(gene_expression_dir, "summary_statistics.csv")
  write_csv(experiment_summary, summary_file)
  
  tissue_file <- file.path(gene_expression_dir, "tissue_statistics.csv")
  write_csv(tissue_summary, tissue_file)
  
  log_message(paste("📊 Arquivos de sumário criados"))
}

# =============================================================================
# SEÇÃO 7: FUNÇÃO PRINCIPAL
# =============================================================================

#' Função principal que orquestra toda a análise
analyze_rnaseq_only <- function() {
  log_message("=== ANÁLISE RNA-seq APENAS (organism_part) COM CLASSIFICAÇÃO CORRIGIDA ===")
  
  # ---------------------------------------------------------------------------
  # SUBSESSÃO 7.1: INICIALIZAÇÃO E CARREGAMENTO DE DADOS
  # ---------------------------------------------------------------------------
  
  # Carregar genes
  target_genes <- load_target_genes()
  if (length(target_genes) == 0) return()
  
  # Carregar classificação de amostras
  classification_df <- load_sample_classification()
  if (is.null(classification_df)) {
    log_message("⚠️  Não foi possível carregar classificação de amostras")
    log_message("⚠️  Continuando análise usando apenas organism_part")
    classification_df <- data.frame()  # DataFrame vazio
  } else {
    log_message("✅ Classificação de amostras carregada com sucesso")
  }
  
  # Listar arquivos de experimentos
  experiment_files <- c(
    list.files(experiment_dir, pattern = "\\.rds$", full.names = TRUE),
    list.files(experiment_dir, pattern = "\\.qs$", full.names = TRUE)
  )
  
  if (length(experiment_files) == 0) {
    log_message("❌ Nenhum experimento encontrado")
    return()
  }
  
  log_message(paste("📂", length(experiment_files), "arquivos encontrados"))
  
  # ---------------------------------------------------------------------------
  # SUBSESSÃO 7.2: PROCESSAMENTO DOS EXPERIMENTOS
  # ---------------------------------------------------------------------------
  
  all_expression_data <- list()
  not_analyzed_list <- list()
  processed_count <- 0
  
  for (file in experiment_files) {
    result <- process_rnaseq_experiment(file, target_genes, classification_df)
    acc <- gsub("\\.(rds|qs)$", "", basename(file))
    
    if (!is.null(result$data) && result$has_organism_part) {
      all_expression_data[[acc]] <- result$data
      processed_count <- processed_count + 1
      log_message(paste("✓", acc, "processado -", result$reason))
    } else {
      not_analyzed_list[[acc]] <- result$reason
      log_message(paste("⏭️ ", acc, "não analisado:", result$reason))
    }
  }
  
  log_message(paste("✅", processed_count, "experimentos RNA-seq com organism_part processados"))
  log_message(paste("⏭️", length(not_analyzed_list), "experimentos não analisados"))
  
  # ---------------------------------------------------------------------------
  # SUBSESSÃO 7.3: SALVAMENTO DE RESULTADOS E RELATÓRIOS
  # ---------------------------------------------------------------------------
  
  # Salvar lista de não analisados
  save_not_analyzed(not_analyzed_list)
  
  # Salvar lista de tecidos encontrados
  save_tissues_list(all_expression_data)
  
  # Criar arquivos CSV no formato wide (0_ e 1_)
  create_gene_csv_files_wide(all_expression_data)
  
  # Criar arquivos 2_ com agrupamentos
  create_gene_csv_files_2()
  
  # Criar arquivos 3_ com médias globais
  create_gene_csv_files_3()
  
  # Criar arquivo final com todos os genes
  create_all_genes_final_file()
  
  # Criar arquivos de sumário
  create_summary_file(all_expression_data)
  
  log_message("=== ANÁLISE RNA-seq COM CLASSIFICAÇÃO CORRIGIDA CONCLUÍDA ===")
}

# =============================================================================
# SEÇÃO 8: EXECUÇÃO PRINCIPAL
# =============================================================================

# Executar análise
log_message("Iniciando análise RNA-seq (apenas organism_part) com classificação corrigida")
analyze_rnaseq_only()
log_message("Análise concluída")

# =============================================================================
# FIM DO SCRIPT
# =============================================================================