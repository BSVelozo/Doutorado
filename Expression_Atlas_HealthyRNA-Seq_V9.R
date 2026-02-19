# Script: Expression_Atlas_RNAseq_Only.R
# Versão que analisa apenas experimentos com "organism_part"

# Carregamento de pacotes
library(readr)
library(tibble)
library(purrr)
library(tidyr)
library(dplyr)
library(SummarizedExperiment)
library(qs)
library(S4Vectors)

# Configurações
experiment_dir <- "/home/bsvelozo/Expression_Atlas/Experiment_Files/Projeto_Saudavel"
gene_expression_dir <- "/home/bsvelozo/Expression_Atlas/Expression_Files_Healthy"
genes_file <- "All_IDs_ENSG.txt"
log_file <- "gene_analysis_log.txt"
not_analyzed_file <- "NotAnalyzed.txt"

# Criar diretório de saída
dir.create(gene_expression_dir, showWarnings = FALSE, recursive = TRUE)

# Função de log
log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  full_msg <- paste0("[", timestamp, "] ", msg)
  cat(full_msg, "\n")
  write(full_msg, file = log_file, append = TRUE)
}

# Carregar genes
load_target_genes <- function() {
  genes <- readLines(genes_file, warn = FALSE) %>% 
    trimws() %>% 
    .[. != ""]
  log_message(paste("✅ Carregados", length(genes), "genes de interesse"))
  return(genes)
}

# Processar APENAS RNA-seq com organism_part
process_rnaseq_experiment <- function(file, target_genes) {
  acc <- gsub("\\.(rds|qs)$", "", basename(file))
  tryCatch({
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
    } else if (is(exp_data, "SummarizedExperiment")) {
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
    
    # Criar dataframe de expressão
    expr_df <- as.data.frame(expr_matrix[common_genes, , drop = FALSE]) %>%
      tibble::rownames_to_column("GeneID") %>%
      pivot_longer(cols = -GeneID, names_to = "Sample", values_to = "Expression") %>%
      mutate(Experiment = acc, Sample_ID = paste(acc, Sample, sep = "_"))
    
    # Adicionar metadados de tecido - APENAS organism_part
    sample_meta <- sample_metadata %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Sample") %>%
      select(Sample, Tissue = organism_part)
    
    expr_df <- left_join(expr_df, sample_meta, by = "Sample")
    log_message(paste("   ✅ Usando coluna de tecido: organism_part"))
    
    return(list(data = expr_df, has_organism_part = TRUE, reason = "Sucesso"))
    
  }, error = function(e) {
    log_message(paste("❌ Erro em", acc, ":", e$message))
    return(list(data = NULL, has_organism_part = FALSE, reason = paste("Erro:", e$message)))
  })
}

# Função para salvar lista de tecidos encontrados
save_tissues_list <- function(all_expression_data) {
  if (length(all_expression_data) == 0) return()
  
  combined_data <- bind_rows(all_expression_data)
  all_tissues <- unique(combined_data$Tissue)
  all_tissues <- sort(all_tissues[!is.na(all_tissues)])
  
  tissues_file <- file.path(gene_expression_dir, "tissues_found.txt")
  writeLines(all_tissues, tissues_file)
  log_message(paste("📝 Lista de", length(all_tissues), "tecidos salva em", tissues_file))
}

# Função para corrigir nomes de tecidos (substituir vírgulas por underscores)
correct_tissue_names <- function(df) {
  # Substituir vírgulas por underscores em todas as colunas (exceto Sample_ID)
  colnames(df) <- gsub(",", "_", colnames(df))
  return(df)
}

# Criar arquivos CSV no formato WIDE por tecido - VERSÃO MODIFICADA COM ARREDONDAMENTO
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
      
      # Criar formato wide com tecidos como colunas - ARQUIVO BRUTO (0_)
      wide_data_raw <- gene_data %>%
        select(Sample_ID, Tissue, Expression) %>%
        distinct(Sample_ID, Tissue, .keep_all = TRUE) %>%
        pivot_wider(
          names_from = Tissue,
          values_from = Expression,
          values_fill = list(Expression = NA)
        )
      
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
      
      # ABORDAGEM COM R BASE - Processar o arquivo 0_ para criar o 1_
      # Ler o arquivo CSV
      raw_data <- read.csv(output_file_raw, stringsAsFactors = FALSE)
      
      # Obter nomes das colunas (tecidos)
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
        mean_val <- round(mean(exp_data$values, na.rm = TRUE))  # ARREDONDAMENTO AQUI
        
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

# Função para criar arquivos 2_ com agrupamentos de tecidos - VERSÃO CORRIGIDA
create_gene_csv_files_2 <- function() {
  # Encontrar todos os diretórios de genes
  gene_dirs <- list.dirs(gene_expression_dir, recursive = FALSE, full.names = TRUE)
  gene_dirs <- gene_dirs[grepl("ENSG", basename(gene_dirs))]
  
  log_message(paste("📁 Criando arquivos 2_ para", length(gene_dirs), "genes"))
  
  # Definir agrupamentos de tecidos - usando pontos como nos arquivos CSV
  tissue_groups <- list(
    adipose = c("adipose", "adipose.tissue"),
    adrenal = c("adrenal", "adrenal.gland"),
    appendix = c("appendix", "vermiform.appendix"),
    bladder = c("bladder", "urinary.bladder"),
    blood = c("blood", "venous.blood", "leukocyte"),
    brain_general = c("brain", "brain.fragment"),
    diencephalon = c("diencephalon", "diencephalon.and.midbrain", "midbrain"),
    forebrain = c("forebrain", "forebrain.and.midbrain", "forebrain.fragment"),
    frontal_lobe = c("frontal.lobe", "prefrontal.cortex"),
    hindbrain = c("hindbrain", "hindbrain.fragment", "hindbrain.without.cerebellum"),
    intestine = c("intestine", "small.intestine", "duodenum"),
    brain_stem = c("medulla.oblongata", "spinal.cord"),
    prostate = c("prostate", "prostate.gland"),
    salivary_gland = c("saliva.secreting.gland", "salivary.gland"),
    skeletal_muscle = c("skeletal.muscle", "skeletal.muscle.tissue"),
    skin = c("skin", "zone.of.skin"),
    thyroid = c("thyroid", "thyroid.gland"),
    umbilical_cord = c("umbilical.cord.blood", "umbilical.vein")
  )
  
  # Tecidos para remover completamente
  tissues_to_remove <- c("animal.ovary", "reproductive.organ")
  
  # Tecidos que não sofrem modificação
  unchanged_tissues <- c(
    "basal.ganglion", "gall.bladder", "bone.marrow", "breast", "cerebellum",
    "cerebral.cortex", "choroid.plexus", "colon", "endometrium", "esophagus",
    "fallopian.tube", "heart", "hippocampus", "kidney", "liver", "lung",
    "lymph.node", "ovary", "pancreas", "pituitary.and.diencephalon", "placenta",
    "pons", "spleen", "stomach", "telencephalon", "temporal.lobe", "testis",
    "thymus", "tonsil"
  )
  
  # Estatísticas globais
  stats_removed_columns <- list()
  stats_created_files <- 0
  
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
      
      # Aplicar agrupamentos de tecidos
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
              return(round(mean(valid_values)))  # ARREDONDAMENTO AQUI
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
      
      # Corrigir nomes de tecidos (substituir vírgulas por underscores)
      data_2 <- correct_tissue_names(data_2)
      
      # Registrar estatísticas
      removed_cols <- names(na_columns)[na_columns]
      if (length(removed_cols) > 0) {
        stats_removed_columns[[gene]] <- removed_cols
      }
      
      # Salvar arquivo 2_
      output_file_2 <- file.path(gene_dir, paste0("2_", gene, ".csv"))
      write.csv(data_2, output_file_2, row.names = FALSE)
      stats_created_files <- stats_created_files + 1
      
      log_message(paste("   ✅ Arquivo 2_ salvo:", output_file_2))
      log_message(paste("   📊 Colunas removidas (só NA):", length(removed_cols)))
      
    }, error = function(e) {
      log_message(paste("❌ Erro ao processar arquivo 2_ para", basename(gene_dir), ":", e$message))
    })
  }
  
  # Salvar estatísticas de colunas removidas - CORRIGIDO
  if (length(stats_removed_columns) > 0) {
    stats_file <- file.path(gene_expression_dir, "removed_columns_stats.txt")
    file_conn <- file(stats_file, "w")
    writeLines("ESTATÍSTICAS DE COLUNAS REMOVIDAS (SÓ NA):", file_conn)
    for (gene in names(stats_removed_columns)) {
      writeLines(paste("\nGene:", gene), file_conn)
      writeLines(paste("Colunas removidas:", paste(stats_removed_columns[[gene]], collapse = ", ")), file_conn)
    }
    close(file_conn)
    log_message(paste("📊 Estatísticas de colunas removidas salvas em:", stats_file))
  }
  
  log_message(paste("✅", stats_created_files, "arquivos 2_ criados com sucesso"))
}

# Função para criar arquivos 3_ com as médias globais por gene - VERSÃO CORRIGIDA
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
          return(round(mean(valid_values)))  # ARREDONDAMENTO AQUI
        } else {
          return(NA)
        }
      })
      
      # Criar data frame com uma linha: Gene e os valores médios
      result_df <- data.frame(Gene = gene, as.list(mean_values), stringsAsFactors = FALSE)
      
      # Corrigir nomes de tecidos (substituir vírgulas por underscores)
      result_df <- correct_tissue_names(result_df)
      
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

# Nova função para criar o arquivo All_Genes_Final.csv
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

# Salvar lista de experimentos não analisados
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

# Função para criar um arquivo de sumário com todos os genes
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

# Função principal
analyze_rnaseq_only <- function() {
  log_message("=== ANÁLISE RNA-seq APENAS (organism_part) ===")
  
  # Carregar genes
  target_genes <- load_target_genes()
  if (length(target_genes) == 0) return()
  
  # Listar arquivos
  experiment_files <- c(
    list.files(experiment_dir, pattern = "\\.rds$", full.names = TRUE),
    list.files(experiment_dir, pattern = "\\.qs$", full.names = TRUE)
  )
  
  if (length(experiment_files) == 0) {
    log_message("❌ Nenhum experimento encontrado")
    return()
  }
  
  log_message(paste("📂", length(experiment_files), "arquivos encontrados"))
  
  # Processar APENAS RNA-seq com organism_part
  all_expression_data <- list()
  not_analyzed_list <- list()
  processed_count <- 0
  
  for (file in experiment_files) {
    result <- process_rnaseq_experiment(file, target_genes)
    acc <- gsub("\\.(rds|qs)$", "", basename(file))
    
    if (!is.null(result$data) && result$has_organism_part && result$reason == "Sucesso") {
      all_expression_data[[acc]] <- result$data
      processed_count <- processed_count + 1
      log_message(paste("✓", acc, "processado com sucesso"))
    } else {
      not_analyzed_list[[acc]] <- result$reason
      log_message(paste("⏭️ ", acc, "não analisado:", result$reason))
    }
  }
  
  log_message(paste("✅", processed_count, "experimentos RNA-seq com organism_part processados"))
  log_message(paste("⏭️", length(not_analyzed_list), "experimentos não analisados"))
  
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
  
  log_message("=== ANÁLISE RNA-seq CONCLUÍDA ===")
}

# Executar
log_message("Iniciando análise RNA-seq (apenas organism_part)")
analyze_rnaseq_only()
log_message("Análise concluída")