# Script: Expression_Atlas_Microarray_Prostata.R
# Versão: Microarray apenas para tecidos da próstata
# Descrição: Analisa experimentos Microarray, classificando amostras de próstata (Cancer/Saudavel/HPB)

library(readr)
library(tibble)
library(purrr)
library(tidyr)
library(dplyr)
library(Biobase)
library(qs)
library(stringr)

# =============================================================================
# CONFIGURAÇÕES
# =============================================================================

experiment_dir <- "/home/bsvelozo/Expression_Atlas/Experiment_Files/Projeto_Cancer_Prostata"
gene_expression_dir <- "/home/bsvelozo/Expression_Atlas/MicroArray_Cancer"
genes_file <- "All_IDs_ENSG.txt"
log_file <- "MIcroArray_Cancer_Log.txt"
mapping_file <- "probe_to_ensg_mapping.csv"
classification_file <- "experimentos_categorizados.csv"

# Criar diretório de saída
dir.create(gene_expression_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# FUNÇÕES AUXILIARES
# =============================================================================

log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  full_msg <- paste0("[", timestamp, "] ", msg)
  cat(full_msg, "\n")
  write(full_msg, file = log_file, append = TRUE)
}

# Carregar genes de interesse
load_target_genes <- function() {
  if (!file.exists(genes_file)) {
    log_message(paste("ERRO: Arquivo de genes não encontrado:", genes_file))
    return(NULL)
  }
  
  genes <- readLines(genes_file, warn = FALSE) %>% 
    trimws() %>% 
    .[. != ""]
  
  log_message(paste("Carregados", length(genes), "genes ENSG de interesse"))
  return(genes)
}

# Carregar mapeamento probe->ENSG
load_probe_mapping <- function() {
  if (!file.exists(mapping_file)) {
    log_message(paste("ERRO: Arquivo de mapeamento não encontrado:", mapping_file))
    return(NULL)
  }
  
  log_message(paste("Carregando mapeamento de:", mapping_file))
  
  mapping <- read_csv(mapping_file, show_col_types = FALSE)
  
  # Verificar colunas essenciais
  required_cols <- c("probe_id", "ensembl_gene_id")
  if (!all(required_cols %in% colnames(mapping))) {
    log_message(paste("ERRO: Colunas necessárias não encontradas. Esperado:", 
                      paste(required_cols, collapse = ", ")))
    return(NULL)
  }
  
  # Limpar e preparar mapeamento
  mapping_clean <- mapping %>%
    select(probe_id, ensembl_gene_id) %>%
    filter(!is.na(probe_id) & !is.na(ensembl_gene_id) &
             probe_id != "" & ensembl_gene_id != "") %>%
    distinct()
  
  log_message(paste("Mapeamento limpo:", nrow(mapping_clean), "mapeamentos"))
  log_message(paste("Probes únicos:", length(unique(mapping_clean$probe_id))))
  log_message(paste("Genes ENSG únicos:", length(unique(mapping_clean$ensembl_gene_id))))
  
  return(mapping_clean)
}

# Carregar classificação de amostras
load_sample_classification <- function() {
  if (!file.exists(classification_file)) {
    log_message(paste("ERRO: Arquivo de classificação não encontrado:", classification_file))
    return(NULL)
  }
  
  classification_df <- read_csv(classification_file, show_col_types = FALSE)
  
  # Verificar colunas
  required_cols <- c("Experimento", "Corrida", "Tipo_Amostra")
  if (!all(required_cols %in% colnames(classification_df))) {
    log_message(paste("ERRO: Colunas necessárias não encontradas. Esperado:", 
                      paste(required_cols, collapse = ", ")))
    return(NULL)
  }
  
  log_message(paste("Classificação carregada:", nrow(classification_df), "amostras"))
  
  # Estatísticas
  stats <- table(classification_df$Tipo_Amostra)
  for (tipo in names(stats)) {
    log_message(paste("  ", tipo, ":", stats[tipo], "amostras"))
  }
  
  log_message(paste("  Experimentos únicos:", length(unique(classification_df$Experimento))))
  
  return(classification_df)
}

# =============================================================================
# PROCESSAMENTO DE EXPERIMENTOS - APENAS PRÓSTATA
# =============================================================================

process_microarray_prostata <- function(file_path, target_genes, probe_mapping, classification_df) {
  acc <- gsub("\\.qs$", "", basename(file_path))
  
  tryCatch({
    log_message(paste("\n=== Processando:", acc, "==="))
    
    # 1. Carregar arquivo .qs
    exp_data <- qs::qread(file_path)
    
    # 2. Extrair ExpressionSet
    eset <- NULL
    if (is(exp_data, "ExpressionSet")) {
      eset <- exp_data
    } else if (inherits(exp_data, "SimpleList")) {
      for (i in seq_along(exp_data)) {
        if (is(exp_data[[i]], "ExpressionSet")) {
          eset <- exp_data[[i]]
          break
        }
      }
    }
    
    if (is.null(eset)) {
      return(list(success = FALSE, data = NULL, reason = "Não é ExpressionSet"))
    }
    
    # 3. Obter dados de expressão
    expr_matrix <- exprs(eset)
    probe_ids <- rownames(expr_matrix)
    pheno_data <- pData(eset)
    
    log_message(paste("  Dimensões:", nrow(expr_matrix), "probes x", ncol(expr_matrix), "amostras"))
    
    # 4. Mapear probes para genes ENSG
    mapping_result <- probe_mapping %>%
      filter(probe_id %in% probe_ids)
    
    if (nrow(mapping_result) == 0) {
      # Tentar matching flexível
      probe_ids_clean <- gsub("_at$|_s_at$|_x_at$|_a_at$", "", probe_ids)
      probe_mapping_clean <- probe_mapping %>%
        mutate(probe_clean = gsub("_at$|_s_at$|_x_at$|_a_at$", "", probe_id))
      
      mapping_result <- data.frame(
        probe_id_experiment = probe_ids,
        probe_clean = probe_ids_clean,
        stringsAsFactors = FALSE
      ) %>%
        inner_join(probe_mapping_clean, by = "probe_clean") %>%
        select(probe_id = probe_id_experiment, ensembl_gene_id) %>%
        distinct()
    }
    
    if (nrow(mapping_result) == 0) {
      return(list(success = FALSE, data = NULL, reason = "Nenhum probe mapeado"))
    }
    
    # 5. Filtrar para genes de interesse
    mapping_target <- mapping_result %>%
      filter(ensembl_gene_id %in% target_genes)
    
    if (nrow(mapping_target) == 0) {
      return(list(success = FALSE, data = NULL, reason = "Nenhum gene de interesse"))
    }
    
    # 6. Selecionar melhor probe por gene (maior variância)
    valid_probes <- mapping_target$probe_id[mapping_target$probe_id %in% rownames(expr_matrix)]
    
    if (length(valid_probes) == 0) {
      return(list(success = FALSE, data = NULL, reason = "Probes não encontrados na matriz"))
    }
    
    expr_subset <- expr_matrix[valid_probes, , drop = FALSE]
    variances <- apply(expr_subset, 1, var, na.rm = TRUE)
    
    mapping_target$variance <- variances[mapping_target$probe_id]
    mapping_target <- mapping_target %>% filter(!is.na(variance))
    
    best_probes <- mapping_target %>%
      group_by(ensembl_gene_id) %>%
      arrange(desc(variance)) %>%
      slice(1) %>%
      ungroup()
    
    log_message(paste("  Genes de interesse mapeados:", nrow(best_probes)))
    
    # 7. Extrair dados de expressão
    expr_best <- expr_matrix[best_probes$probe_id, , drop = FALSE]
    
    # 8. Criar dataframe com expressão
    expr_df <- as.data.frame(expr_best)
    expr_df$probe_id <- rownames(expr_df)
    
    expr_long <- expr_df %>%
      pivot_longer(
        cols = -probe_id,
        names_to = "Sample",
        values_to = "Expression"
      ) %>%
      left_join(best_probes %>% select(probe_id, ensembl_gene_id), by = "probe_id") %>%
      rename(GeneID = ensembl_gene_id) %>%
      select(-probe_id) %>%
      mutate(
        Experiment = acc,
        Sample_ID = paste(acc, Sample, sep = "_")
      )
    
    # 9. APLICAR CLASSIFICAÇÃO DE PRÓSTATA - VERSÃO SIMPLIFICADA
    log_message("  Aplicando classificação de próstata...")
    
    # Filtrar classificação para este experimento
    classification_subset <- classification_df %>%
      filter(Experimento == acc)
    
    if (nrow(classification_subset) > 0) {
      log_message(paste("  Encontradas", nrow(classification_subset), "classificações para", acc))
      
      # Limpar nomes para matching
      classification_clean <- classification_subset %>%
        mutate(
          Corrida_Clean = gsub("[^a-zA-Z0-9_.-]", "", Corrida),
          Corrida_Simple = gsub("^[A-Za-z]+\\.\\.", "", Corrida_Clean)
        )
      
      # Preparar dados de amostras
      samples_clean <- data.frame(
        Sample = colnames(expr_matrix),
        Sample_Clean = gsub("[^a-zA-Z0-9_.-]", "", colnames(expr_matrix)),
        Sample_Simple = gsub("^[A-Za-z]+\\.\\.", "", gsub("[^a-zA-Z0-9_.-]", "", colnames(expr_matrix))),
        stringsAsFactors = FALSE
      )
      
      # Tentar diferentes estratégias de matching
      matched_samples <- NULL
      
      # Estratégia 1: junção direta
      temp_match <- samples_clean %>%
        left_join(classification_clean, by = c("Sample" = "Corrida"))
      
      if (sum(!is.na(temp_match$Tipo_Amostra)) > 0) {
        matched_samples <- temp_match
        log_message("    Matching: estratégia direta")
      } else {
        # Estratégia 2: Sample_Clean
        temp_match <- samples_clean %>%
          left_join(classification_clean, by = c("Sample_Clean" = "Corrida"))
        
        if (sum(!is.na(temp_match$Tipo_Amostra)) > 0) {
          matched_samples <- temp_match
          log_message("    Matching: Sample_Clean")
        } else {
          # Estratégia 3: Sample_Simple
          temp_match <- samples_clean %>%
            left_join(classification_clean, by = c("Sample_Simple" = "Corrida"))
          
          if (sum(!is.na(temp_match$Tipo_Amostra)) > 0) {
            matched_samples <- temp_match
            log_message("    Matching: Sample_Simple")
          }
        }
      }
      
      if (!is.null(matched_samples)) {
        # Juntar com dados de expressão
        expr_long <- expr_long %>%
          left_join(matched_samples %>% select(Sample, Tipo_Amostra), by = "Sample")
        
        # Aplicar classificação de tecido
        expr_long <- expr_long %>%
          mutate(
            Tissue = case_when(
              Tipo_Amostra == "Cancer" ~ "Prostate_Cancer",
              Tipo_Amostra == "Saudavel" ~ "Prostate_Healthy",
              Tipo_Amostra == "HPB" ~ "Prostate_HPB",
              TRUE ~ NA_character_
            )
          ) %>%
          filter(!is.na(Tissue))  # Manter apenas amostras classificadas
        
        # Remover coluna auxiliar
        expr_long <- expr_long %>% select(-Tipo_Amostra)
      } else {
        log_message(paste("  Nenhuma amostra classificada para", acc))
        return(list(success = FALSE, data = NULL, reason = "Nenhuma amostra classificada"))
      }
    } else {
      log_message(paste("  Nenhuma classificação encontrada para", acc))
      return(list(success = FALSE, data = NULL, reason = "Sem classificação"))
    }
    
    # 10. FILTRAR - APENAS AMOSTRAS DE PRÓSTATA CLASSIFICADAS
    # (Já feito acima, mas garantindo)
    expr_long <- expr_long %>%
      filter(!is.na(Tissue) & Tissue %in% c("Prostate_Cancer", "Prostate_Healthy", "Prostate_HPB"))
    
    if (nrow(expr_long) == 0) {
      return(list(success = FALSE, data = NULL, reason = "Nenhuma amostra de próstata após filtro"))
    }
    
    log_message(paste("  ✓ Finalizado:", nrow(expr_long), "linhas,",
                      length(unique(expr_long$GeneID)), "genes,",
                      length(unique(expr_long$Tissue)), "tipos de próstata"))
    
    # Contar tipos de próstata
    tissue_counts <- table(expr_long$Tissue)
    for (tissue_name in names(tissue_counts)) {
      log_message(paste("    ", tissue_name, ":", tissue_counts[tissue_name], "amostras"))
    }
    
    return(list(success = TRUE, data = expr_long))
    
  }, error = function(e) {
    log_message(paste("  ✗ ERRO:", e$message))
    return(list(success = FALSE, data = NULL, reason = paste("Erro:", e$message)))
  })
}

# =============================================================================
# FUNÇÕES PARA CRIAR ARQUIVOS - APENAS PRÓSTATA
# =============================================================================

# Criar arquivo 0_ (dados brutos) - apenas próstata
create_file_0_prostata <- function(gene_data, gene_dir, gene) {
  # Verificar se há dados
  if (nrow(gene_data) == 0) {
    wide_data <- data.frame(Sample_ID = character(), stringsAsFactors = FALSE)
  } else {
    # Garantir que temos apenas os 3 tipos de próstata
    gene_data <- gene_data %>%
      filter(Tissue %in% c("Prostate_Cancer", "Prostate_Healthy", "Prostate_HPB"))
    
    # Criar formato wide
    wide_data <- gene_data %>%
      select(Sample_ID, Tissue, Expression) %>%
      distinct(Sample_ID, Tissue, .keep_all = TRUE) %>%
      pivot_wider(
        names_from = Tissue,
        values_from = Expression,
        values_fill = NA
      )
    
    # Garantir que temos todas as 3 colunas (mesmo que vazias)
    prostate_cols <- c("Prostate_Cancer", "Prostate_Healthy", "Prostate_HPB")
    missing_cols <- setdiff(prostate_cols, colnames(wide_data))
    
    for (col in missing_cols) {
      wide_data[[col]] <- NA_real_
    }
    
    # Ordenar colunas: Sample_ID primeiro, depois as 3 colunas de próstata
    wide_data <- wide_data %>%
      select(Sample_ID, Prostate_Cancer, Prostate_Healthy, Prostate_HPB)
  }
  
  # Salvar
  output_file <- file.path(gene_dir, paste0("0_", gene, ".csv"))
  write_csv(wide_data, output_file, na = "")
  
  return(wide_data)
}

# Criar arquivo 1_ (média por experimento) - apenas próstata
create_file_1_prostata <- function(wide_0_data, gene_dir, gene) {
  if (nrow(wide_0_data) == 0 || ncol(wide_0_data) <= 1) {
    mean_by_experiment <- data.frame(Sample_ID = character(), stringsAsFactors = FALSE)
  } else {
    # Extrair código do experimento
    experiment_codes <- sapply(strsplit(wide_0_data$Sample_ID, "_"), function(x) x[1])
    
    # Adicionar coluna de experimento
    data_with_exp <- wide_0_data
    data_with_exp$Experiment <- experiment_codes
    
    # Garantir que temos as 3 colunas de próstata
    prostate_cols <- c("Prostate_Cancer", "Prostate_Healthy", "Prostate_HPB")
    for (col in prostate_cols) {
      if (!col %in% colnames(data_with_exp)) {
        data_with_exp[[col]] <- NA_real_
      }
    }
    
    # Calcular média por experimento - APENAS PARA PRÓSTATA
    mean_by_experiment <- data_with_exp %>%
      select(Experiment, Prostate_Cancer, Prostate_Healthy, Prostate_HPB) %>%
      group_by(Experiment) %>%
      summarise(
        across(all_of(prostate_cols), 
               ~ {
                 vals <- .[!is.na(.)]
                 if (length(vals) > 0) {
                   return(round(mean(vals)))
                 } else {
                   return(NA_real_)
                 }
               }, 
               .names = "{.col}")
      )
  }
  
  # Salvar
  output_file <- file.path(gene_dir, paste0("1_", gene, ".csv"))
  write_csv(mean_by_experiment, output_file, na = "")
  
  return(mean_by_experiment)
}

# Criar arquivo 2_ (média global) - apenas próstata
create_file_2_prostata <- function(mean_by_experiment, gene) {
  if (nrow(mean_by_experiment) == 0) {
    result_df <- data.frame(
      Gene = gene,
      Prostate_Cancer = NA_real_,
      Prostate_Healthy = NA_real_,
      Prostate_HPB = NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    # Calcular média para cada tipo de próstata
    prostate_cols <- c("Prostate_Cancer", "Prostate_Healthy", "Prostate_HPB")
    
    mean_tissues <- sapply(prostate_cols, function(tissue) {
      col_values <- mean_by_experiment[[tissue]]
      valid_values <- col_values[!is.na(col_values)]
      
      if (length(valid_values) > 0) {
        return(round(mean(valid_values)))
      } else {
        return(NA_real_)
      }
    }, USE.NAMES = TRUE)
    
    result_df <- data.frame(
      Gene = gene,
      Prostate_Cancer = mean_tissues["Prostate_Cancer"],
      Prostate_Healthy = mean_tissues["Prostate_Healthy"],
      Prostate_HPB = mean_tissues["Prostate_HPB"],
      stringsAsFactors = FALSE
    )
  }
  
  return(result_df)
}

# Criar arquivos por gene - apenas próstata
create_gene_files_prostata <- function(all_data) {
  genes <- unique(all_data$GeneID)
  log_message(paste("\nCriando arquivos para", length(genes), "genes (apenas próstata)"))
  
  all_gene_2_files <- list()
  
  for (gene in genes) {
    tryCatch({
      # Criar diretório para o gene
      gene_dir <- file.path(gene_expression_dir, gene)
      dir.create(gene_dir, showWarnings = FALSE, recursive = TRUE)
      
      # Filtrar dados do gene
      gene_data <- all_data %>% filter(GeneID == gene)
      
      # Criar arquivo 0_ (dados brutos)
      wide_0 <- create_file_0_prostata(gene_data, gene_dir, gene)
      
      # Criar arquivo 1_ (média por experimento)
      mean_by_exp <- create_file_1_prostata(wide_0, gene_dir, gene)
      
      # Criar arquivo 2_ (média global)
      mean_global <- create_file_2_prostata(mean_by_exp, gene)
      
      # Salvar arquivo 2_
      output_file_2 <- file.path(gene_dir, paste0("2_", gene, ".csv"))
      write_csv(mean_global, output_file_2, na = "")
      
      # Armazenar para criar All_Genes.csv depois
      all_gene_2_files[[gene]] <- mean_global
      
      # Log
      non_na_counts <- sum(!is.na(mean_global[-1]))  # Contar valores não-NA (excluindo Gene)
      log_message(paste("  ✓", gene, "-", non_na_counts, "tipo(s) de próstata com expressão"))
      
    }, error = function(e) {
      log_message(paste("  ✗ Erro com gene", gene, ":", e$message))
    })
  }
  
  return(all_gene_2_files)
}

# Criar arquivo All_Genes.csv - apenas próstata
create_all_genes_file_prostata <- function(all_gene_2_files) {
  if (length(all_gene_2_files) == 0) {
    log_message("Nenhum arquivo 2_ para consolidar")
    return(NULL)
  }
  
  # Combinar todos os dataframes
  all_genes_df <- bind_rows(all_gene_2_files) %>% arrange(Gene)
  
  # Garantir que temos as 3 colunas
  prostate_cols <- c("Prostate_Cancer", "Prostate_Healthy", "Prostate_HPB")
  for (col in prostate_cols) {
    if (!col %in% colnames(all_genes_df)) {
      all_genes_df[[col]] <- NA_real_
    }
  }
  
  # Reordenar colunas
  all_genes_df <- all_genes_df %>%
    select(Gene, Prostate_Cancer, Prostate_Healthy, Prostate_HPB)
  
  # Salvar arquivo
  output_file <- file.path(gene_expression_dir, "All_Genes.csv")
  write_csv(all_genes_df, output_file, na = "")
  
  # Estatísticas
  cancer_count <- sum(!is.na(all_genes_df$Prostate_Cancer))
  healthy_count <- sum(!is.na(all_genes_df$Prostate_Healthy))
  hpb_count <- sum(!is.na(all_genes_df$Prostate_HPB))
  
  log_message(paste("\nArquivo All_Genes.csv criado com", nrow(all_genes_df), "genes"))
  log_message(paste("  Prostate_Cancer:", cancer_count, "genes com expressão"))
  log_message(paste("  Prostate_Healthy:", healthy_count, "genes com expressão"))
  log_message(paste("  Prostate_HPB:", hpb_count, "genes com expressão"))
  
  return(all_genes_df)
}

# =============================================================================
# FUNÇÃO PRINCIPAL - APENAS PRÓSTATA
# =============================================================================

main_prostata <- function() {
  log_message("=========================================")
  log_message("MICROARRAY - APENAS PRÓSTATA")
  log_message("=========================================")
  
  # 1. Carregar genes de interesse
  log_message("\n1. CARREGANDO GENES DE INTERESSE")
  target_genes <- load_target_genes()
  if (is.null(target_genes) || length(target_genes) == 0) {
    log_message("ERRO: Nenhum gene carregado")
    return()
  }
  
  # 2. Carregar mapeamento
  log_message("\n2. CARREGANDO MAPEAMENTO")
  probe_mapping <- load_probe_mapping()
  if (is.null(probe_mapping)) {
    log_message("ERRO: Mapeamento não carregado")
    return()
  }
  
  # 3. Carregar classificação
  log_message("\n3. CARREGANDO CLASSIFICAÇÃO")
  classification_df <- load_sample_classification()
  if (is.null(classification_df)) {
    log_message("ERRO: Classificação não carregada - necessário para análise de próstata")
    return()
  }
  
  # 4. Listar arquivos
  log_message("\n4. LOCALIZANDO ARQUIVOS")
  files <- list.files(experiment_dir, pattern = "\\.qs$", full.names = TRUE)
  log_message(paste("Encontrados", length(files), "arquivos"))
  
  # 5. Processar arquivos - APENAS PRÓSTATA
  log_message("\n5. PROCESSANDO ARQUIVOS (apenas próstata)")
  all_results <- list()
  successful <- 0
  failed <- 0
  failures <- list()
  
  for (i in seq_along(files)) {
    file <- files[i]
    log_message(paste("\n[", i, "/", length(files), "]", basename(file)))
    
    result <- process_microarray_prostata(file, target_genes, probe_mapping, classification_df)
    
    if (result$success && !is.null(result$data) && nrow(result$data) > 0) {
      acc <- gsub("\\.qs$", "", basename(file))
      all_results[[acc]] <- result$data
      successful <- successful + 1
    } else {
      failed <- failed + 1
      failures[[basename(file)]] <- result$reason
      log_message(paste("  ✗ FALHA:", result$reason))
    }
  }
  
  # 6. Relatório
  log_message("\n=========================================")
  log_message("RELATÓRIO FINAL - APENAS PRÓSTATA")
  log_message("=========================================")
  log_message(paste("Sucesso:", successful, "/", length(files)))
  log_message(paste("Falha:", failed, "/", length(files)))
  
  if (failed > 0) {
    log_message("\nArquivos com falha:")
    for (f in names(failures)) {
      log_message(paste("  -", f, ":", failures[[f]]))
    }
  }
  
  # 7. Salvar resultados
  if (successful > 0) {
    log_message("\n6. CRIANDO ARQUIVOS POR GENE")
    
    # Combinar todos os dados
    all_data <- bind_rows(all_results)
    
    # Criar arquivos por gene
    all_gene_2_files <- create_gene_files_prostata(all_data)
    
    # Criar arquivo consolidado
    log_message("\n7. CRIANDO ARQUIVO CONSOLIDADO")
    all_genes_df <- create_all_genes_file_prostata(all_gene_2_files)
    
    # Salvar dados brutos
    log_message("\n8. SALVANDO DADOS BRUTOS")
    write_csv(all_data, file.path(gene_expression_dir, "all_expression_data_prostata.csv"))
    log_message("Dados brutos salvos em: all_expression_data_prostata.csv")
    
    # Estatísticas finais
    if (!is.null(all_genes_df)) {
      log_message("\n9. ESTATÍSTICAS FINAIS")
      
      # Contar genes por tipo de próstata
      genes_with_data <- apply(all_genes_df[, -1], 2, function(col) sum(!is.na(col)))
      
      for (tissue in names(genes_with_data)) {
        log_message(paste("  ", tissue, ":", genes_with_data[tissue], "genes com expressão"))
      }
      
      # Genes com expressão em múltiplos tipos
      genes_multiple <- sum(rowSums(!is.na(all_genes_df[, -1])) > 1)
      log_message(paste("  Genes com expressão em múltiplos tipos:", genes_multiple))
    }
  }
  
  log_message("\n=========================================")
  log_message("ANÁLISE CONCLUÍDA")
  log_message("=========================================")
}

# =============================================================================
# EXECUÇÃO
# =============================================================================

log_message("Iniciando análise de microarray (apenas próstata)...")
main_prostata()