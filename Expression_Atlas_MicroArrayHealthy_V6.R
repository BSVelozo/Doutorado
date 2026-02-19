# Script: Expression_Atlas_Microarray_Organizado_v2.R
# Versão organizada com NA para tecidos sem expressão

library(readr)
library(tibble)
library(purrr)
library(tidyr)
library(dplyr)
library(Biobase)
library(qs)
library(stringr)

# Configurações
experiment_dir <- "/home/bsvelozo/Expression_Atlas/Experiment_Files/MicroArray_TESTE"
gene_expression_dir <- "/home/bsvelozo/Expression_Atlas/TESTE_MicroArray"
genes_file <- "Some_IDs_ENSG.txt"
log_file <- "gene_analysis_organizado_v2.txt"
mapping_file <- "probe_to_ensg_mapping.csv"

# Criar diretório de saída
dir.create(gene_expression_dir, showWarnings = FALSE, recursive = TRUE)

# Função de log
log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  full_msg <- paste0("[", timestamp, "] ", msg)
  cat(full_msg, "\n")
  write(full_msg, file = log_file, append = TRUE)
}

# Carregar genes de interesse (ENSG)
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

# Função para processar um arquivo .qs
process_expression_set <- function(file_path, target_genes, probe_mapping) {
  acc <- gsub("\\.qs$", "", basename(file_path))
  
  tryCatch({
    log_message(paste("\n=== Processando:", acc, "==="))
    
    # Carregar arquivo .qs
    exp_data <- qs::qread(file_path)
    
    # Extrair ExpressionSet
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
      return(list(success = FALSE, data = NULL, reason = "Não foi possível extrair ExpressionSet"))
    }
    
    # Obter dados de expressão
    expr_matrix <- exprs(eset)
    probe_ids <- rownames(expr_matrix)
    
    # Obter metadados
    pheno_data <- pData(eset)
    
    # Mapear probes para genes ENSG
    # Estratégia 1: Matching exato
    mapping_result <- probe_mapping %>%
      filter(probe_id %in% probe_ids)
    
    if (nrow(mapping_result) == 0) {
      # Estratégia 2: Matching flexível
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
      return(list(success = FALSE, data = NULL, reason = "Nenhum probe mapeado para ENSG"))
    }
    
    # Filtrar para genes de interesse
    mapping_target <- mapping_result %>%
      filter(ensembl_gene_id %in% target_genes)
    
    if (nrow(mapping_target) == 0) {
      return(list(success = FALSE, data = NULL, reason = "Nenhum gene de interesse após mapeamento"))
    }
    
    # Selecionar melhor probe por gene (maior variância)
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
    
    # Extrair dados de expressão
    expr_best <- expr_matrix[best_probes$probe_id, , drop = FALSE]
    
    # Criar dataframe final
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
        SampleID = paste(acc, Sample, sep = "_")
      )
    
    # Adicionar metadados de tecido
    tissue_cols <- c("organism_part", "cell_type", "tissue", "cell_line", "source_name")
    tissue_col <- NULL
    
    for (col in tissue_cols) {
      if (col %in% colnames(pheno_data)) {
        tissue_col <- col
        break
      }
    }
    
    if (!is.null(tissue_col)) {
      sample_meta <- pheno_data %>%
        as.data.frame() %>%
        rownames_to_column("Sample") %>%
        select(Sample, Tissue = !!tissue_col)
      
      expr_long <- expr_long %>%
        left_join(sample_meta, by = "Sample")
    } else {
      expr_long$Tissue <- "Desconhecido"
    }
    
    # Limpar nome do tecido (substituir caracteres problemáticos)
    expr_long <- expr_long %>%
      mutate(Tissue = gsub("\\s+", "_", Tissue),
             Tissue = gsub("[(),]", "", Tissue),
             Tissue = gsub("\\.", "_", Tissue))
    
    # REMOVER: Filtrar amostras sem tecido identificado
    expr_long <- expr_long %>% 
      filter(!is.na(Tissue) & Tissue != "Desconhecido" & Tissue != "")
    
    log_message(paste("  ✓ Sucesso:", length(unique(expr_long$GeneID)), "genes,",
                     length(unique(expr_long$SampleID)), "amostras,",
                     length(unique(expr_long$Tissue)), "tecidos"))
    
    return(list(success = TRUE, data = expr_long))
    
  }, error = function(e) {
    log_message(paste("  ✗ ERRO:", e$message))
    return(list(success = FALSE, data = NULL, reason = paste("Erro:", e$message)))
  })
}

# Função para criar arquivo 0_NomeDoGene.csv (dados brutos)
create_file_0 <- function(gene_data, gene_dir, gene) {
  # Verificar se há dados
  if (nrow(gene_data) == 0) {
    # Criar um dataframe vazio
    wide_data <- data.frame(SampleID = character(), stringsAsFactors = FALSE)
  } else {
    # Criar formato wide com tecidos como colunas
    wide_data <- gene_data %>%
      select(SampleID, Tissue, Expression) %>%
      distinct(SampleID, Tissue, .keep_all = TRUE) %>%
      pivot_wider(
        names_from = Tissue,
        values_from = Expression,
        values_fill = NA  # Garantir que valores faltantes sejam NA
      )
    
    # Manter apenas colunas que são tecidos (excluir SampleID)
    tissue_cols <- setdiff(colnames(wide_data), "SampleID")
    
    # Se não houver tecidos, criar dataframe vazio
    if (length(tissue_cols) == 0) {
      wide_data <- data.frame(SampleID = unique(gene_data$SampleID), stringsAsFactors = FALSE)
    }
  }
  
  # Salvar arquivo 0_
  output_file <- file.path(gene_dir, paste0("0_", gene, ".csv"))
  write_csv(wide_data, output_file, na = "")
  
  return(wide_data)
}

# Função para criar arquivo 1_NomeDoGene.csv (média entre corridas)
create_file_1 <- function(wide_0_data, gene_dir, gene) {
  # Verificar se há dados
  if (nrow(wide_0_data) == 0 || ncol(wide_0_data) <= 1) {
    # Criar dataframe vazio
    mean_by_experiment <- data.frame(Experiment = character(), stringsAsFactors = FALSE)
  } else {
    # Extrair código do experimento do SampleID (parte antes do primeiro _)
    experiment_codes <- sapply(strsplit(wide_0_data$SampleID, "_"), function(x) x[1])
    
    # Adicionar coluna de experimento
    data_with_exp <- wide_0_data
    data_with_exp$Experiment <- experiment_codes
    
    # Identificar colunas que são tecidos (não são SampleID ou Experiment)
    tissue_cols <- setdiff(colnames(data_with_exp), c("SampleID", "Experiment"))
    
    # Calcular média por experimento para cada tecido - ARREDONDANDO
    if (length(tissue_cols) > 0) {
      mean_by_experiment <- data_with_exp %>%
        select(-SampleID) %>%
        group_by(Experiment) %>%
        summarise(across(all_of(tissue_cols), 
                        ~ {
                          # Remover NAs e calcular média apenas se houver pelo menos 1 valor
                          vals <- .[!is.na(.)]
                          if (length(vals) > 0) {
                            return(round(mean(vals)))
                          } else {
                            return(NA_real_)
                          }
                        }, 
                        .names = "{.col}"))
    } else {
      mean_by_experiment <- data_with_exp %>%
        select(-SampleID) %>%
        distinct(Experiment) %>%
        as.data.frame()
    }
  }
  
  # Salvar arquivo 1_
  output_file <- file.path(gene_dir, paste0("1_", gene, ".csv"))
  write_csv(mean_by_experiment, output_file, na = "")
  
  return(mean_by_experiment)
}

# Função para criar arquivo 2_NomeDoGene.csv (média entre experimentos)
create_file_2 <- function(mean_by_experiment, gene) {
  # Verificar se há dados além da coluna Experiment
  if (nrow(mean_by_experiment) == 0 || ncol(mean_by_experiment) <= 1) {
    # Criar dataframe vazio com apenas o Gene
    result_df <- data.frame(Gene = gene, stringsAsFactors = FALSE)
    return(result_df)
  }
  
  # Remover coluna Experiment
  tissue_cols <- setdiff(colnames(mean_by_experiment), "Experiment")
  
  if (length(tissue_cols) == 0) {
    # Criar dataframe com apenas o Gene
    result_df <- data.frame(Gene = gene, stringsAsFactors = FALSE)
    return(result_df)
  }
  
  # Calcular média para cada tecido, considerando apenas valores válidos (>0 e não NA)
  mean_tissues <- sapply(tissue_cols, function(tissue) {
    col_values <- mean_by_experiment[[tissue]]
    
    # Filtrar valores não-NA
    valid_values <- col_values[!is.na(col_values)]
    
    if (length(valid_values) > 0) {
      # Calcular média e arredondar
      return(round(mean(valid_values)))
    } else {
      # Se todos os valores são NA, retornar NA
      return(NA_real_)
    }
  }, USE.NAMES = TRUE)
  
  # Criar dataframe com uma linha
  result_df <- data.frame(Gene = gene, as.list(mean_tissues), stringsAsFactors = FALSE)
  
  return(result_df)
}

# Função para criar arquivos por gene (0_, 1_, 2_)
create_gene_files <- function(all_data) {
  genes <- unique(all_data$GeneID)
  log_message(paste("\nCriando arquivos para", length(genes), "genes"))
  
  all_gene_2_files <- list()
  
  for (gene in genes) {
    tryCatch({
      # Criar diretório para o gene
      gene_dir <- file.path(gene_expression_dir, gene)
      dir.create(gene_dir, showWarnings = FALSE, recursive = TRUE)
      
      # Filtrar dados do gene
      gene_data <- all_data %>% filter(GeneID == gene)
      
      # Criar arquivo 0_ (dados brutos)
      wide_0 <- create_file_0(gene_data, gene_dir, gene)
      
      # Criar arquivo 1_ (média por experimento)
      mean_by_exp <- create_file_1(wide_0, gene_dir, gene)
      
      # Criar arquivo 2_ (média global)
      mean_global <- create_file_2(mean_by_exp, gene)
      
      # Salvar arquivo 2_
      output_file_2 <- file.path(gene_dir, paste0("2_", gene, ".csv"))
      write_csv(mean_global, output_file_2, na = "")
      
      # Armazenar para criar All_Genes.csv depois
      all_gene_2_files[[gene]] <- mean_global
      
      # Log detalhado
      tissue_cols_2 <- setdiff(colnames(mean_global), "Gene")
      log_message(paste("  ✓", gene, "-", length(tissue_cols_2), 
                       "tecido(s) com expressão"))
      
    }, error = function(e) {
      log_message(paste("  ✗ Erro com gene", gene, ":", e$message))
    })
  }
  
  return(all_gene_2_files)
}

# Função para criar arquivo All_Genes.csv
create_all_genes_file <- function(all_gene_2_files) {
  if (length(all_gene_2_files) == 0) {
    log_message("Nenhum arquivo 2_ para consolidar")
    return(NULL)
  }
  
  # Encontrar todos os tecidos únicos em todos os genes
  all_tissues <- unique(unlist(lapply(all_gene_2_files, function(df) {
    setdiff(colnames(df), "Gene")
  })))
  
  log_message(paste("Encontrados", length(all_tissues), "tecidos únicos"))
  
  # Processar cada gene para garantir que tenha todas as colunas de tecidos
  processed_files <- list()
  
  for (gene in names(all_gene_2_files)) {
    df <- all_gene_2_files[[gene]]
    
    # Garantir que todas as colunas de tecidos estejam presentes
    missing_tissues <- setdiff(all_tissues, colnames(df))
    
    if (length(missing_tissues) > 0) {
      # Adicionar colunas faltantes com NA
      for (tissue in missing_tissues) {
        df[[tissue]] <- NA_real_
      }
    }
    
    # Reordenar colunas: Gene primeiro, depois tecidos em ordem alfabética
    tissue_cols_sorted <- sort(all_tissues)
    df <- df %>% select(Gene, all_of(tissue_cols_sorted))
    
    processed_files[[gene]] <- df
  }
  
  # Combinar todos os dataframes
  all_genes_df <- bind_rows(processed_files)
  
  # Ordenar por Gene
  all_genes_df <- all_genes_df %>% arrange(Gene)
  
  # Salvar arquivo
  output_file <- file.path(gene_expression_dir, "All_Genes.csv")
  write_csv(all_genes_df, output_file, na = "")
  
  log_message(paste("Arquivo All_Genes.csv criado com", 
                   nrow(all_genes_df), "genes e",
                   ncol(all_genes_df) - 1, "tecidos"))
  
  # Salvar lista de tecidos
  tissues_file <- file.path(gene_expression_dir, "tecidos_encontrados.txt")
  writeLines(sort(all_tissues), tissues_file)
  log_message(paste("Lista de tecidos salva em:", tissues_file))
  
  return(all_genes_df)
}

# Função principal
main <- function() {
  log_message("=========================================")
  log_message("ANÁLISE DE MICROARRAY - ESTRUTURA ORGANIZADA v2")
  log_message("=========================================")
  
  # Carregar genes de interesse
  log_message("\n1. CARREGANDO GENES DE INTERESSE")
  target_genes <- load_target_genes()
  if (is.null(target_genes) || length(target_genes) == 0) {
    log_message("ERRO: Nenhum gene carregado")
    return()
  }
  
  # Carregar mapeamento
  log_message("\n2. CARREGANDO MAPEAMENTO")
  probe_mapping <- load_probe_mapping()
  if (is.null(probe_mapping)) {
    log_message("ERRO: Mapeamento não carregado")
    return()
  }
  
  # Listar arquivos
  log_message("\n3. LOCALIZANDO ARQUIVOS")
  files <- list.files(experiment_dir, pattern = "\\.qs$", full.names = TRUE)
  log_message(paste("Encontrados", length(files), "arquivos"))
  
  # Processar arquivos
  log_message("\n4. PROCESSANDO ARQUIVOS")
  all_results <- list()
  successful <- 0
  failed <- 0
  failures <- list()
  
  for (i in seq_along(files)) {
    file <- files[i]
    log_message(paste("\n[", i, "/", length(files), "]", basename(file)))
    
    result <- process_expression_set(file, target_genes, probe_mapping)
    
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
  
  # Relatório
  log_message("\n=========================================")
  log_message("RELATÓRIO FINAL")
  log_message("=========================================")
  log_message(paste("Sucesso:", successful, "/", length(files)))
  log_message(paste("Falha:", failed, "/", length(files)))
  
  if (failed > 0) {
    log_message("\nArquivos com falha:")
    for (f in names(failures)) {
      log_message(paste("  -", f, ":", failures[[f]]))
    }
  }
  
  # Salvar resultados
  if (successful > 0) {
    log_message("\n5. CRIANDO ARQUIVOS POR GENE")
    
    # Combinar todos os dados
    all_data <- bind_rows(all_results)
    
    # Criar arquivos por gene (0_, 1_, 2_)
    all_gene_2_files <- create_gene_files(all_data)
    
    # Criar arquivo consolidado All_Genes.csv
    log_message("\n6. CRIANDO ARQUIVO CONSOLIDADO")
    all_genes_df <- create_all_genes_file(all_gene_2_files)
    
    if (!is.null(all_genes_df)) {
      # Resumo estatístico
      log_message("\n7. RESUMO ESTATÍSTICO")
      
      # Contar tecidos com expressão por gene
      tissues_with_expression <- apply(all_genes_df[, -1], 1, function(row) {
        sum(!is.na(row))
      })
      
      log_message(paste("Média de tecidos com expressão por gene:", 
                       round(mean(tissues_with_expression), 2)))
      log_message(paste("Genes sem expressão em nenhum tecido:", 
                       sum(tissues_with_expression == 0)))
      log_message(paste("Genes com expressão em todos os tecidos:", 
                       sum(tissues_with_expression == (ncol(all_genes_df) - 1))))
      
      # Contar genes por tecido
      genes_per_tissue <- sapply(all_genes_df[, -1], function(col) {
        sum(!is.na(col))
      })
      
      top_10_tissues <- sort(genes_per_tissue, decreasing = TRUE)[1:min(10, length(genes_per_tissue))]
      log_message("\nTop 10 tecidos com mais genes:")
      for (i in seq_along(top_10_tissues)) {
        log_message(paste("  ", names(top_10_tissues)[i], ":", top_10_tissues[i], "genes"))
      }
    }
    
    # Salvar dados brutos consolidados
    log_message("\n8. SALVANDO DADOS BRUTOS CONSOLIDADOS")
    write_csv(all_data, file.path(gene_expression_dir, "all_expression_data_raw.csv"))
    log_message("Dados brutos salvos em: all_expression_data_raw.csv")
  }
  
  log_message("\n=========================================")
  log_message("ANÁLISE CONCLUÍDA")
  log_message("=========================================")
}

# Executar
log_message("Iniciando análise...")
main()
