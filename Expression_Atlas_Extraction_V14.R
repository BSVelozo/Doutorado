# Script: Expression_Atlas_Pipeline.R
# Descrição: Pipeline para download e análise de dados de expressão gênica do EMBL-EBI Expression Atlas
# Versão melhorada com tratamento robusto de erros e uso do pacote qs

# Carregamento de pacotes necessários
library(readr)        # Leitura/escrita de arquivos CSV
library(stringr)      # Manipulação de strings
library(ExpressionAtlas) # Pacote principal para acesso ao Expression Atlas
library(SummarizedExperiment) # Manipulação de dados de sequenciamento
library(tidyr)        # Manipulação de dados (pivotagem)
library(dplyr)        # Manipulação de dados (gramática de dados)
library(S4Vectors)    # Estruturas de dados fundamentais do Bioconductor
library(GenomicRanges) # Representação de intervalos genômicos
library(qs)           # Serialização rápida e eficiente de objetos (alternativa ao saveRDS)

# Configurações iniciais de diretórios
experiment_dir <- "Experiment_Files"      # Diretório para armazenar dados brutos
gene_expression_dir <- "Gene_Expression_Files" # Diretório para dados processados
log_file <- "pipeline_log.txt"            # Arquivo de log do pipeline
failed_experiments_file <- "Experiments_Failed.txt" # Registro de experimentos falhos

# Criar diretórios se não existirem
dir.create(experiment_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(gene_expression_dir, showWarnings = FALSE, recursive = TRUE)

# Função de log com timestamp
log_message <- function(msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  full_msg <- paste0("[", timestamp, "] ", msg)
  cat(full_msg, "\n")
  write(full_msg, file = log_file, append = TRUE)
}

# CONFIGURAÇÕES PARA UBUNTU
log_message("Iniciando pipeline do Expression Atlas para Ubuntu 24.04")
options(expressions = 100000)  # Aumenta limite de expressões aninhadas
options(timeout = 300)        # Aumenta timeout para downloads (5 minutos)

# Função para verificar experimentos já baixados
get_downloaded_experiments <- function() {
  if (!dir.exists(experiment_dir)) {
    return(character(0))
  }
  
  # Procurar por arquivos .rds e .qs (suporte a ambos os formatos)
  downloaded_files <- c(
    list.files(experiment_dir, pattern = "\\.rds$", full.names = TRUE),
    list.files(experiment_dir, pattern = "\\.qs$", full.names = TRUE)
  )
  
  if (length(downloaded_files) == 0) {
    return(character(0))
  }
  
  # Verificar se os arquivos são válidos
  valid_files <- character(0)
  for (file in downloaded_files) {
    tryCatch({
      # Verificar extensão do arquivo para usar o método de leitura correto
      if (grepl("\\.rds$", file)) {
        obj <- readRDS(file)
      } else if (grepl("\\.qs$", file)) {
        obj <- qs::qread(file)
      }
      
      # Verifica se é um objeto SimpleList (estrutura esperada do ExpressionAtlas)
      if (inherits(obj, "SimpleList") && length(obj) > 0) {
        valid_files <- c(valid_files, file)
      } else {
        log_message(paste("⚠️  Arquivo inválido encontrado:", file, "- Será ignorado"))
      }
    }, error = function(e) {
      log_message(paste("⚠️  Arquivo corrompido encontrado:", file, "- Será ignorado"))
    })
  }
  
  # Extrair apenas os IDs dos arquivos válidos (remover extensão)
  downloaded_ids <- gsub("\\.(rds|qs)$", "", basename(valid_files))
  return(downloaded_ids)
}

# Função para salvar o objeto complexo de forma robusta usando qs
save_complex_object <- function(obj, file_path) {
  max_attempts <- 3
  file_path_qs <- sub("\\.rds$", ".qs", file_path)  # Usar extensão .qs
  
  for (attempt in 1:max_attempts) {
    tryCatch({
      # Criar diretório se não existir
      dir.create(dirname(file_path_qs), showWarnings = FALSE, recursive = TRUE)
      
      # Salvar o objeto usando qs (mais eficiente para objetos grandes)
      qs::qsave(obj, file = file_path_qs, preset = "high")
      
      # Verificar se o arquivo foi criado e é válido
      if (file.exists(file_path_qs)) {
        file_info <- file.info(file_path_qs)
        if (file_info$size > 1024) {  # Pelo menos 1KB
          # Verificar se pode ser lido novamente
          test_obj <- qs::qread(file_path_qs)
          if (inherits(test_obj, class(obj)[1])) {
            # Se também queremos manter uma versão .rds, podemos converter
            # Mas vamos usar apenas .qs para economizar espaço e evitar problemas
            return(TRUE)
          } else {
            file.remove(file_path_qs)
            return(FALSE)
          }
        } else {
          file.remove(file_path_qs)
          return(FALSE)
        }
      } else {
        return(FALSE)
      }
    }, error = function(e) {
      log_message(paste("❌ Erro ao salvar objeto (tentativa", attempt, "):", e$message))
      if (file.exists(file_path_qs)) {
        file.remove(file_path_qs)
      }
      
      # Tentar com saveRDS como fallback na última tentativa
      if (attempt == max_attempts) {
        tryCatch({
          log_message("Tentando saveRDS como fallback...")
          saveRDS(obj, file = file_path, compress = "gzip")
          if (file.exists(file_path)) {
            return(TRUE)
          }
        }, error = function(e2) {
          log_message(paste("❌ Falha também com saveRDS:", e2$message))
        })
      }
      
      # Pausa exponencial entre tentativas
      Sys.sleep(2 ^ attempt)
      return(FALSE)
    })
  }
  return(FALSE)
}

# Função para baixar UM experimento com verificação robusta
download_single_experiment <- function(acc, max_attempts = 3) {
  for (attempt in 1:max_attempts) {
    tryCatch({
      final_file <- file.path(experiment_dir, paste0(acc, ".rds"))
      
      log_message(paste("Tentativa", attempt, "de", max_attempts, "para baixar", acc))
      
      # Baixar experimento usando getAtlasExperiment (função principal do pacote)
      exp_data <- getAtlasExperiment(acc)
      
      # Verificar se o objeto tem a estrutura esperada
      if (!inherits(exp_data, "SimpleList") || length(exp_data) == 0) {
        log_message(paste("❌ Estrutura de objeto inválida para", acc))
        return(FALSE)
      }
      
      # Tentar salvar o arquivo
      if (save_complex_object(exp_data, final_file)) {
        log_message(paste("✓ Sucesso ao baixar e salvar", acc))
        return(TRUE)
      } else {
        log_message(paste("❌ Falha ao salvar o arquivo para", acc))
        return(FALSE)
      }
      
    }, error = function(e) {
      log_message(paste("❌ Erro na tentativa", attempt, "para", acc, ":", e$message))
      
      # Pausa entre tentativas (backoff exponencial)
      wait_time <- 15 * attempt
      log_message(paste("Aguardando", wait_time, "segundos antes da próxima tentativa"))
      Sys.sleep(wait_time)
    })
  }
  return(FALSE)
}

# PARTE 1: DOWNLOAD DOS EXPERIMENTOS
download_experiments <- function() {
  log_message("=== INICIANDO FASE DE DOWNLOAD ===")
  
  # Verificar experimentos já baixados (com validação)
  downloaded_experiments <- get_downloaded_experiments()
  log_message(paste("Encontrados", length(downloaded_experiments), "experimentos já baixados e válidos"))
  
  # Buscar experimentos Homo sapiens usando searchAtlasExperiments
  log_message("Buscando experimentos Homo sapiens no Expression Atlas...")
  atlas_data <- searchAtlasExperiments(properties = character(0), species = "homo sapiens")
  log_message(paste("Total de experimentos disponíveis:", nrow(atlas_data)))
  
  # Identificar experimentos que precisam ser baixados
  experiments_to_download <- setdiff(atlas_data$Accession, downloaded_experiments)
  log_message(paste(length(experiments_to_download), "experimentos precisam ser baixados"))
  
  # Se não há experimentos para baixar, encerrar
  if (length(experiments_to_download) == 0) {
    log_message("✅ Todos os experimentos já estão baixados")
    return(list(success = 0, failed = 0, actual_files = length(downloaded_experiments)))
  }
  
  # Inicializar lista de experimentos falhos
  failed_experiments <- character(0)
  success_count <- 0
  fail_count <- 0
  
  # Processar experimentos sequencialmente
  for (acc in experiments_to_download) {
    log_message(paste("Processando experimento:", acc))
    
    # Baixar experimento
    if (download_single_experiment(acc, 3)) {
      success_count <- success_count + 1
      log_message(paste("✓", acc, "baixado com sucesso"))
    } else {
      fail_count <- fail_count + 1
      failed_experiments <- c(failed_experiments, acc)
      log_message(paste("❌ Falha ao baixar", acc, "após 3 tentativas"))
    }
    
    # Pausa entre downloads para não sobrecarregar o servidor
    Sys.sleep(2)
  }
  
  # Salvar lista de experimentos falhos
  if (length(failed_experiments) > 0) {
    write_lines(failed_experiments, failed_experiments_file)
    log_message(paste("Lista de", length(failed_experiments), "experimentos falhos salva em", failed_experiments_file))
  }
  
  # Verificar consistência dos downloads
  downloaded_final <- get_downloaded_experiments()
  log_message(paste("✅ Total de arquivos válidos após download:", length(downloaded_final)))
  
  # Relatório final
  log_message("=== RELATÓRIO DE DOWNLOAD ===")
  log_message(paste("Experimentros já baixados anteriormente:", length(downloaded_experiments)))
  log_message(paste("Novos experimentos baixados com sucesso:", success_count))
  log_message(paste("Novos experimentos com falha:", fail_count))
  log_message(paste("Total de arquivos no diretório:", length(downloaded_final)))
  log_message("=== FASE DE DOWNLOAD CONCLUÍDA ===")
  
  return(list(
    previously_downloaded = length(downloaded_experiments),
    success = success_count, 
    failed = fail_count, 
    actual_files = length(downloaded_final)
  ))
}

# PARTE 2: ANÁLISE DE EXPRESSÃO GÊNICA
analyze_expression <- function() {
  log_message("=== INICIANDO FASE DE ANÁLISE DE EXPRESSÃO ===")
  
  # Carregar experimentos baixados
  log_message("Carregando experimentos baixados")
  experiment_files <- c(
    list.files(experiment_dir, pattern = "\\.rds$", full.names = TRUE),
    list.files(experiment_dir, pattern = "\\.qs$", full.names = TRUE)
  )
  
  if (length(experiment_files) == 0) {
    log_message("❌ Nenhum experimento encontrado. Execute primeiro a fase de download.")
    return()
  }
  
  # Processar experimentos sequencialmente
  all_expression_data <- list()
  
  for (file in experiment_files) {
    acc <- gsub("\\.(rds|qs)$", "", basename(file))
    log_message(paste("Processando experimento:", acc))
    
    tryCatch({
      # Carregar o objeto (suporte a .rds e .qs)
      if (grepl("\\.rds$", file)) {
        exp_data <- readRDS(file)
      } else if (grepl("\\.qs$", file)) {
        exp_data <- qs::qread(file)
      }
      
      # Verificar se é um objeto SimpleList válido
      if (!inherits(exp_data, "SimpleList") || length(exp_data) == 0) {
        log_message(paste("⚠️  Objeto inválido em", file, "- Pulando"))
        next
      }
      
      # Extrair dados de RNA-seq (SummarizedExperiment)
      if ("rnaseq" %in% names(exp_data)) {
        rse <- exp_data$rnaseq
        
        # Extrair matriz de expressão (counts)
        if ("counts" %in% assayNames(rse)) {
          expr_matrix <- assay(rse, "counts")
          
          # Extrair metadados das amostras (colData)
          sample_metadata <- as.data.frame(colData(rse))
          
          # Preparar dados para saída (formato tidy)
          expr_df <- as.data.frame(expr_matrix) %>%
            rownames_to_column("GeneID") %>%
            pivot_longer(cols = -GeneID, names_to = "Sample", values_to = "Expression") %>%
            mutate(Experiment = acc, ValueType = "counts")
          
          # Adicionar metadados das amostras se disponíveis
          if ("organism_part" %in% colnames(sample_metadata)) {
            sample_meta <- sample_metadata %>%
              rownames_to_column("Sample") %>%
              select(Sample, organism_part)
            
            expr_df <- left_join(expr_df, sample_meta, by = "Sample")
          } else {
            expr_df$organism_part <- "Desconhecido"
          }
          
          all_expression_data[[acc]] <- expr_df
        }
      }
      
      log_message(paste("✓ Experimentos", acc, "processado"))
      
    }, error = function(e) {
      log_message(paste("❌ Erro ao processar experimento", acc, ":", e$message))
    })
    
    # Liberar memória
    gc(full = TRUE)
  }
  
  # Combinar todos os dados de expressão
  if (length(all_expression_data) > 0) {
    combined_expression <- bind_rows(all_expression_data)
    
    # Salvar dados consolidados
    log_message("Salvando dados consolidados de expressão")
    write_csv(combined_expression, "expression_all_genes.csv")
    
    log_message(paste("✅ Dados de expressão salvos para", nrow(combined_expression), "registros"))
  } else {
    log_message("❌ Nenhum dado de expressão encontrado para os genes de interesse")
  }
  
  log_message("=== FASE DE ANÁLISE CONCLUÍDA ===")
}

# EXECUÇÃO DO PIPELINE
log_message("Iniciando pipeline do Expression Atlas para Ubuntu 24.04")

# Fase 1: Download dos experimentos
resultados <- download_experiments()

# Fase 2: Análise de expressão (executar apenas após download concluído)
#analyze_expression()

log_message("Pipeline concluído com sucesso")