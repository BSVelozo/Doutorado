# Script de Classificação de Amostras do Expression Atlas
# Classificação de amostras de próstata em Cancer, Saudavel, HPB, Descartado e Indeciso

# Instalação e carregamento de pacotes necessários
if (!require("qs", quietly = TRUE)) {
  install.packages("qs")
}
if (!require("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!require("stringr", quietly = TRUE)) {
  install.packages("stringr")
}
if (!require("Biobase", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("Biobase")
}
if (!require("SummarizedExperiment", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("SummarizedExperiment")
}

library(qs)
library(dplyr)
library(stringr)
library(Biobase)
library(SummarizedExperiment)

# DEFINIÇÃO DOS PARÂMETROS
pasta_dados <- "."  # Pasta onde estão os arquivos .qs
arquivo_saida_principal <- "experimentos_categorizados.csv"
arquivo_contagem <- "contagem_categorias.csv"
arquivo_frases <- "frases_descritivas.csv"
arquivo_indecisos <- "amostras_indecisas.csv"

# DEFINIÇÃO DOS BANCOS DE PALAVRAS PARA CATEGORIZAÇÃO
keyword_cancer <- c("cancer", "câncer", "adenocarcinoma", "cancro", "carcinoma", 
                   "neoplasia", "tumor", "malign", "maligno", "neoplastic",
                   "prostate cancer", "cancer of prostate", "pc", "pca", "tumour",
                   "malignant", "neoplasm", "prostate carcinoma", "metastatic prostate cancer",
                   "aggressive androgen negative prostate cancer", "androgen independent prostate cancer",
                   "castration-resistant prostate cancer", "prostate adenocarcinoma", "bone metastasis of grade IV prostate cancer")

keyword_saudavel <- c("normal", "saudável", "healthy", "wild type", "controle", 
                     "não tumor", "non.tumor", "benign", "non.tumoral", "non.neoplastic",
                     "healthy prostate", "normal prostate", "non.cancer", "non.malignant",
                     "wildtype", "wt", "control", "non.tumor", "non.malignant", "morphologically normal tissue",
                     "adjacent normal prostate")

keyword_hpb <- c("hpb", "hyperplasia", "hiperplasia", "bph", "benign.prostatic.hyperplasia",
                "benign hyperplasia", "hyperplastic", "hiperplásico", "benign prostate tumor", 
                "benign prostate tumour", "benign tumour", "benign tumor")

# Lista de tecidos não-próstata para descartar
tecidos_nao_prostata <- c("synovial membrane", "cardiac atrium", "cardiac ventricle", 
                         "heart", "lung", "liver", "kidney", "brain", "skin", "breast",
                         "colon", "stomach", "intestine", "pancreas", "spleen", "thyroid",
                         "adrenal", "bladder", "testis", "ovary", "uterus", "placenta",
                         "blood", "bone marrow", "lymph node", "skeletal muscle", "smooth muscle",
                         "nerve", "spinal cord", "cerebellum", "cerebrum", "hippocampus",
                         "retina", "cornea", "lens", "salivary gland", "pituitary", "pineal",
                         "thymus", "tonsil", "appendix", "esophagus", "trachea", "bronchus",
                         "alveolar", "epididymis", "seminal vesicle", "vas deferens", "penis")

# FUNÇÃO PARA CATEGORIZAR UMA DESCRIÇÃO DE AMOSTRA (COM TODAS AS NOVAS REGRAS)
categorize_sample <- function(description) {
  if (is.null(description) || is.na(description) || description == "") {
    return("Descartado")
  }
  
  desc <- as.character(description)
  desc_lower <- tolower(desc)
  
  # REGRA 1: "benign prostatic hyperplasia" -> HPB (com alta prioridade)
  if (grepl("benign prostatic hyperplasia", desc_lower, ignore.case = TRUE)) {
    return("HPB")
  }
  
  # REGRA 2: "benign prostate tumor" -> HPB (com alta prioridade)
  if (grepl("benign prostate tumor", desc_lower, ignore.case = TRUE)) {
    return("HPB")
  }
  
  # REGRA 3: "adjacent non-tumor tissue" -> Saudavel (alta prioridade)
  if (grepl("adjacent non-tumor tissue", desc_lower, ignore.case = TRUE)) {
    return("Saudavel")
  }
  
  # REGRA 4: "| control siRNA |" ou "| control |" entre pipes -> Saudavel
  if (grepl("\\|\\s*control\\s*(sirna)?\\s*\\|", desc_lower, ignore.case = TRUE)) {
    return("Saudavel")
  }
  
  # REGRA 5: "| female |" entre pipes -> Descartado
  if (grepl("\\|\\s*female\\s*\\|", desc_lower, ignore.case = TRUE)) {
    return("Descartado")
  }
  
  # REGRA 6: Amostras com "prostate carcinoma" e "control" juntos -> Saudavel
  if (grepl("prostate carcinoma", desc_lower, ignore.case = TRUE) && 
      grepl("control", desc_lower, ignore.case = TRUE)) {
    return("Saudavel")
  }
  
  # REGRA 7: Células estromais da próstata -> Saudavel
  if (grepl("primary prostate stromal cell", desc_lower, ignore.case = TRUE)) {
    return("Saudavel")
  }
  
  # REGRA 8: Verificar se é tecido não-próstata
  is_tecido_nao_prostata <- any(sapply(tecidos_nao_prostata, 
                                       function(x) grepl(x, desc_lower, ignore.case = TRUE)))
  
  # Verificar se é tecido de próstata
  is_prostata <- grepl("prostate|prostatic|próstata", desc_lower, ignore.case = TRUE)
  
  # Se for tecido não-próstata e não for próstata, descartar
  if (is_tecido_nao_prostata && !is_prostata) {
    return("Descartado")
  }
  
  # Agora aplicar as regras gerais de classificação
  has_cancer <- any(sapply(keyword_cancer, function(x) grepl(x, desc_lower, ignore.case = TRUE)))
  has_saudavel <- any(sapply(keyword_saudavel, function(x) grepl(x, desc_lower, ignore.case = TRUE)))
  has_hpb <- any(sapply(keyword_hpb, function(x) grepl(x, desc_lower, ignore.case = TRUE)))
  
  # REGRA ESPECIAL: Para amostras saudáveis, verificar se são de próstata
  # Se for saudável mas não for de próstata, descartar
  if (has_saudavel && !has_cancer && !has_hpb) {
    if (!is_prostata) {
      return("Descartado")
    }
  }
  
  category_count <- sum(has_cancer, has_saudavel, has_hpb)
  
  if (category_count > 1) {
    return("Indeciso")
  } else if (has_cancer) {
    return("Cancer")
  } else if (has_saudavel) {
    return("Saudavel")
  } else if (has_hpb) {
    return("HPB")
  } else {
    return("Descartado")
  }
}

# FUNÇÕES DE EXTRAÇÃO (MANTIDAS)
extrair_dados_expression_set <- function(expression_set, nome_experimento) {
  pheno_data <- pData(expression_set)
  if (nrow(pheno_data) == 0) return(NULL)
  
  colunas_texto <- sapply(pheno_data, function(x) is.character(x) || is.factor(x))
  if (any(colunas_texto)) {
    pheno_data_texto <- pheno_data
    for (col in which(colunas_texto)) {
      pheno_data_texto[[col]] <- as.character(pheno_data_texto[[col]])
    }
    
    descricoes <- apply(pheno_data_texto[, colunas_texto, drop = FALSE], 1, 
                       function(x) {
                         valores <- na.omit(x)
                         if (length(valores) > 0) paste(valores, collapse = " | ") else ""
                       })
    
    return(data.frame(
      Experimento = nome_experimento,
      Corrida = rownames(pheno_data),
      Frase_Descritiva = descricoes,
      stringsAsFactors = FALSE
    ))
  }
  return(NULL)
}

extrair_dados_ranged_summarized_experiment <- function(rse, nome_experimento) {
  col_data <- colData(rse)
  if (nrow(col_data) == 0) return(NULL)
  
  colunas_texto <- sapply(col_data, function(x) is.character(x) || is.factor(x))
  if (any(colunas_texto)) {
    col_data_texto <- col_data
    for (col in which(colunas_texto)) {
      col_data_texto[[col]] <- as.character(col_data_texto[[col]])
    }
    
    descricoes <- apply(col_data_texto[, colunas_texto, drop = FALSE], 1, 
                       function(x) {
                         valores <- na.omit(x)
                         if (length(valores) > 0) paste(valores, collapse = " | ") else ""
                       })
    
    return(data.frame(
      Experimento = nome_experimento,
      Corrida = rownames(col_data),
      Frase_Descritiva = descricoes,
      stringsAsFactors = FALSE
    ))
  }
  return(NULL)
}

# FUNÇÃO PRINCIPAL
processar_e_classificar_dados <- function(pasta_dados) {
  cat("Iniciando processamento e classificação dos dados do Expression Atlas...\n")
  
  qs_files <- list.files(path = pasta_dados, pattern = "\\.qs$", full.names = TRUE)
  if (length(qs_files) == 0) stop("Nenhum arquivo .qs encontrado")
  
  cat("Encontrados", length(qs_files), "arquivos .qs\n")
  
  todos_dados <- list()
  
  for (file in qs_files) {
    nome_arquivo <- basename(file)
    nome_experimento <- gsub("\\.qs$", "", nome_arquivo)
    
    cat("Processando:", nome_arquivo, "\n")
    
    tryCatch({
      dados <- qs::qread(file)
      dados_extraidos <- NULL
      
      if (inherits(dados, "SimpleList")) {
        for (elemento_name in names(dados)) {
          elemento <- dados[[elemento_name]]
          if (inherits(elemento, "ExpressionSet")) {
            cat("  -> ExpressionSet:", elemento_name, "\n")
            dados_temp <- extrair_dados_expression_set(elemento, nome_experimento)
          } else if (inherits(elemento, "RangedSummarizedExperiment")) {
            cat("  -> RangedSummarizedExperiment:", elemento_name, "\n")
            dados_temp <- extrair_dados_ranged_summarized_experiment(elemento, nome_experimento)
          } else {
            next
          }
          
          if (!is.null(dados_temp)) {
            if (is.null(dados_extraidos)) {
              dados_extraidos <- dados_temp
            } else {
              dados_extraidos <- bind_rows(dados_extraidos, dados_temp)
            }
          }
        }
      } else if (inherits(dados, "ExpressionSet")) {
        cat("  -> ExpressionSet direto\n")
        dados_extraidos <- extrair_dados_expression_set(dados, nome_experimento)
      } else if (inherits(dados, "RangedSummarizedExperiment")) {
        cat("  -> RangedSummarizedExperiment direto\n")
        dados_extraidos <- extrair_dados_ranged_summarized_experiment(dados, nome_experimento)
      }
      
      if (!is.null(dados_extraidos)) {
        todos_dados[[nome_experimento]] <- dados_extraidos
        cat("  -> Extraídas", nrow(dados_extraidos), "amostras\n")
      }
      
    }, error = function(e) {
      cat("  -> Erro:", e$message, "\n")
    })
  }
  
  if (length(todos_dados) == 0) stop("Nenhum dado extraído")
  
  # Combinar dados
  dados_combinados <- bind_rows(todos_dados)
  cat("Total de amostras extraídas:", nrow(dados_combinados), "\n")
  
  # Categorização
  cat("Categorizando amostras com as novas regras...\n")
  dados_combinados$Tipo_Amostra <- vapply(
    dados_combinados$Frase_Descritiva, 
    categorize_sample, 
    character(1),
    USE.NAMES = FALSE
  )
  
  # Garantir que é vetor atômico
  if (!is.atomic(dados_combinados$Tipo_Amostra)) {
    cat("Convertendo Tipo_Amostra para vetor atômico...\n")
    dados_combinados$Tipo_Amostra <- as.character(dados_combinados$Tipo_Amostra)
  }
  
  # Contagem
  cat("Calculando contagens...\n")
  contagem_table <- table(dados_combinados$Tipo_Amostra)
  contagem <- data.frame(
    Tipo_Amostra = names(contagem_table),
    Contagem = as.numeric(contagem_table),
    stringsAsFactors = FALSE
  ) %>% arrange(desc(Contagem))
  
  # CRIAR ARQUIVOS DE SAÍDA
  cat("Criando arquivos de saída...\n")
  
  # 1. Arquivo principal
  write.csv(dados_combinados %>% select(Experimento, Corrida, Tipo_Amostra), 
            arquivo_saida_principal, row.names = FALSE, fileEncoding = "UTF-8")
  
  # 2. Arquivo de contagem
  write.csv(contagem, arquivo_contagem, row.names = FALSE, fileEncoding = "UTF-8")
  
  # 3. Arquivo com frases descritivas
  write.csv(dados_combinados, arquivo_frases, row.names = FALSE, fileEncoding = "UTF-8")
  
  # 4. Separar por categoria e criar arquivos individuais
  categorias <- unique(dados_combinados$Tipo_Amostra)
  for (categoria in categorias) {
    dados_categoria <- dados_combinados %>% filter(Tipo_Amostra == categoria)
    nome_arquivo <- paste0(tolower(categoria), ".csv")
    write.csv(dados_categoria, nome_arquivo, row.names = FALSE, fileEncoding = "UTF-8")
    cat("  ->", nome_arquivo, ":", nrow(dados_categoria), "amostras\n")
  }
  
  # 5. Arquivo específico para indecisos se existir
  if ("Indeciso" %in% categorias) {
    indecisos <- dados_combinados %>% filter(Tipo_Amostra == "Indeciso")
    write.csv(indecisos, arquivo_indecisos, row.names = FALSE, fileEncoding = "UTF-8")
  }
  
  # RELATÓRIO FINAL
  cat("\n=== RELATÓRIO FINAL ===\n")
  cat("Total de experimentos:", length(unique(dados_combinados$Experimento)), "\n")
  cat("Total de amostras:", nrow(dados_combinados), "\n")
  cat("\nDistribuição por categoria:\n")
  print(contagem)
  
  # Mostrar exemplos de cada categoria
  cat("\n=== EXEMPLOS POR CATEGORIA ===\n")
  for (categoria in categorias) {
    cat("\n", categoria, " (primeiras 2 amostras):\n", sep = "")
    exemplos <- dados_combinados %>% 
      filter(Tipo_Amostra == categoria) %>%
      head(2) %>%
      select(Experimento, Corrida, Frase_Descritiva)
    print(exemplos)
  }
  
  cat("\nArquivos gerados:\n")
  cat("-", arquivo_saida_principal, "\n")
  cat("-", arquivo_contagem, "\n") 
  cat("-", arquivo_frases, "\n")
  for (categoria in categorias) {
    cat("-", paste0(tolower(categoria), ".csv"), "\n")
  }
  if ("Indeciso" %in% categorias) {
    cat("-", arquivo_indecisos, "\n")
  }
  
  return(list(
    dados_categorizados = dados_combinados %>% select(Experimento, Corrida, Tipo_Amostra),
    contagem = contagem,
    frases = dados_combinados
  ))
}

# EXECUTAR
tryCatch({
  resultados <- processar_e_classificar_dados(pasta_dados)
  cat("\n=== PROCESSAMENTO CONCLUÍDO ===\n")
  print(resultados$contagem)
}, error = function(e) {
  cat("ERRO:", e$message, "\n")
})