# Carregar o pacote necessário
library(qs)

# Função para explorar a estrutura do arquivo .qs
explore_expression_data <- function(file_path) {
  # Ler o arquivo .qs
  exp_data <- qs::qread(file_path)
  
  cat("=== ESTRUTURA DO OBJETO ===\n")
  cat("Classe do objeto:", paste(class(exp_data), collapse = ", "), "\n")
  cat("Tipo do objeto:", typeof(exp_data), "\n")
  cat("Comprimento do objeto:", length(exp_data), "\n\n")
  
  # Verificar se é um ExpressionSet
  if (is(exp_data, "ExpressionSet")) {
    cat("✅ É um ExpressionSet\n")
    
    # Dimensões da matriz de expressão
    cat("Dimensões da matriz de expressão:", paste(dim(exp_data), collapse = " x "), "\n\n")
    
    # Acessar a matriz de expressão
    if (require(Biobase)) {
      expr_matrix <- Biobase::exprs(exp_data)
      cat("=== MATRIZ DE EXPRESSÃO ===\n")
      cat("Dimensões:", paste(dim(expr_matrix), collapse = " x "), "\n")
      cat("Primeiras 10 linhas e 5 colunas:\n")
      print(expr_matrix[1:10, 1:5])
      cat("\n")
      
      # Verificar IDs do ENSEMBL nos nomes das linhas
      cat("=== IDs DO ENSEMBL (rownames) ===\n")
      ensembl_ids <- rownames(expr_matrix)
      cat("Total de IDs:", length(ensembl_ids), "\n")
      cat("Primeiros 20 IDs:\n")
      print(head(ensembl_ids, 20))
      cat("\n")
      
      # Procurar padrão ENSEMBL nos IDs
      ensembl_pattern <- "ENSG[0-9]+"  # Padrão para IDs de gene humano do ENSEMBL
      ensembl_matches <- grepl(ensembl_pattern, ensembl_ids)
      cat("IDs que correspondem ao padrão ENSEMBL:", sum(ensembl_matches), "\n")
      
      if (sum(ensembl_matches) > 0) {
        cat("Exemplos de IDs do ENSEMBL encontrados:\n")
        print(head(ensembl_ids[ensembl_matches], 10))
      }
      cat("\n")
    }
    
    # Acessar featureData (informações sobre os genes)
    if (require(Biobase)) {
      feature_data <- Biobase::fData(exp_data)
      if (!is.null(feature_data) && ncol(feature_data) > 0) {
        cat("=== FEATURE DATA (informações dos genes) ===\n")
        cat("Dimensões do featureData:", paste(dim(feature_data), collapse = " x "), "\n")
        cat("Colunas disponíveis:", paste(colnames(feature_data), collapse = ", "), "\n\n")
        
        # Procurar colunas que possam conter IDs do ENSEMBL
        ensembl_cols <- character(0)
        for (col_name in colnames(feature_data)) {
          if (any(grepl(ensembl_pattern, feature_data[[col_name]]))) {
            ensembl_cols <- c(ensembl_cols, col_name)
          }
        }
        
        if (length(ensembl_cols) > 0) {
          cat("Colunas que contêm IDs do ENSEMBL:", paste(ensembl_cols, collapse = ", "), "\n")
          for (col in ensembl_cols) {
            cat("Coluna:", col, "\n")
            ensembl_values <- feature_data[[col]][grepl(ensembl_pattern, feature_data[[col]])]
            cat("Primeiros 10 valores do ENSEMBL:\n")
            print(head(ensembl_values, 10))
            cat("\n")
          }
        }
      }
    }
    
    # Acessar phenoData (informações sobre as amostras)
    if (require(Biobase)) {
      pheno_data <- Biobase::pData(exp_data)
      if (!is.null(pheno_data) && ncol(pheno_data) > 0) {
        cat("=== PHENO DATA (informações das amostras) ===\n")
        cat("Dimensões do phenoData:", paste(dim(pheno_data), collapse = " x "), "\n")
        cat("Colunas disponíveis:", paste(colnames(pheno_data), collapse = ", "), "\n")
        cat("Primeiras 5 amostras:\n")
        print(head(pheno_data, 5))
        cat("\n")
      }
    }
    
  } else if (inherits(exp_data, "SimpleList")) {
    cat("✅ É uma SimpleList\n")
    cat("Nomes dos elementos:", paste(names(exp_data), collapse = ", "), "\n\n")
    
    # Explorar cada elemento da SimpleList
    for (elem_name in names(exp_data)) {
      cat("=== ELEMENTO:", elem_name, "===\n")
      elem <- exp_data[[elem_name]]
      cat("Classe:", paste(class(elem), collapse = ", "), "\n")
      cat("Tipo:", typeof(elem), "\n")
      
      if (is(elem, "ExpressionSet")) {
        cat("🎯 Este elemento é um ExpressionSet!\n")
        # Chamar recursivamente a função para este ExpressionSet
        # Criar um ambiente temporário para não modificar o objeto original
        temp_env <- new.env()
        temp_env$temp_data <- elem
        explore_expression_data("temp_data")
        rm(temp_env)
      } else if (is.matrix(elem) || is.data.frame(elem)) {
        cat("Dimensões:", paste(dim(elem), collapse = " x "), "\n")
        cat("Primeiras 5 linhas e 3 colunas:\n")
        print(elem[1:5, 1:3])
        cat("\n")
        
        # Verificar IDs do ENSEMBL
        ensembl_pattern <- "ENSG[0-9]+"
        if (!is.null(rownames(elem))) {
          ensembl_matches <- grepl(ensembl_pattern, rownames(elem))
          cat("IDs do ENSEMBL nos rownames:", sum(ensembl_matches), "\n")
          if (sum(ensembl_matches) > 0) {
            cat("Exemplos:\n")
            print(head(rownames(elem)[ensembl_matches], 10))
          }
        }
      }
      cat("\n")
    }
  }
  
  return(exp_data)
}

# Função específica para extrair IDs do ENSEMBL
extract_ensembl_ids <- function(file_path) {
  exp_data <- qs::qread(file_path)
  ensembl_ids <- character(0)
  
  if (is(exp_data, "ExpressionSet")) {
    # IDs dos rownames da matriz de expressão
    expr_matrix <- Biobase::exprs(exp_data)
    potential_ids <- rownames(expr_matrix)
    
    # Padrão para IDs do ENSEMBL (gene humano)
    ensembl_pattern <- "^ENSG[0-9]+"
    ensembl_ids <- grep(ensembl_pattern, potential_ids, value = TRUE)
    
    # Se não encontrou nos rownames, procurar no featureData
    if (length(ensembl_ids) == 0) {
      feature_data <- Biobase::fData(exp_data)
      for (col in colnames(feature_data)) {
        if (any(grepl(ensembl_pattern, feature_data[[col]]))) {
          ensembl_ids <- c(ensembl_ids, feature_data[[col]][grepl(ensembl_pattern, feature_data[[col]])])
        }
      }
    }
    
  } else if (inherits(exp_data, "SimpleList")) {
    # Procurar em todos os elementos da SimpleList
    for (elem_name in names(exp_data)) {
      elem <- exp_data[[elem_name]]
      
      if (is(elem, "ExpressionSet")) {
        expr_matrix <- Biobase::exprs(elem)
        potential_ids <- rownames(expr_matrix)
        ensembl_pattern <- "^ENSG[0-9]+"
        ensembl_ids <- c(ensembl_ids, grep(ensembl_pattern, potential_ids, value = TRUE))
      }
    }
  }
  
  cat("=== RESUMO DOS IDs DO ENSEMBL ===\n")
  cat("Total de IDs únicos encontrados:", length(unique(ensembl_ids)), "\n")
  cat("Primeiros 20 IDs:\n")
  print(head(unique(ensembl_ids), 20))
  
  return(unique(ensembl_ids))
}

# USO DO CÓDIGO:
# Substitua "caminho/para/seu/arquivo.qs" pelo caminho real do seu arquivo

# Explorar a estrutura completa do arquivo
file_path <- "/home/bsvelozo/Expression_Atlas/Experiment_Files/E-GEOD-11100.qs"  # ALTERE ESTE CAMINHO
exp_data <- explore_expression_data(file_path)

# Extrair apenas os IDs do ENSEMBL
ensembl_ids <- extract_ensembl_ids(file_path)