# Carregar pacotes
library(biomaRt)
library(dplyr)

# Configurações
options(timeout = 600)
Sys.setenv(TZ = "UTC")

# 1. Ler o arquivo .txt com os IDs UniProt
uniprot_ids <- readLines("3_IDs_uniprot.txt")
total_ids <- length(uniprot_ids)

cat("Iniciando mapeamento de", total_ids, "IDs UniProt para Ensembl...\n")

# 2. Conectar ao Ensembl
connect_ensembl_fast <- function() {
  tryCatch({
    ensembl <- useEnsembl(biomart = "genes", 
                         dataset = "hsapiens_gene_ensembl",
                         host = "https://www.ensembl.org")
    return(ensembl)
  }, error = function(e) {
    cat("Erro de conexão:", e$message, "\n")
    return(NULL)
  })
}

ensembl <- connect_ensembl_fast()

# 3. Função para buscar gene, transcrito e proteína
search_ensembl_complete <- function(uniprot_batch, mart) {
  if (is.null(mart)) return(NULL)
  
  tryCatch({
    result <- getBM(
      attributes = c("uniprotswissprot", 
                    "ensembl_gene_id", 
                    "ensembl_transcript_id", 
                    "ensembl_peptide_id"),
      filters = "uniprotswissprot",
      values = uniprot_batch,
      mart = mart
    )
    
    if (nrow(result) > 0) {
      result_df <- data.frame(
        UNIPROT = result$uniprotswissprot,
        `ENSEMBL(Gene)` = result$ensembl_gene_id,
        `ENSEMBL(Transcrito)` = result$ensembl_transcript_id,
        `ENSEMBL(Proteína)` = result$ensembl_peptide_id,
        Origem = "ENSEMBL",
        stringsAsFactors = FALSE,
        check.names = FALSE  # Mantém os nomes com parênteses
      )
      return(result_df)
    }
    return(NULL)
  }, error = function(e) {
    cat("Erro no lote:", e$message, "\n")
    return(NULL)
  })
}

# 4. Função para buscar com filtro uniprotsptrembl
search_ensembl_trembl <- function(uniprot_batch, mart) {
  if (is.null(mart)) return(NULL)
  
  tryCatch({
    result <- getBM(
      attributes = c("uniprotsptrembl", 
                    "ensembl_gene_id", 
                    "ensembl_transcript_id", 
                    "ensembl_peptide_id"),
      filters = "uniprotsptrembl",
      values = uniprot_batch,
      mart = mart
    )
    
    if (nrow(result) > 0) {
      result_df <- data.frame(
        UNIPROT = result$uniprotsptrembl,
        `ENSEMBL(Gene)` = result$ensembl_gene_id,
        `ENSEMBL(Transcrito)` = result$ensembl_transcript_id,
        `ENSEMBL(Proteína)` = result$ensembl_peptide_id,
        Origem = "ENSEMBL",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      return(result_df)
    }
    return(NULL)
  }, error = function(e) {
    cat("Erro no lote TrEMBL:", e$message, "\n")
    return(NULL)
  })
}

# 5. Processar todos os IDs
cat("Processando IDs...\n")
all_results <- data.frame()

# Primeiro tenta SwissProt
result_swissprot <- search_ensembl_complete(uniprot_ids, ensembl)
if (!is.null(result_swissprot)) {
  all_results <- bind_rows(all_results, result_swissprot)
}

# IDs restantes tentar com TrEMBL
found_ids <- all_results$UNIPROT
remaining_ids <- setdiff(uniprot_ids, found_ids)

if (length(remaining_ids) > 0) {
  result_trembl <- search_ensembl_trembl(remaining_ids, ensembl)
  if (!is.null(result_trembl)) {
    all_results <- bind_rows(all_results, result_trembl)
  }
}

# 6. Processar resultados
if (nrow(all_results) > 0) {
  # Remover duplicatas
  all_results <- all_results[!duplicated(all_results$UNIPROT), ]
  
  # Ordenar
  all_results <- all_results[order(all_results$UNIPROT), ]
  
  # Salvar com nomes de colunas exatos
  colnames(all_results) <- c("UNIPROT", "ENSEMBL(Gene)", "ENSEMBL(Transcrito)", 
                            "ENSEMBL(Proteína)", "Origem")
  
  write.table(all_results, "uniprot_ensembl_found.tsv", sep = "\t", 
              row.names = FALSE, quote = FALSE, na = "")
  
  cat("\n=== PROCESSAMENTO CONCLUÍDO ===\n")
  cat("Arquivo gerado: uniprot_ensembl_found.tsv\n")
  cat("Total processado:", total_ids, "\n")
  cat("Encontrados:", nrow(all_results), "\n")
  cat("Taxa de sucesso:", round(nrow(all_results)/total_ids*100, 1), "%\n")
  
} else {
  # Criar arquivo vazio com cabeçalho correto
  empty_df <- data.frame(
    UNIPROT = character(),
    `ENSEMBL(Gene)` = character(),
    `ENSEMBL(Transcrito)` = character(),
    `ENSEMBL(Proteína)` = character(),
    Origem = character(),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  colnames(empty_df) <- c("UNIPROT", "ENSEMBL(Gene)", "ENSEMBL(Transcrito)", 
                         "ENSEMBL(Proteína)", "Origem")
  
  write.table(empty_df, "uniprot_ensembl_found.tsv", sep = "\t", 
              row.names = FALSE, quote = FALSE, na = "")
  
  cat("\n=== PROCESSAMENTO CONCLUÍDO ===\n")
  cat("Nenhuma correspondência encontrada.\n")
}