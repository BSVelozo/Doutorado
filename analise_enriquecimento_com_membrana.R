# Instalar pacotes necessários (se ainda não estiverem instalados)
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("readr")) install.packages("readr")
if (!require("stringr")) install.packages("stringr")

# Carregar pacotes
library(dplyr)
library(tidyr)
library(readr)
library(stringr)

# Função para criar relatório de estatísticas
criar_relatorio_estatisticas <- function(resultados, dados_completos, genes_classificados) {
  
  cat("Gerando relatório de estatísticas...\n")
  
  # Criar data frame para estatísticas
  estatisticas <- data.frame(
    Categoria = character(),
    Gene = character(),
    Tecido_Enriquecido = character(),
    stringsAsFactors = FALSE
  )
  
  # Adicionar dados para cada categoria
  categorias <- c("Enriquecido10X", "Enriquecido5X", "Enriquecido2X")
  
  for (cat in categorias) {
    genes <- resultados[[cat]]$genes
    tecidos <- resultados[[cat]]$tecidos
    
    if (length(genes) > 0) {
      # Criar data frame para cada categoria
      df_cat <- data.frame(
        Categoria = rep(cat, length(genes)),
        Gene = genes,
        Tecido_Enriquecido = tecidos,
        stringsAsFactors = FALSE
      )
      
      estatisticas <- rbind(estatisticas, df_cat)
    }
  }
  
  # Salvar estatísticas detalhadas
  if (nrow(estatisticas) > 0) {
    write_csv(estatisticas, "estatisticas_enriquecimento_detalhado.csv")
    
    # Criar resumo por categoria
    resumo <- estatisticas %>%
      group_by(Categoria) %>%
      summarise(
        Numero_Genes = n(),
        Lista_Genes = paste(Gene, collapse = "; "),
        Tecidos_Unicos = paste(unique(Tecido_Enriquecido), collapse = "; ")
      )
    
    write_csv(resumo, "resumo_estatisticas_enriquecimento.csv")
  }
  
  # Criar relatório em formato texto
  sink("relatorio_enriquecimento.txt")
  cat("===========================================\n")
  cat("RELATÓRIO DE ANÁLISE DE ENRIQUECIMENTO\n")
  cat("===========================================\n\n")
  
  for (cat in categorias) {
    cat("\n", cat, ":\n")
    cat("-------------------\n")
    genes <- resultados[[cat]]$genes
    tecidos <- resultados[[cat]]$tecidos
    
    cat("Número de genes:", length(genes), "\n")
    
    if (length(genes) > 0) {
      cat("\nDistribuição por tecido:\n")
      tab <- table(tecidos)
      for (tecido_nome in names(tab)) {
        cat("  ", tecido_nome, ": ", tab[tecido_nome], " genes\n", sep = "")
      }
    } else {
      cat("Nenhum gene encontrado nesta categoria.\n")
    }
  }
  
  cat("\n===========================================\n")
  cat("RESUMO GERAL\n")
  cat("===========================================\n")
  cat("Total de genes analisados:", nrow(dados_completos), "\n")
  cat("Genes Enriquecido10X:", length(resultados$Enriquecido10X$genes), "\n")
  cat("Genes Enriquecido5X:", length(resultados$Enriquecido5X$genes), "\n")
  cat("Genes Enriquecido2X:", length(resultados$Enriquecido2X$genes), "\n")
  cat("Genes não classificados:", nrow(dados_completos) - length(genes_classificados), "\n")
  
  sink()
  
  cat("\nRelatórios salvos:\n")
  cat("- estatisticas_enriquecimento_detalhado.csv\n")
  cat("- resumo_estatisticas_enriquecimento.csv\n")
  cat("- relatorio_enriquecimento.txt\n")
}

# Função para criar diretórios e salvar dados por tecido
criar_estrutura_tecido <- function(categoria_nome, genes_categoria, tecidos_categoria, dados_completos) {
  if (length(genes_categoria) == 0) {
    return()
  }
  
  # Criar diretório principal da categoria
  dir_principal <- categoria_nome
  if (!dir.exists(dir_principal)) {
    dir.create(dir_principal, recursive = TRUE)
    cat("Diretório criado:", dir_principal, "\n")
  }
  
  # Criar dataframe com todos os genes da categoria
  df_categoria <- dados_completos %>%
    filter(Gene %in% genes_categoria) %>%
    mutate(Categoria = categoria_nome,
           Tecido_Enriquecido = tecidos_categoria[match(Gene, genes_categoria)])
  
  # Salvar arquivo com todos os genes da categoria
  write_csv(df_categoria, file.path(dir_principal, paste0("todos_genes_", tolower(categoria_nome), ".csv")))
  
  # Organizar por tecido
  tecidos_unicos <- unique(tecidos_categoria)
  
  for (tecido in tecidos_unicos) {
    # Limpar nome do tecido para usar em diretórios
    tecido_limpo <- gsub("[^a-zA-Z0-9._-]", "_", tecido)
    
    # Criar diretório do tecido
    dir_tecido <- file.path(dir_principal, tecido_limpo)
    if (!dir.exists(dir_tecido)) {
      dir.create(dir_tecido, recursive = TRUE)
    }
    
    # Filtrar genes deste tecido
    genes_tecido <- genes_categoria[tecidos_categoria == tecido]
    df_tecido <- dados_completos %>%
      filter(Gene %in% genes_tecido)
    
    # Salvar arquivo do tecido
    if (nrow(df_tecido) > 0) {
      nome_arquivo <- paste0("genes_", tolower(categoria_nome), "_", tecido_limpo, ".csv")
      write_csv(df_tecido, file.path(dir_tecido, nome_arquivo))
    }
  }
}

# Função para filtrar por proteínas de membrana e criar estrutura separada
filtrar_por_membrana <- function(categoria_nome, genes_categoria, tecidos_categoria, dados_completos, genes_membrana) {
  if (length(genes_categoria) == 0) {
    return(list(filtrados = 0, total = 0, distribuicao = NULL))
  }
  
  # Filtrar genes que estão na lista de proteínas de membrana
  indices_membrana <- which(genes_categoria %in% genes_membrana)
  genes_membrana_cat <- genes_categoria[indices_membrana]
  tecidos_membrana_cat <- tecidos_categoria[indices_membrana]
  
  if (length(genes_membrana_cat) == 0) {
    return(list(filtrados = 0, total = length(genes_categoria), distribuicao = NULL))
  }
  
  # Criar diretório principal da categoria de membrana
  dir_principal_membrana <- paste0(categoria_nome, "_Membrane")
  if (!dir.exists(dir_principal_membrana)) {
    dir.create(dir_principal_membrana, recursive = TRUE)
    cat("Diretório criado:", dir_principal_membrana, "\n")
  }
  
  # Criar dataframe com todos os genes da categoria (membrana)
  df_categoria_membrana <- dados_completos %>%
    filter(Gene %in% genes_membrana_cat) %>%
    mutate(Categoria = paste0(categoria_nome, "_Membrane"),
           Tecido_Enriquecido = tecidos_membrana_cat[match(Gene, genes_membrana_cat)])
  
  # Salvar arquivo com todos os genes da categoria (membrana)
  write_csv(df_categoria_membrana, file.path(dir_principal_membrana, 
                                            paste0("todos_genes_", tolower(categoria_nome), "_membrane.csv")))
  
  # Organizar por tecido
  tecidos_unicos_membrana <- unique(tecidos_membrana_cat)
  
  for (tecido in tecidos_unicos_membrana) {
    # Limpar nome do tecido para usar em diretórios
    tecido_limpo <- gsub("[^a-zA-Z0-9._-]", "_", tecido)
    
    # Criar diretório do tecido
    dir_tecido <- file.path(dir_principal_membrana, tecido_limpo)
    if (!dir.exists(dir_tecido)) {
      dir.create(dir_tecido, recursive = TRUE)
    }
    
    # Filtrar genes deste tecido
    genes_tecido <- genes_membrana_cat[tecidos_membrana_cat == tecido]
    df_tecido <- dados_completos %>%
      filter(Gene %in% genes_tecido)
    
    # Salvar arquivo do tecido
    if (nrow(df_tecido) > 0) {
      nome_arquivo <- paste0("genes_", tolower(categoria_nome), "_membrane_", tecido_limpo, ".csv")
      write_csv(df_tecido, file.path(dir_tecido, nome_arquivo))
    }
  }
  
  # Calcular distribuição por tecido
  distribuicao <- table(tecidos_membrana_cat)
  
  return(list(
    filtrados = length(genes_membrana_cat),
    total = length(genes_categoria),
    distribuicao = distribuicao
  ))
}

# Criar função para calcular enriquecimento
analisar_enriquecimento <- function() {
  
  # Carregar dados
  cat("Carregando dados...\n")
  dados_saudaveis <- read_csv("All_Genes_Healthy_Tissues.csv")
  dados_cancer <- read_csv("All_Genes_Final_Cancer_Tissue.csv")
  
  # Verificar se os arquivos foram carregados corretamente
  cat("Dimensões dos dados saudáveis:", dim(dados_saudaveis), "\n")
  cat("Dimensões dos dados de câncer:", dim(dados_cancer), "\n")
  
  # Juntar os dados
  dados_completos <- dados_saudaveis %>%
    full_join(dados_cancer, by = "Gene") %>%
    distinct(Gene, .keep_all = TRUE)
  
  # Obter lista de todos os tecidos (excluindo a coluna Gene)
  tecidos <- colnames(dados_completos)[-1]
  
  # Criar listas para armazenar resultados
  resultados <- list(
    Enriquecido10X = list(genes = character(), tecidos = character()),
    Enriquecido5X = list(genes = character(), tecidos = character()),
    Enriquecido2X = list(genes = character(), tecidos = character())
  )
  
  # Vetor para controlar genes já classificados
  genes_classificados <- character()
  
  # Função para calcular se um gene é enriquecido em um tecido
  verificar_enriquecimento <- function(gene_row, tecido_alvo) {
    valor_alvo <- gene_row[[tecido_alvo]]
    
    # Se o valor for NA, não é enriquecido
    if (is.na(valor_alvo)) return(0)
    
    # Obter valores dos outros tecidos
    outros_tecidos <- tecidos[tecidos != tecido_alvo]
    outros_valores <- gene_row[outros_tecidos]
    
    # Remover NAs e converter para numérico
    outros_valores <- suppressWarnings(as.numeric(outros_valores))
    outros_valores <- outros_valores[!is.na(outros_valores)]
    
    # Se não houver outros valores para comparar
    if (length(outros_valores) == 0) return(0)
    
    # Calcular razão mínima em relação a outros tecidos
    razoes <- valor_alvo / outros_valores
    razao_minima <- min(razoes, na.rm = TRUE)
    
    # Verificar categoria
    if (razao_minima >= 10) {
      return(10)
    } else if (razao_minima >= 5) {
      return(5)
    } else if (razao_minima >= 2) {
      return(2)
    } else {
      return(0)
    }
  }
  
  cat("Analisando enriquecimento para", nrow(dados_completos), "genes...\n")
  pb <- txtProgressBar(min = 0, max = nrow(dados_completos), style = 3)
  
  # Analisar cada gene
  for (i in 1:nrow(dados_completos)) {
    setTxtProgressBar(pb, i)
    
    gene <- dados_completos$Gene[i]
    gene_row <- dados_completos[i, , drop = FALSE]
    
    # Pular se o gene já foi classificado
    if (gene %in% genes_classificados) next
    
    # Verificar enriquecimento em cada tecido
    enriquecimentos <- list()
    
    for (tecido in tecidos) {
      categoria <- verificar_enriquecimento(gene_row, tecido)
      if (categoria > 0) {
        enriquecimentos[[tecido]] <- categoria
      }
    }
    
    # Se o gene tem algum enriquecimento
    if (length(enriquecimentos) > 0) {
      # Encontrar a categoria mais alta (prioridade: 10X > 5X > 2X)
      categorias <- unlist(enriquecimentos)
      max_categoria <- max(categorias)
      
      # Encontrar o tecido com a categoria mais alta
      tecidos_max <- names(enriquecimentos)[categorias == max_categoria]
      tecido_max <- tecidos_max[1]  # Em caso de empate, escolhe o primeiro
      
      # Classificar o gene
      if (max_categoria == 10) {
        resultados$Enriquecido10X$genes <- c(resultados$Enriquecido10X$genes, gene)
        resultados$Enriquecido10X$tecidos <- c(resultados$Enriquecido10X$tecidos, tecido_max)
      } else if (max_categoria == 5) {
        resultados$Enriquecido5X$genes <- c(resultados$Enriquecido5X$genes, gene)
        resultados$Enriquecido5X$tecidos <- c(resultados$Enriquecido5X$tecidos, tecido_max)
      } else if (max_categoria == 2) {
        resultados$Enriquecido2X$genes <- c(resultados$Enriquecido2X$genes, gene)
        resultados$Enriquecido2X$tecidos <- c(resultados$Enriquecido2X$tecidos, tecido_max)
      }
      
      genes_classificados <- c(genes_classificados, gene)
    }
  }
  
  close(pb)
  
  # Criar estrutura de pastas por tecido para cada categoria
  cat("\nCriando estrutura de pastas por tecido...\n")
  
  criar_estrutura_tecido("Enriquecido10X", 
                         resultados$Enriquecido10X$genes, 
                         resultados$Enriquecido10X$tecidos, 
                         dados_completos)
  
  criar_estrutura_tecido("Enriquecido5X", 
                         resultados$Enriquecido5X$genes, 
                         resultados$Enriquecido5X$tecidos, 
                         dados_completos)
  
  criar_estrutura_tecido("Enriquecido2X", 
                         resultados$Enriquecido2X$genes, 
                         resultados$Enriquecido2X$tecidos, 
                         dados_completos)
  
  # Carregar lista de proteínas de membrana
  cat("\nCarregando lista de proteínas de membrana...\n")
  if (file.exists("ENSG_Membrane_Proteins.txt")) {
    genes_membrana <- read_lines("ENSG_Membrane_Proteins.txt")
    genes_membrana <- genes_membrana[genes_membrana != ""]  # Remover linhas vazias
    cat("Total de proteínas de membrana carregadas:", length(genes_membrana), "\n")
    
    # Filtrar por proteínas de membrana para cada categoria
    cat("\nFiltrando por proteínas de membrana...\n")
    
    resultados_membrana <- list()
    
    resultados_membrana$Enriquecido10X <- filtrar_por_membrana(
      "Enriquecido10X", 
      resultados$Enriquecido10X$genes, 
      resultados$Enriquecido10X$tecidos, 
      dados_completos, 
      genes_membrana
    )
    
    resultados_membrana$Enriquecido5X <- filtrar_por_membrana(
      "Enriquecido5X", 
      resultados$Enriquecido5X$genes, 
      resultados$Enriquecido5X$tecidos, 
      dados_completos, 
      genes_membrana
    )
    
    resultados_membrana$Enriquecido2X <- filtrar_por_membrana(
      "Enriquecido2X", 
      resultados$Enriquecido2X$genes, 
      resultados$Enriquecido2X$tecidos, 
      dados_completos, 
      genes_membrana
    )
    
    # Criar relatório para proteínas de membrana
    sink("relatorio_membrana_enriquecimento.txt")
    cat("===========================================\n")
    cat("RELATÓRIO DE PROTEÍNAS DE MEMBRANA ENRIQUECIDAS\n")
    cat("===========================================\n\n")
    
    for (categoria in c("Enriquecido10X", "Enriquecido5X", "Enriquecido2X")) {
      cat("\n", categoria, ":\n")
      cat("-------------------\n")
      
      res <- resultados_membrana[[categoria]]
      cat("Total de genes na categoria:", res$total, "\n")
      cat("Proteínas de membrana encontradas:", res$filtrados, "\n")
      cat("Porcentagem: ", round(res$filtrados/res$total * 100, 2), "%\n", sep = "")
      
      if (!is.null(res$distribuicao) && length(res$distribuicao) > 0) {
        cat("\nDistribuição por tecido das proteínas de membrana:\n")
        for (tecido_nome in names(res$distribuicao)) {
          cat("  ", tecido_nome, ": ", res$distribuicao[tecido_nome], " genes\n", sep = "")
        }
      } else if (res$filtrados > 0) {
        cat("\nNenhuma distribuição por tecido disponível.\n")
      }
    }
    
    # Resumo geral
    cat("\n===========================================\n")
    cat("RESUMO GERAL - PROTEÍNAS DE MEMBRANA\n")
    cat("===========================================\n")
    
    total_membrana <- sum(
      resultados_membrana$Enriquecido10X$filtrados,
      resultados_membrana$Enriquecido5X$filtrados,
      resultados_membrana$Enriquecido2X$filtrados
    )
    
    cat("Total de proteínas de membrana enriquecidas:", total_membrana, "\n")
    cat("Proteínas de membrana Enriquecido10X:", resultados_membrana$Enriquecido10X$filtrados, "\n")
    cat("Proteínas de membrana Enriquecido5X:", resultados_membrana$Enriquecido5X$filtrados, "\n")
    cat("Proteínas de membrana Enriquecido2X:", resultados_membrana$Enriquecido2X$filtrados, "\n")
    
    sink()
    
    cat("\nRelatório de proteínas de membrana salvo: relatorio_membrana_enriquecimento.txt\n")
    
  } else {
    cat("Arquivo 'ENSG_Membrane_Proteins.txt' não encontrado. Pulando filtragem por proteínas de membrana.\n")
  }
  
  # Criar relatório de estatísticas
  criar_relatorio_estatisticas(resultados, dados_completos, genes_classificados)
  
  cat("\nAnálise concluída!\n")
  cat("Genes Enriquecido10X:", length(resultados$Enriquecido10X$genes), "\n")
  cat("Genes Enriquecido5X:", length(resultados$Enriquecido5X$genes), "\n")
  cat("Genes Enriquecido2X:", length(resultados$Enriquecido2X$genes), "\n")
  cat("Total de genes classificados:", length(genes_classificados), "\n")
  
  return(list(resultados = resultados, dados_completos = dados_completos))
}

# Função principal para executar a análise
executar_analise <- function() {
  cat("===========================================\n")
  cat("ANÁLISE DE EXPRESSÃO ENRIQUECIDA\n")
  cat("===========================================\n")
  
  # Verificar se os arquivos existem
  if (!file.exists("All_Genes_Healthy_Tissues.csv")) {
    cat("ERRO: Arquivo 'All_Genes_Healthy_Tissues.csv' não encontrado!\n")
    return(NULL)
  }
  
  if (!file.exists("All_Genes_Final_Cancer_Tissue.csv")) {
    cat("ERRO: Arquivo 'All_Genes_Final_Cancer_Tissue.csv' não encontrado!\n")
    return(NULL)
  }
  
  # Executar análise
  resultado_final <- analisar_enriquecimento()
  
  return(resultado_final)
}

# Executar a análise
resultado_final <- executar_analise()

# Função para visualizar a estrutura criada
visualizar_estrutura <- function() {
  cat("\n\nESTRUTURA DE PASTAS CRIADA:\n")
  cat("============================\n")
  
  categorias <- c("Enriquecido10X", "Enriquecido5X", "Enriquecido2X", 
                  "Enriquecido10X_Membrane", "Enriquecido5X_Membrane", "Enriquecido2X_Membrane")
  
  for (categoria in categorias) {
    if (dir.exists(categoria)) {
      cat("\n", categoria, ":\n")
      subdirs <- list.dirs(categoria, recursive = FALSE, full.names = FALSE)
      if (length(subdirs) > 0) {
        for (subdir in subdirs) {
          arquivos <- list.files(file.path(categoria, subdir), pattern = "\\.csv$")
          cat("  ", subdir, "/ (", length(arquivos), " arquivo(s))\n", sep = "")
        }
      }
      arquivos_raiz <- list.files(categoria, pattern = "\\.csv$")
      cat("  Arquivos na raiz:", length(arquivos_raiz), "\n")
    }
  }
}

# Visualizar estrutura
visualizar_estrutura()

cat("\n===========================================\n")
cat("PROCESSO CONCLUÍDO COM SUCESSO!\n")
cat("===========================================\n")