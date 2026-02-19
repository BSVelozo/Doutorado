#!/bin/bash

# Configurações
DELAY=1
INPUT_FILE="$1"
OUTPUT_FILE="uniprot_ensembl_found.tsv"
HTML_DIR="Gene_HTML"

# Verifica se o EDirect está instalado
if ! command -v esearch &> /dev/null || ! command -v efetch &> /dev/null; then
    echo "Erro: EDirect não encontrado. Instale de: https://www.ncbi.nlm.nih.gov/books/NBK179288/"
    exit 1
fi

# Verifica argumentos
if [ $# -ne 1 ]; then
    echo "Uso: $0 <arquivo_de_ids.txt>"
    exit 1
fi

if [ ! -f "$INPUT_FILE" ]; then
    echo "Erro: Arquivo '$INPUT_FILE' não encontrado!"
    exit 1
fi

# Cria diretório para armazenar páginas HTML
mkdir -p "$HTML_DIR"

# Inicializa arquivo de saída
echo -e "UNIPROT\tENSEMBL(Gene)\tENSEMBL(Transcrito)\tENSEMBL(Proteína)\tOrigem" > "$OUTPUT_FILE"

# Lê IDs
mapfile -t UNIPROT_IDS < "$INPUT_FILE"
TOTAL_IDS=${#UNIPROT_IDS[@]}
PROCESSED=0
FOUND_COUNT=0

echo "Iniciando processamento de $TOTAL_IDS IDs..."
echo "Páginas serão salvas APENAS para IDs com correspondência em: $HTML_DIR/"
START_TIME=$(date +%s)

# Processa cada ID
for UNIPROT_ID in "${UNIPROT_IDS[@]}"; do
    UNIPROT_ID=$(echo "$UNIPROT_ID" | tr -d '[:space:]' | tr -d '\r')
    
    if [ -z "$UNIPROT_ID" ]; then
        continue
    fi

    PROCESSED=$((PROCESSED + 1))

    # Arquivo temporário para processamento
    TEMP_FILE=$(mktemp)
    
    # Consulta NCBI
    echo "Consultando: $UNIPROT_ID"
    esearch -db protein -query "$UNIPROT_ID" 2>/dev/null | \
    efetch -format xrefs 2>/dev/null > "$TEMP_FILE"
    
    # Extrai os diferentes tipos de IDs do Ensembl
    ENSEMBL_GENE=""
    ENSEMBL_TRANSCRIPT=""
    ENSEMBL_PROTEIN=""
    
    if [ -s "$TEMP_FILE" ]; then
        # Extrai gene (ENSG)
        ENSEMBL_GENE=$(grep -o "Ensembl:ENSG[^,]*" "$TEMP_FILE" | head -1 | sed 's/Ensembl://; s/\.[0-9]*//' 2>/dev/null)
        
        # Extrai transcrito (ENST)
        ENSEMBL_TRANSCRIPT=$(grep -o "Ensembl:ENST[^,]*" "$TEMP_FILE" | head -1 | sed 's/Ensembl://; s/\.[0-9]*//' 2>/dev/null)
        
        # Extrai proteína (ENSP)
        ENSEMBL_PROTEIN=$(grep -o "Ensembl:ENSP[^,]*" "$TEMP_FILE" | head -1 | sed 's/Ensembl://; s/\.[0-9]*//' 2>/dev/null)
    fi

    # Se encontrou pelo menos um ID do Ensembl, salva a página e escreve no arquivo de correlação
    if [ -n "$ENSEMBL_GENE" ] || [ -n "$ENSEMBL_TRANSCRIPT" ] || [ -n "$ENSEMBL_PROTEIN" ]; then
        # Salva a página no diretório Gene_HTML
        HTML_FILE="$HTML_DIR/${UNIPROT_ID}.txt"
        cp "$TEMP_FILE" "$HTML_FILE"
        
        # Substitui valores vazios por "-"
        [ -z "$ENSEMBL_GENE" ] && ENSEMBL_GENE="-"
        [ -z "$ENSEMBL_TRANSCRIPT" ] && ENSEMBL_TRANSCRIPT="-"
        [ -z "$ENSEMBL_PROTEIN" ] && ENSEMBL_PROTEIN="-"
        
        # Escreve no arquivo de correlação
        echo -e "$UNIPROT_ID\t$ENSEMBL_GENE\t$ENSEMBL_TRANSCRIPT\t$ENSEMBL_PROTEIN\tNCBI" >> "$OUTPUT_FILE"
        FOUND_COUNT=$((FOUND_COUNT + 1))
        
        # Mostra o que foi encontrado
        FOUND_ITEMS=""
        [ "$ENSEMBL_GENE" != "-" ] && FOUND_ITEMS="$FOUND_ITEMS Gene:$ENSEMBL_GENE"
        [ "$ENSEMBL_TRANSCRIPT" != "-" ] && FOUND_ITEMS="$FOUND_ITEMS Transcrito:$ENSEMBL_TRANSCRIPT"
        [ "$ENSEMBL_PROTEIN" != "-" ] && FOUND_ITEMS="$FOUND_ITEMS Proteína:$ENSEMBL_PROTEIN"
        echo "  ✓ Encontrado: $UNIPROT_ID ->$FOUND_ITEMS (página salva)"
    else
        echo "  ✗ Não encontrado: $UNIPROT_ID (página não salva)"
    fi

    # Remove arquivo temporário
    rm -f "$TEMP_FILE"

    # Delay
    sleep "$DELAY"

    # Progresso a cada 10 IDs
    if [ $((PROCESSED % 10)) -eq 0 ]; then
        PERCENT=$((PROCESSED * 100 / TOTAL_IDS))
        CURRENT_TIME=$(date +%s)
        ELAPSED=$((CURRENT_TIME - START_TIME))
        if [ $PROCESSED -gt 0 ]; then
            REMAINING=$(( (TOTAL_IDS - PROCESSED) * ELAPSED / PROCESSED ))
            echo "=== Progresso: $PERCENT% | Processados: $PROCESSED/$TOTAL_IDS | Encontrados: $FOUND_COUNT | Restante: ~${REMAINING}s ==="
        else
            echo "=== Progresso: $PERCENT% | Processados: $PROCESSED/$TOTAL_IDS | Encontrados: $FOUND_COUNT ==="
        fi
    fi
done

# Estatísticas finais
END_TIME=$(date +%s)
TOTAL_TIME=$((END_TIME - START_TIME))
NOT_FOUND_COUNT=$((PROCESSED - FOUND_COUNT))

echo ""
echo "================================================"
echo "         ESTATÍSTICAS FINAIS"
echo "================================================"
echo "Arquivo de entrada: $INPUT_FILE"
echo "Total de IDs processadas: $PROCESSED"
echo "IDs com correspondência encontrada: $FOUND_COUNT"
echo "IDs sem correspondência: $NOT_FOUND_COUNT"
echo "Taxa de sucesso: $((FOUND_COUNT * 100 / PROCESSED))%"
echo "Tempo total de processamento: ${TOTAL_TIME}s"
echo ""
echo "Arquivos gerados:"
echo "- $OUTPUT_FILE ($FOUND_COUNT correlações encontradas)"
echo "- $HTML_DIR/ ($FOUND_COUNT páginas salvas)"
echo ""
echo "Concluído!"