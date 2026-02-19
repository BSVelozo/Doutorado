import pandas as pd

# Configurações
UNIPROT_IDS_FILE = "3_IDs_uniprot.txt"
HPA_DATA_FILE = "All_Data_HPA_Brute.tsv"
OUTPUT_MATCHES = "UniProt_ensembl_found_HPA.tsv"
OUTPUT_STATS = "log.txt"
ORIGEM = "Human Protein Atlas"

# Lê os IDs UniProt do arquivo de entrada
with open(UNIPROT_IDS_FILE, 'r') as f:
    uniprot_ids = [line.strip() for line in f if line.strip()]

total_ids = len(uniprot_ids)

# Lê o arquivo HPA
hpa_df = pd.read_csv(HPA_DATA_FILE, sep='\t', usecols=['Uniprot', 'Ensembl'])
# Converte para dicionário para busca rápida
uniprot_to_ensembl_dict = hpa_df.dropna(subset=['Uniprot', 'Ensembl']).set_index('Uniprot')['Ensembl'].to_dict()

matches = []
not_found_count = 0

for uniprot_id in uniprot_ids:
    ensembl_gene = uniprot_to_ensembl_dict.get(uniprot_id)
    
    if ensembl_gene:
        matches.append({
            'UNIPROT': uniprot_id,
            'ENSEMBL(Gene)': ensembl_gene,
            'Origem': ORIGEM
        })
    else:
        not_found_count += 1

# Cria DataFrame com os matches
matches_df = pd.DataFrame(matches)

# Salva arquivo com correspondências
matches_df.to_csv(OUTPUT_MATCHES, sep='\t', index=False)

# Calcula estatísticas
found_count = len(matches)
success_rate = (found_count / total_ids) * 100

# Gera relatório de estatísticas
with open(OUTPUT_STATS, 'w') as f:
    f.write(f"Total de IDs processadas: {total_ids}\n")
    f.write(f"IDs com correspondência encontrada: {found_count}\n")
    f.write(f"IDs sem correspondência: {not_found_count}\n")
    f.write(f"Taxa de sucesso: {success_rate:.0f}%\n")

print("Processamento concluído!")
print(f"Arquivo com correspondências salvo como: {OUTPUT_MATCHES}")
print(f"Estatísticas salvas como: {OUTPUT_STATS}")