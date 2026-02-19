import pandas as pd
from datetime import datetime
import sys

def remove_version(ensembl_id):
    """Remove a versão dos IDs ENSEMBL"""
    if pd.isna(ensembl_id):
        return ensembl_id
    # Remove tudo após o ponto, se houver
    return str(ensembl_id).split('.')[0]

def is_valid_ensembl_id(ensembl_id, prefix):
    """Verifica se o ID ENSEMBL é válido (começa com o prefixo correto)"""
    if pd.isna(ensembl_id):
        return False
    return str(ensembl_id).startswith(prefix)

def processar_arquivos():
    """Processa os arquivos e gera o arquivo combinado e estatísticas"""
    
    print("Lendo arquivos...")
    
    # Carregar os arquivos com os nomes corretos
    try:
        # Ler apenas as duas primeiras colunas de cada arquivo
        gene_df = pd.read_csv('UniProt_Ensembl_Genome.tsv', sep='\t', usecols=[0, 1], dtype=str)
        transcript_df = pd.read_csv('UniProt_Ensembl_Transcrito.tsv', sep='\t', usecols=[0, 1], dtype=str)
        protein_df = pd.read_csv('UniProt_Ensembl_Protein.tsv', sep='\t', usecols=[0, 1], dtype=str)
    except FileNotFoundError as e:
        print(f"Erro: Arquivo não encontrado - {e}")
        sys.exit(1)
    
    # Renomear colunas para consistência
    gene_df.columns = ['UNIPROT', 'ENSEMBL(Gene)']
    transcript_df.columns = ['UNIPROT', 'ENSEMBL(Transcrito)']
    protein_df.columns = ['UNIPROT', 'ENSEMBL(Proteina)']
    
    print(f"Total de linhas no arquivo de gene: {len(gene_df)}")
    print(f"Total de linhas no arquivo de transcrito: {len(transcript_df)}")
    print(f"Total de linhas no arquivo de proteína: {len(protein_df)}")
    
    # Remover duplicatas mantendo apenas o primeiro registro
    # Para transcritos e proteínas, alguns UNIPROT têm múltiplas entradas
    transcript_df_sem_dups = transcript_df.drop_duplicates(subset='UNIPROT', keep='first')
    protein_df_sem_dups = protein_df.drop_duplicates(subset='UNIPROT', keep='first')
    
    print(f"Transcritos únicos: {len(transcript_df_sem_dups)}")
    print(f"Proteínas únicas: {len(protein_df_sem_dups)}")
    
    # Combinar os dataframes
    merged_df = gene_df.merge(
        transcript_df_sem_dups, 
        on='UNIPROT', 
        how='left'
    ).merge(
        protein_df_sem_dups, 
        on='UNIPROT', 
        how='left'
    )
    
    # Aplicar função para remover versões
    for col in ['ENSEMBL(Gene)', 'ENSEMBL(Transcrito)', 'ENSEMBL(Proteina)']:
        merged_df[col] = merged_df[col].apply(remove_version)
    
    # Adicionar coluna de origem
    merged_df['Origem'] = 'UNIPROT'
    
    # Filtrar apenas os IDs válidos
    print("\nFiltrando IDs válidos...")
    
    # Contar antes da filtragem
    total_antes = len(merged_df)
    
    # Filtrar linhas onde todos os IDs começam com os prefixos corretos
    # Gene deve começar com ENSG
    mask_gene = merged_df['ENSEMBL(Gene)'].apply(lambda x: is_valid_ensembl_id(x, 'ENSG'))
    
    # Transcrito deve começar com ENST (mas pode ser NaN)
    mask_transcript = merged_df['ENSEMBL(Transcrito)'].apply(
        lambda x: pd.isna(x) or is_valid_ensembl_id(x, 'ENST')
    )
    
    # Proteína deve começar com ENSP (mas pode ser NaN)
    mask_protein = merged_df['ENSEMBL(Proteina)'].apply(
        lambda x: pd.isna(x) or is_valid_ensembl_id(x, 'ENSP')
    )
    
    # Aplicar todas as máscaras
    filtered_df = merged_df[mask_gene & mask_transcript & mask_protein].copy()
    
    # Contar após a filtragem
    total_depois = len(filtered_df)
    removidos = total_antes - total_depois
    
    print(f"Total antes da filtragem: {total_antes}")
    print(f"Total após filtragem: {total_depois}")
    print(f"Linhas removidas: {removidos}")
    
    # Ordenar por UNIPROT
    filtered_df = filtered_df.sort_values('UNIPROT')
    
    # Calcular estatísticas APÓS a filtragem
    total_uniprots = len(filtered_df)
    
    # Contar UNIPROTs com transcrito (válido)
    uniprots_com_transcrito = filtered_df['ENSEMBL(Transcrito)'].notna().sum()
    
    # Contar UNIPROTs com proteína (válida)
    uniprots_com_proteina = filtered_df['ENSEMBL(Proteina)'].notna().sum()
    
    # Contar UNIPROTs com correspondência completa
    uniprots_completo = filtered_df.dropna(subset=['ENSEMBL(Transcrito)', 'ENSEMBL(Proteina)']).shape[0]
    
    # UNIPROTs sem correspondência completa
    uniprots_incompletos = total_uniprots - uniprots_completo
    
    # Salvar arquivo final
    filtered_df.to_csv('3_UniProt_ensembl_found_UNIPROT.tsv', sep='\t', index=False)
    
    # Gerar estatísticas
    gerar_estatisticas(
        total_uniprots=total_uniprots,
        uniprots_com_transcrito=uniprots_com_transcrito,
        uniprots_com_proteina=uniprots_com_proteina,
        uniprots_completo=uniprots_completo,
        uniprots_incompletos=uniprots_incompletos,
        removidos=removidos
    )
    
    print(f"\nProcessamento concluído!")
    print(f"Arquivo gerado: '3_UniProt_ensembl_found_UNIPROT.tsv'")
    print(f"Estatísticas geradas: 'estatisticas_correspondencia.txt'")
    
    return filtered_df

def gerar_estatisticas(total_uniprots, uniprots_com_transcrito, 
                       uniprots_com_proteina, uniprots_completo, 
                       uniprots_incompletos, removidos):
    """Gera arquivo de estatísticas"""
    
    data_processamento = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Calcular porcentagens
    taxa_gene = 100.0  # Todos têm gene válido após filtragem
    taxa_transcrito = (uniprots_com_transcrito / total_uniprots * 100) if total_uniprots > 0 else 0
    taxa_proteina = (uniprots_com_proteina / total_uniprots * 100) if total_uniprots > 0 else 0
    taxa_completa = (uniprots_completo / total_uniprots * 100) if total_uniprots > 0 else 0
    
    # Formatar estatísticas
    estatisticas = f"""ESTATÍSTICAS DE CORRESPONDÊNCIA UNIPROT->ENSEMBL
========================================================
Data de processamento: {data_processamento}
Arquivos de entrada: 
  - UniProt_Ensembl_Genome.tsv
  - UniProt_Ensembl_Transcrito.tsv  
  - UniProt_Ensembl_Protein.tsv
Total de IDs UNIPROT processados: {total_uniprots}
IDs com correspondência completa (gene + transcrito + proteína): {uniprots_completo}
IDs com correspondência parcial: {uniprots_incompletos}
Taxa de sucesso "Genoma": {taxa_gene:.1f} %
Taxa de sucesso "Transcriptoma": {taxa_transcrito:.1f} %
Taxa de sucesso "Proteínas": {taxa_proteina:.1f} %
Taxa de correspondência completa: {taxa_completa:.1f} %
IDs removidos (não ENSG/ENST/ENSP): {removidos}
CONFIGURAÇÃO:
- Método: Mapeamento UniProt para ENSEMBL
- Versões removidas: Sim
- Filtro aplicado: Apenas ENSG/ENST/ENSP
- Arquivo de saída: 3_UniProt_ensembl_found_UNIPROT.tsv
- Origem: UNIPROT

DETALHAMENTO:
- UNIPROTs processados (após filtro): {total_uniprots}
- UNIPROTs com gene válido (ENSG): {total_uniprots} (100%)
- UNIPROTs com transcrito válido (ENST): {uniprots_com_transcrito} ({taxa_transcrito:.1f}%)
- UNIPROTs com proteína válida (ENSP): {uniprots_com_proteina} ({taxa_proteina:.1f}%)
- UNIPROTs com todos os três: {uniprots_completo} ({taxa_completa:.1f}%)
- UNIPROTs sem transcrito válido: {total_uniprots - uniprots_com_transcrito}
- UNIPROTs sem proteína válida: {total_uniprots - uniprots_com_proteina}
"""
    
    # Salvar estatísticas em arquivo
    with open('estatisticas_correspondencia.txt', 'w', encoding='utf-8') as f:
        f.write(estatisticas)
    
    # Imprimir resumo no console
    print("\n" + "="*60)
    print("RESUMO DAS ESTATÍSTICAS:")
    print("="*60)
    print(f"Total de UNIPROTs após filtro: {total_uniprots}")
    print(f"IDs removidos (não ENSG/ENST/ENSP): {removidos}")
    print(f"Com transcrito (ENST): {uniprots_com_transcrito} ({taxa_transcrito:.1f}%)")
    print(f"Com proteína (ENSP): {uniprots_com_proteina} ({taxa_proteina:.1f}%)")
    print(f"Correspondência completa: {uniprots_completo} ({taxa_completa:.1f}%)")
    print("="*60)

def main():
    """Função principal"""
    print("Processando arquivos UniProt-ENSEMBL...")
    print("-" * 40)
    
    # Processar arquivos
    df_resultado = processar_arquivos()
    
    # Mostrar preview do resultado
    print("\n" + "="*60)
    print("PREVIEW DO ARQUIVO GERADO (primeiras 10 linhas):")
    print("="*60)
    print(df_resultado.head(10).to_string(index=False))
    
    # Mostrar informações adicionais
    print(f"\nTotal de linhas no arquivo final: {len(df_resultado)}")
    print(f"Colunas: {', '.join(df_resultado.columns)}")
    
    # Mostrar informações sobre valores ausentes
    print(f"\nANÁLISE DE VALORES AUSENTES:")
    print("-" * 30)
    for col in ['ENSEMBL(Gene)', 'ENSEMBL(Transcrito)', 'ENSEMBL(Proteina)']:
        missing = df_resultado[col].isna().sum()
        percent_missing = missing/len(df_resultado)*100 if len(df_resultado) > 0 else 0
        print(f"{col}: {missing} ausentes ({percent_missing:.1f}%)")
    
    # Mostrar exemplos de IDs válidos
    print(f"\nEXEMPLOS DE IDs VÁLIDOS:")
    print("-" * 30)
    print(f"Genes (ENSG): {df_resultado['ENSEMBL(Gene)'].iloc[0] if len(df_resultado) > 0 else 'N/A'}")
    print(f"Transcritos (ENST): {df_resultado['ENSEMBL(Transcrito)'].iloc[0] if len(df_resultado) > 0 else 'N/A'}")
    print(f"Proteínas (ENSP): {df_resultado['ENSEMBL(Proteina)'].iloc[0] if len(df_resultado) > 0 else 'N/A'}")

if __name__ == "__main__":
    main()