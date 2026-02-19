import os
import pandas as pd
import numpy as np
import re

def clean_tissue_name(tissue_name):
    """
    Limpa o nome do tecido removendo números no final
    """
    # Remove números e espaços no final do nome do tecido
    cleaned = re.sub(r'\s*\d+$', '', tissue_name.strip())
    return cleaned

def process_hpa_expression(input_file, output_dir):
    """
    Processa o arquivo HPA e cria tabelas de expressão por tecido
    """
    
    # Criar pasta de saída
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, "healthy_individual_genes"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "cancer_individual_genes"), exist_ok=True)
    
    # Ler o arquivo
    try:
        df = pd.read_csv(input_file, sep='\t', encoding='utf-8')
    except:
        try:
            df = pd.read_csv(input_file, sep='\t', encoding='latin-1')
        except Exception as e:
            print(f"Erro ao ler arquivo: {e}")
            return
    
    print(f"Total de genes no arquivo: {len(df)}")
    print(f"Colunas disponíveis: {list(df.columns)}")
    
    # Verificar se a coluna Ensembl existe
    ensembl_col = None
    for col in df.columns:
        if 'ensembl' in col.lower():
            ensembl_col = col
            break
    
    if not ensembl_col:
        print("AVISO: Coluna Ensembl não encontrada. Usando a primeira coluna como identificador.")
        ensembl_col = df.columns[0]
    
    # Coletar todos os tecidos únicos das colunas nTPM (saudáveis)
    healthy_tissues = set()
    
    # Coletar todos os tecidos únicos das colunas FPKM (cancer)
    cancer_tissues = set()
    
    # Processar colunas de tecido saudável (nTPM)
    nTPM_columns = [col for col in df.columns if 'nTPM' in col and 'cancer' not in col.lower()]
    print(f"Colunas nTPM saudáveis encontradas: {nTPM_columns}")
    
    # Processar colunas de câncer (FPKM)
    FPKM_columns = [col for col in df.columns if 'FPKM' in col and 'cancer' in col.lower()]
    print(f"Colunas FPKM de câncer encontradas: {FPKM_columns}")
    
    # Coletar todos os tecidos únicos
    for _, row in df.iterrows():
        # Tecidos saudáveis
        for col in nTPM_columns:
            if pd.notna(row[col]) and str(row[col]) != 'NA':
                tissues = str(row[col]).split('\t')
                for tissue in tissues:
                    if ':' in tissue:
                        tissue_name = tissue.split(':')[0].strip()
                        # Limpar nome do tecido
                        clean_name = clean_tissue_name(tissue_name)
                        healthy_tissues.add(clean_name)
        
        # Tecidos cancerígenos
        for col in FPKM_columns:
            if pd.notna(row[col]) and str(row[col]) != 'NA':
                cancers = str(row[col]).split('\t')
                for cancer in cancers:
                    if ':' in cancer:
                        cancer_name = cancer.split(':')[0].strip()
                        # Limpar nome do tecido cancer
                        clean_name = clean_tissue_name(cancer_name)
                        cancer_tissues.add(clean_name)
    
    # Ordenar tecidos
    healthy_tissues = sorted(list(healthy_tissues))
    cancer_tissues = sorted(list(cancer_tissues))
    
    print(f"Total de tecidos saudáveis únicos encontrados: {len(healthy_tissues)}")
    print(f"Total de tecidos cancerígenos únicos encontrados: {len(cancer_tissues)}")
    
    # Criar estruturas para as novas tabelas
    healthy_data = []
    cancer_data = []
    
    # Processar cada gene
    for idx, row in df.iterrows():
        gene_id = str(row[ensembl_col]).strip()
        if not gene_id or gene_id == 'nan':
            continue
            
        # Criar dicionários para este gene
        healthy_expression = {'Ensembl': gene_id}
        cancer_expression = {'Ensembl': gene_id}
        
        # Inicializar todos os tecidos saudáveis com NA
        for tissue in healthy_tissues:
            healthy_expression[tissue] = 'NA'
        
        # Inicializar todos os tecidos cancerígenos com NA
        for tissue in cancer_tissues:
            cancer_expression[tissue] = 'NA'
        
        # Preencher tecidos saudáveis (nTPM)
        for col in nTPM_columns:
            if pd.notna(row[col]) and str(row[col]) != 'NA':
                tissues = str(row[col]).split('\t')
                for tissue_expr in tissues:
                    if ':' in tissue_expr:
                        tissue_name, expr_value = tissue_expr.split(':', 1)
                        tissue_name = clean_tissue_name(tissue_name.strip())
                        expr_value = expr_value.strip()
                        if tissue_name in healthy_tissues:
                            healthy_expression[tissue_name] = expr_value
        
        # Preencher tecidos cancerígenos (FPKM)
        for col in FPKM_columns:
            if pd.notna(row[col]) and str(row[col]) != 'NA':
                cancers = str(row[col]).split('\t')
                for cancer_expr in cancers:
                    if ':' in cancer_expr:
                        cancer_name, expr_value = cancer_expr.split(':', 1)
                        cancer_name = clean_tissue_name(cancer_name.strip())
                        expr_value = expr_value.strip()
                        if cancer_name in cancer_tissues:
                            cancer_expression[cancer_name] = expr_value
        
        healthy_data.append(healthy_expression)
        cancer_data.append(cancer_expression)
    
    # Criar DataFrames finais
    healthy_df = pd.DataFrame(healthy_data)
    cancer_df = pd.DataFrame(cancer_data)
    
    # Reordenar colunas: Ensembl, depois tecidos em ordem alfabética
    healthy_columns = ['Ensembl'] + sorted([col for col in healthy_df.columns if col != 'Ensembl'])
    cancer_columns = ['Ensembl'] + sorted([col for col in cancer_df.columns if col != 'Ensembl'])
    
    healthy_df = healthy_df[healthy_columns]
    cancer_df = cancer_df[cancer_columns]
    
    # Salvar tabelas gerais
    healthy_output_file = os.path.join(output_dir, "Healthy_Genes_Expression.tsv")
    cancer_output_file = os.path.join(output_dir, "Cancer_Genes_Expression.tsv")
    
    healthy_df.to_csv(healthy_output_file, sep='\t', index=False)
    cancer_df.to_csv(cancer_output_file, sep='\t', index=False)
    
    print(f"Tabela saudável salva: {healthy_output_file} ({len(healthy_df)} genes)")
    print(f"Tabela cancer salva: {cancer_output_file} ({len(cancer_df)} genes)")
    
    # Salvar tabelas individuais para cada gene (formato correto)
    healthy_individual_count = 0
    cancer_individual_count = 0
    
    for idx, row in healthy_df.iterrows():
        gene_id = row['Ensembl']
        if gene_id and gene_id != 'nan':
            # Criar DataFrame com apenas este gene (formato wide)
            # Remover a coluna Ensembl e usar apenas os valores de expressão
            expression_values = {tissue: row[tissue] for tissue in healthy_columns[1:]}
            gene_df = pd.DataFrame([expression_values])
            
            # Adicionar linha de cabeçalho com "Expression Value (nTPM)"
            gene_df.index = ["Expression Value (nTPM)"]
            
            # Salvar arquivo individual
            gene_file = os.path.join(output_dir, "healthy_individual_genes", f"{gene_id}.tsv")
            gene_df.to_csv(gene_file, sep='\t', index=True)
            healthy_individual_count += 1
    
    for idx, row in cancer_df.iterrows():
        gene_id = row['Ensembl']
        if gene_id and gene_id != 'nan':
            # Criar DataFrame com apenas este gene (formato wide)
            # Remover a coluna Ensembl e usar apenas os valores de expressão
            expression_values = {tissue: row[tissue] for tissue in cancer_columns[1:]}
            gene_df = pd.DataFrame([expression_values])
            
            # Adicionar linha de cabeçalho com "Expression Value (FPKM)"
            gene_df.index = ["Expression Value (FPKM)"]
            
            # Salvar arquivo individual
            gene_file = os.path.join(output_dir, "cancer_individual_genes", f"{gene_id}.tsv")
            gene_df.to_csv(gene_file, sep='\t', index=True)
            cancer_individual_count += 1
    
    print(f"Tabelas individuais saudáveis salvas: {healthy_individual_count} arquivos")
    print(f"Tabelas individuais cancer salvas: {cancer_individual_count} arquivos")
    
    # Gerar relatório
    print("\n=== RELATÓRIO FINAL ===")
    print(f"Total de genes processados: {len(healthy_df)}")
    print(f"Total de tecidos saudáveis: {len(healthy_tissues)}")
    print(f"Total de tecidos cancerígenos: {len(cancer_tissues)}")
    
    # Mostrar alguns tecidos como exemplo
    print(f"\nExemplos de tecidos saudáveis encontrados:")
    for i, tissue in enumerate(healthy_tissues[:10]):
        print(f"  {i+1}. {tissue}")
    if len(healthy_tissues) > 10:
        print(f"  ... e mais {len(healthy_tissues) - 10} tecidos")
    
    print(f"\nExemplos de tecidos cancerígenos encontrados:")
    for i, tissue in enumerate(cancer_tissues[:10]):
        print(f"  {i+1}. {tissue}")
    if len(cancer_tissues) > 10:
        print(f"  ... e mais {len(cancer_tissues) - 10} tecidos")

def main():
    # Configurações
    input_file = "/home/single1/Doutorado_Bernardo/3.4_Identificação_proteínas_expressão_enriquecida/Human_Expression_Atlas/Genes_Prostate/Shared_Prostate.tsv"
    output_dir = "/home/single1/Doutorado_Bernardo/3.4_Identificação_proteínas_expressão_enriquecida/Human_Expression_Atlas/Genes_Prostate/TissueExpression_SharedProstate"
    
    # Processar dados
    process_hpa_expression(input_file, output_dir)
    
    print(f"\nProcessamento concluído! Arquivos salvos em: {output_dir}/")

if __name__ == "__main__":
    main()