import pandas as pd
import os

def filter_hpa_by_ensembl_ids(hpa_file, ensembl_ids_file, output_base_dir):
    """
    Filtra o arquivo HPA baseado em uma lista de IDs do ENSEMBL e salva em CSV + diretórios individuais + arquivos consolidados
    """
    
    # Criar diretório base
    os.makedirs(output_base_dir, exist_ok=True)
    
    # Ler a lista de IDs do ENSEMBL
    try:
        with open(ensembl_ids_file, 'r') as f:
            target_ensembl_ids = {line.strip() for line in f if line.strip()}
        print(f"IDs do ENSEMBL carregados do arquivo: {len(target_ensembl_ids)}")
    except Exception as e:
        print(f"Erro ao ler arquivo de IDs do ENSEMBL: {e}")
        return
    
    # Ler arquivo HPA
    try:
        df_hpa = pd.read_csv(hpa_file, sep='\t', encoding='utf-8')
    except:
        try:
            df_hpa = pd.read_csv(hpa_file, sep='\t', encoding='latin-1')
        except Exception as e:
            print(f"Erro ao ler arquivo HPA: {e}")
            return
    
    print(f"Total de genes no arquivo HPA: {len(df_hpa)}")
    print(f"Colunas disponíveis: {list(df_hpa.columns)}")
    
    # Encontrar coluna Ensembl
    ensembl_col = None
    for col in df_hpa.columns:
        if 'ensembl' in col.lower():
            ensembl_col = col
            break
    
    if not ensembl_col:
        print("ERRO: Coluna Ensembl não encontrada no arquivo HPA")
        print("Colunas disponíveis:", df_hpa.columns.tolist())
        return
    
    print(f"Coluna de identificador usada: {ensembl_col}")
    
    # Filtrar genes
    df_filtered = df_hpa[df_hpa[ensembl_col].isin(target_ensembl_ids)]
    
    print(f"Genes encontrados no HPA: {len(df_filtered)}")
    
    if len(df_filtered) == 0:
        print("AVISO: Nenhum gene da lista foi encontrado no arquivo HPA")
        print("Verifique se os IDs do ENSEMBL no arquivo correspondem à coluna usada no HPA")
        return
    
    # Salvar arquivo CSV consolidado com todos os genes
    output_csv = os.path.join(output_base_dir, "Selected_Genes_HPA.csv")
    df_filtered.to_csv(output_csv, index=False)
    print(f"Arquivo consolidado salvo: {output_csv}")
    
    # Criar diretórios para categorias
    categories = {
        'All_Genes': 'Todos os genes',
        'Single_Prostate': 'Expressão exclusiva em próstata',
        'Genes_Prostate': 'Expressão compartilhada em próstata', 
        'Single_ProstateCancer': 'Expressão exclusiva em adenocarcinoma',
        'Genes_ProstateCancer': 'Expressão compartilhada em adenocarcinoma'
    }
    
    for category in categories.keys():
        os.makedirs(os.path.join(output_base_dir, category), exist_ok=True)
    
    # Listas para armazenar genes de cada categoria (para arquivos consolidados)
    single_prostate_genes = []
    genes_prostate_genes = []
    single_prostate_cancer_genes = []
    genes_prostate_cancer_genes = []
    
    # Contadores para classificação
    counters = {category: 0 for category in categories.keys()}
    
    # Processar cada gene e classificar
    for idx, row in df_filtered.iterrows():
        ensembl_id = str(row[ensembl_col]).strip()
        
        # Salvar em All_Genes (todos os genes)
        save_individual_gene(row, output_base_dir, 'All_Genes', ensembl_id)
        counters['All_Genes'] += 1
        
        # Obter distribuições
        tissue_distribution = str(row.get('RNA tissue distribution', '')).strip()
        cancer_distribution = str(row.get('RNA cancer distribution', '')).strip()
        
        # Verificar expressão prostática
        prostate_expr_found = check_prostate_expression(row)
        cancer_expr_found = check_cancer_expression(row)
        
        # Classificar nas categorias específicas
        if tissue_distribution == 'Detected in single' and prostate_expr_found:
            save_individual_gene(row, output_base_dir, 'Single_Prostate', ensembl_id)
            single_prostate_genes.append(row)
            counters['Single_Prostate'] += 1
        
        elif prostate_expr_found:
            save_individual_gene(row, output_base_dir, 'Genes_Prostate', ensembl_id)
            genes_prostate_genes.append(row)
            counters['Genes_Prostate'] += 1
        
        if cancer_distribution == 'Detected in single' and cancer_expr_found:
            save_individual_gene(row, output_base_dir, 'Single_ProstateCancer', ensembl_id)
            single_prostate_cancer_genes.append(row)
            counters['Single_ProstateCancer'] += 1
        
        elif cancer_expr_found:
            save_individual_gene(row, output_base_dir, 'Genes_ProstateCancer', ensembl_id)
            genes_prostate_cancer_genes.append(row)
            counters['Genes_ProstateCancer'] += 1
    
    # Criar arquivos consolidados para cada categoria
    create_consolidated_files(
        output_base_dir,
        single_prostate_genes,
        genes_prostate_genes,
        single_prostate_cancer_genes,
        genes_prostate_cancer_genes
    )
    
    # Gerar relatório
    print("\n=== ESTATÍSTICAS DETALHADAS ===")
    print(f"IDs do ENSEMBL fornecidos: {len(target_ensembl_ids)}")
    print(f"Genes encontrados no HPA: {len(df_filtered)}")
    print(f"Genes não encontrados: {len(target_ensembl_ids) - len(df_filtered)}")
    
    print("\n=== CLASSIFICAÇÃO POR EXPRESSÃO ===")
    for category, description in categories.items():
        print(f"{description}: {counters[category]}")
    
    # Mostrar primeiros genes encontrados
    print(f"\nPrimeiros 10 genes encontrados:")
    found_ids = df_filtered[ensembl_col].head(10).tolist()
    for i, ensembl_id in enumerate(found_ids, 1):
        print(f"  {i}. {ensembl_id}")
    
    return df_filtered, counters

def check_prostate_expression(row):
    """Verifica se há expressão prostática em qualquer coluna nTPM"""
    for i in range(1, 6):
        nTPM_col = f'RNA tissue specific nTPM ({i})'
        if nTPM_col in row and pd.notna(row[nTPM_col]):
            value = str(row[nTPM_col])
            if 'prostate:' in value.lower():
                return True
    return False

def check_cancer_expression(row):
    """Verifica se há expressão de câncer em qualquer coluna FPKM"""
    for i in range(1, 6):
        FPKM_col = f'RNA cancer specific FPKM ({i})'
        if FPKM_col in row and pd.notna(row[FPKM_col]):
            value = str(row[FPKM_col])
            if 'Prostate Adenocarcinoma (TCGA):' in value:
                return True
    return False

def save_individual_gene(row, base_dir, category, ensembl_id):
    """Salva um gene individual na categoria especificada"""
    output_dir = os.path.join(base_dir, category)
    filename = f"{ensembl_id}.tsv"
    filepath = os.path.join(output_dir, filename)
    
    # Criar DataFrame com uma linha
    gene_df = pd.DataFrame([row])
    gene_df.to_csv(filepath, sep='\t', index=False)

def create_consolidated_files(output_base_dir, single_prostate, genes_prostate, single_cancer, genes_cancer):
    """Cria arquivos consolidados para cada categoria"""
    
    # Single Prostate -> Unique_Prostate.tsv
    if single_prostate:
        df_single_prostate = pd.DataFrame(single_prostate)
        output_file = os.path.join(output_base_dir, "Unique_Prostate.tsv")
        df_single_prostate.to_csv(output_file, sep='\t', index=False)
        print(f"Arquivo consolidado criado: Unique_Prostate.tsv ({len(single_prostate)} genes)")
    
    # Genes Prostate -> Shared_Prostate.tsv
    if genes_prostate:
        df_genes_prostate = pd.DataFrame(genes_prostate)
        output_file = os.path.join(output_base_dir, "Shared_Prostate.tsv")
        df_genes_prostate.to_csv(output_file, sep='\t', index=False)
        print(f"Arquivo consolidado criado: Shared_Prostate.tsv ({len(genes_prostate)} genes)")
    
    # Single Prostate Cancer -> Unique_Prostate_Cancer.tsv
    if single_cancer:
        df_single_cancer = pd.DataFrame(single_cancer)
        output_file = os.path.join(output_base_dir, "Unique_Prostate_Cancer.tsv")
        df_single_cancer.to_csv(output_file, sep='\t', index=False)
        print(f"Arquivo consolidado criado: Unique_Prostate_Cancer.tsv ({len(single_cancer)} genes)")
    
    # Genes Prostate Cancer -> Shared_Prostate_Cancer.tsv
    if genes_cancer:
        df_genes_cancer = pd.DataFrame(genes_cancer)
        output_file = os.path.join(output_base_dir, "Shared_Prostate_Cancer.tsv")
        df_genes_cancer.to_csv(output_file, sep='\t', index=False)
        print(f"Arquivo consolidado criado: Shared_Prostate_Cancer.tsv ({len(genes_cancer)} genes)")

def analyze_expression_patterns(df, ensembl_col, output_base_dir):
    """
    Analisa padrões de expressão prostática nos genes filtrados
    """
    
    # Contadores para análise
    counters = {
        'total_genes': len(df),
        'prostate_expr_any': 0,
        'prostate_expr_single': 0,
        'cancer_expr_any': 0,
        'cancer_expr_single': 0
    }
    
    # Verificar padrões de expressão
    for idx, row in df.iterrows():
        tissue_distribution = str(row.get('RNA tissue distribution', '')).strip()
        cancer_distribution = str(row.get('RNA cancer distribution', '')).strip()
        
        # Verificar expressão prostática
        prostate_found = check_prostate_expression(row)
        cancer_found = check_cancer_expression(row)
        
        # Contar padrões
        if prostate_found:
            counters['prostate_expr_any'] += 1
            if tissue_distribution == 'Detected in single':
                counters['prostate_expr_single'] += 1
        
        if cancer_found:
            counters['cancer_expr_any'] += 1
            if cancer_distribution == 'Detected in single':
                counters['cancer_expr_single'] += 1
    
    # Salvar análise
    analysis_data = {
        'Categoria': [
            'Total de genes',
            'Expressão em próstata (qualquer)',
            'Expressão exclusiva em próstata',
            'Expressão em adenocarcinoma (qualquer)',
            'Expressão exclusiva em adenocarcinoma'
        ],
        'Quantidade': [
            counters['total_genes'],
            counters['prostate_expr_any'],
            counters['prostate_expr_single'],
            counters['cancer_expr_any'],
            counters['cancer_expr_single']
        ]
    }
    
    df_analysis = pd.DataFrame(analysis_data)
    analysis_csv = os.path.join(output_base_dir, "Expression_Analysis_Selected_Genes.csv")
    df_analysis.to_csv(analysis_csv, index=False)
    
    print(f"Análise de expressão salva: {analysis_csv}")
    
    # Mostrar resultados
    print("\n=== ANÁLISE DETALHADA DE EXPRESSÃO ===")
    for i, row in df_analysis.iterrows():
        print(f"{row['Categoria']}: {row['Quantidade']}")
    
    return df_analysis

def main():
    # Configurações
    hpa_file = "HPA_Clean.tsv"      # Arquivo completo do HPA
    ensembl_ids_file = "All_IDs_ENSG.txt"    # Arquivo com lista de IDs do ENSEMBL
    output_base_dir = "Results_All_Ensembl"
    
    print("=== FILTRO E CLASSIFICAÇÃO DE GENES DO HPA ===")
    print("Usando IDs do ENSEMBL para filtragem")
    
    # Filtrar e classificar genes
    df_selected, classification_counters = filter_hpa_by_ensembl_ids(hpa_file, ensembl_ids_file, output_base_dir)
    
    if df_selected is not None and len(df_selected) > 0:
        # Encontrar coluna Ensembl novamente para análise
        ensembl_col = None
        for col in df_selected.columns:
            if 'ensembl' in col.lower():
                ensembl_col = col
                break
        
        # Analisar padrões de expressão
        analyze_expression_patterns(df_selected, ensembl_col, output_base_dir)
        
        print(f"\n✅ Processamento concluído!")
        print(f"📁 Diretório base: {output_base_dir}")
        print(f"\n📊 ARQUIVOS CONSOLIDADOS:")
        print(f"   - Selected_Genes_HPA.csv (todos os genes)")
        print(f"   - Unique_Prostate.tsv (expressão exclusiva em próstata)")
        print(f"   - Shared_Prostate.tsv (expressão compartilhada em próstata)")
        print(f"   - Unique_Prostate_Cancer.tsv (expressão exclusiva em adenocarcinoma)")
        print(f"   - Shared_Prostate_Cancer.tsv (expressão compartilhada em adenocarcinoma)")
        print(f"   - Expression_Analysis_Selected_Genes.csv (análise de expressão)")
        
        print(f"\n📂 DIRETÓRIOS COM GENES INDIVIDUAIS:")
        print(f"   - All_Genes/ (todos os genes)")
        print(f"   - Single_Prostate/ (expressão exclusiva em próstata)")
        print(f"   - Genes_Prostate/ (expressão compartilhada em próstata)")
        print(f"   - Single_ProstateCancer/ (expressão exclusiva em adenocarcinoma)")
        print(f"   - Genes_ProstateCancer/ (expressão compartilhada em adenocarcinoma)")
    else:
        print("❌ Nenhum dado foi processado devido a erros.")

if __name__ == "__main__":
    main()