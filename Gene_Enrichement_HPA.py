import pandas as pd
import numpy as np
import os
import re

def extract_expression_values(row, nTPM_cols, FPKM_cols):
    """
    Extrai valores de expressão das colunas nTPM e FPKM
    Retorna: dicionário {tecido: valor}
    """
    expressions = {}
    
    # Processar colunas nTPM
    for col in nTPM_cols:
        if col in row and pd.notna(row[col]) and str(row[col]) != 'NA':
            parts = str(row[col]).strip().split(':')
            if len(parts) == 2:
                tissue = parts[0].strip()
                try:
                    value = float(parts[1].strip())
                    # Se o tecido já existe, manter o maior valor
                    if tissue not in expressions or value > expressions[tissue]:
                        expressions[tissue] = value
                except ValueError:
                    pass
    
    # Processar colunas FPKM
    for col in FPKM_cols:
        if col in row and pd.notna(row[col]) and str(row[col]) != 'NA':
            parts = str(row[col]).strip().split(':')
            if len(parts) == 2:
                tissue = parts[0].strip()
                try:
                    value = float(parts[1].strip())
                    # Se o tecido já existe, manter o maior valor
                    if tissue not in expressions or value > expressions[tissue]:
                        expressions[tissue] = value
                except ValueError:
                    pass
    
    return expressions

def calculate_enrichment(expressions):
    """
    Calcula o enriquecimento de um gene
    Retorna: (categoria, tecido_enriquecido, ratio, max_expression)
    """
    if len(expressions) < 2:
        return None, None, None, None
    
    # Encontrar o tecido com maior expressão
    tissues = list(expressions.keys())
    values = list(expressions.values())
    
    max_idx = np.argmax(values)
    max_tissue = tissues[max_idx]
    max_value = values[max_idx]
    
    # Calcular média dos outros tecidos
    other_values = [values[i] for i in range(len(values)) if i != max_idx]
    mean_other = np.mean(other_values)
    
    if mean_other == 0:
        return None, None, None, None
    
    ratio = max_value / mean_other
    
    # Classificar (apenas categoria mais alta)
    if ratio >= 10:
        return 'Enriquecido10X', max_tissue, ratio, max_value
    elif ratio >= 5:
        return 'Enriquecido5X', max_tissue, ratio, max_value
    elif ratio >= 2:
        return 'Enriquecido2X', max_tissue, ratio, max_value
    else:
        return None, None, None, None

def process_files(file_list):
    """
    Processa todos os arquivos e extrai dados de expressão
    """
    all_genes = []
    nTPM_cols = [f'RNA tissue specific nTPM ({i})' for i in range(1, 6)]
    FPKM_cols = [f'RNA cancer specific FPKM ({i})' for i in range(1, 6)]
    
    for file_path in file_list:
        filename = os.path.basename(file_path)
        source = filename.replace('.tsv', '')
        
        try:
            df = pd.read_csv(file_path, sep='\t', encoding='utf-8')
        except:
            try:
                df = pd.read_csv(file_path, sep='\t', encoding='latin-1')
            except Exception as e:
                print(f"Erro ao ler {filename}: {e}")
                continue
        
        print(f"Processando {filename}: {len(df)} genes")
        
        for idx, row in df.iterrows():
            ensembl = str(row.get('Ensembl', '')).strip()
            if not ensembl or ensembl == 'nan':
                continue
            
            # Extrair valores de expressão
            expressions = extract_expression_values(row, nTPM_cols, FPKM_cols)
            
            if expressions:
                gene_data = {
                    'Ensembl': ensembl,
                    'Gene': row.get('Gene', ''),
                    'Uniprot': row.get('Uniprot', ''),
                    'Source': source,
                    'Expressions': expressions,
                    'Full_Row': row.to_dict()
                }
                all_genes.append(gene_data)
    
    print(f"\nTotal de genes processados: {len(all_genes)}")
    return all_genes

def classify_genes(all_genes):
    """
    Classifica genes nas categorias de enriquecimento
    """
    categories = {
        'Enriquecido10X': [],
        'Enriquecido5X': [],
        'Enriquecido2X': []
    }
    
    categorized_info = {cat: [] for cat in categories}
    all_enriched_genes = []  # Lista para todos os genes enriquecidos
    
    for gene_data in all_genes:
        category, tissue, ratio, max_expr = calculate_enrichment(gene_data['Expressions'])
        
        if category:
            # Adicionar informações adicionais
            gene_info = {
                'Ensembl': gene_data['Ensembl'],
                'Gene': gene_data['Gene'],
                'Uniprot': gene_data['Uniprot'],
                'Source': gene_data['Source'],
                'Tissue': tissue,
                'Ratio': ratio,
                'Max_Expression': max_expr,
                'Num_Tissues': len(gene_data['Expressions']),
                'Category': category,  # Adicionado para o arquivo consolidado
                'Expressions': gene_data['Expressions']
            }
            
            categories[category].append(gene_info)
            categorized_info[category].append(gene_info)
            all_enriched_genes.append(gene_info)
    
    # Estatísticas
    print("\n=== CLASSIFICAÇÃO DE ENRIQUECIMENTO ===")
    for category in categories:
        print(f"{category}: {len(categories[category])} genes")
    
    print(f"Total de genes enriquecidos: {len(all_enriched_genes)}")
    
    return categories, categorized_info, all_enriched_genes

def load_filter_files(membrane_file, interactors_file):
    """
    Carrega os arquivos de filtro
    """
    membrane_genes = set()
    interactors_genes = set()
    
    # Carregar genes de membrana
    if os.path.exists(membrane_file):
        with open(membrane_file, 'r') as f:
            for line in f:
                ensembl = line.strip()
                if ensembl:
                    membrane_genes.add(ensembl)
        print(f"Genes de membrana carregados: {len(membrane_genes)}")
    
    # Carregar genes interatores
    if os.path.exists(interactors_file):
        with open(interactors_file, 'r') as f:
            for line in f:
                ensembl = line.strip()
                if ensembl:
                    interactors_genes.add(ensembl)
        print(f"Genes interatores carregados: {len(interactors_genes)}")
    
    return membrane_genes, interactors_genes

def save_all_proteins_results(categories, all_enriched_genes, output_base_dir):
    """
    Salva todos os resultados gerais na pasta All_Proteins
    """
    all_proteins_dir = os.path.join(output_base_dir, "All_Proteins")
    os.makedirs(all_proteins_dir, exist_ok=True)
    
    # 1. Salvar arquivo consolidado com todas as proteínas enriquecidas
    if all_enriched_genes:
        # Criar DataFrame com as colunas solicitadas
        consolidated_data = []
        for gene_info in all_enriched_genes:
            consolidated_data.append({
                'ID_Ensembl': gene_info['Ensembl'],
                'Tecido_enriquecido': gene_info['Tissue'],
                'Categoria': gene_info['Category']
            })
        
        df_consolidated = pd.DataFrame(consolidated_data)
        consolidated_file = os.path.join(all_proteins_dir, "all_enriched_proteins.csv")
        df_consolidated.to_csv(consolidated_file, sep=';', index=False, encoding='utf-8')
        print(f"Arquivo consolidado salvo: {consolidated_file}")
    
    # 2. Salvar resultados por categoria
    for category in categories.keys():
        cat_dir = os.path.join(all_proteins_dir, category)
        os.makedirs(cat_dir, exist_ok=True)
        
        # Salvar arquivo consolidado por categoria
        if categories[category]:
            df = pd.DataFrame(categories[category])
            # Remover coluna 'Expressions' para salvar como CSV
            df_to_save = df.drop(columns=['Expressions'], errors='ignore')
            output_file = os.path.join(cat_dir, f"{category}_All.csv")
            df_to_save.to_csv(output_file, index=False, encoding='utf-8')
            print(f"Salvo: {output_file}")
    
    return all_proteins_dir

def save_membrane_interactors_results(categories, membrane_genes, interactors_genes, output_base_dir):
    """
    Salva resultados filtrados na pasta Membrane_Interactors_Proteins
    """
    membrane_interactors_dir = os.path.join(output_base_dir, "Membrane_Interactors_Proteins")
    os.makedirs(membrane_interactors_dir, exist_ok=True)
    
    filtered_stats = {
        'membrane': {'Enriquecido10X': 0, 'Enriquecido5X': 0, 'Enriquecido2X': 0},
        'interactors': {'Enriquecido10X': 0, 'Enriquecido5X': 0, 'Enriquecido2X': 0}
    }
    
    tissue_distribution = {
        'membrane': {},
        'interactors': {}
    }
    
    all_membrane_genes = []
    all_interactors_genes = []
    
    for category in categories:
        if categories[category]:
            # Criar subdiretórios para esta categoria
            category_dir = os.path.join(membrane_interactors_dir, category)
            os.makedirs(category_dir, exist_ok=True)
            
            # Filtrar genes de membrana
            membrane_filtered = [g for g in categories[category] if g['Ensembl'] in membrane_genes]
            filtered_stats['membrane'][category] = len(membrane_filtered)
            
            if membrane_filtered:
                df_membrane = pd.DataFrame(membrane_filtered)
                df_membrane = df_membrane.drop(columns=['Expressions'], errors='ignore')
                output_file = os.path.join(category_dir, f"{category}_Membrane.csv")
                df_membrane.to_csv(output_file, index=False, encoding='utf-8')
                
                # Adicionar à lista total de genes de membrana
                all_membrane_genes.extend(membrane_filtered)
                
                # Estatísticas de distribuição por tecido
                for gene in membrane_filtered:
                    tissue = gene['Tissue']
                    if tissue not in tissue_distribution['membrane']:
                        tissue_distribution['membrane'][tissue] = {cat: 0 for cat in categories}
                    tissue_distribution['membrane'][tissue][category] += 1
            
            # Filtrar genes interatores
            interactors_filtered = [g for g in categories[category] if g['Ensembl'] in interactors_genes]
            filtered_stats['interactors'][category] = len(interactors_filtered)
            
            if interactors_filtered:
                df_interactors = pd.DataFrame(interactors_filtered)
                df_interactors = df_interactors.drop(columns=['Expressions'], errors='ignore')
                output_file = os.path.join(category_dir, f"{category}_Interactors.csv")
                df_interactors.to_csv(output_file, index=False, encoding='utf-8')
                
                # Adicionar à lista total de genes interatores
                all_interactors_genes.extend(interactors_filtered)
                
                # Estatísticas de distribuição por tecido
                for gene in interactors_filtered:
                    tissue = gene['Tissue']
                    if tissue not in tissue_distribution['interactors']:
                        tissue_distribution['interactors'][tissue] = {cat: 0 for cat in categories}
                    tissue_distribution['interactors'][tissue][category] += 1
    
    # Salvar arquivos consolidados para membrana e interatores
    if all_membrane_genes:
        # Criar DataFrame consolidado para membrana
        membrane_consolidated = []
        for gene_info in all_membrane_genes:
            membrane_consolidated.append({
                'ID_Ensembl': gene_info['Ensembl'],
                'Tecido_enriquecido': gene_info['Tissue'],
                'Categoria': gene_info['Category']
            })
        
        df_membrane_consolidated = pd.DataFrame(membrane_consolidated)
        membrane_file = os.path.join(membrane_interactors_dir, "all_membrane_proteins.csv")
        df_membrane_consolidated.to_csv(membrane_file, sep=';', index=False, encoding='utf-8')
        print(f"Arquivo consolidado de membrana salvo: {membrane_file}")
    
    if all_interactors_genes:
        # Criar DataFrame consolidado para interatores
        interactors_consolidated = []
        for gene_info in all_interactors_genes:
            interactors_consolidated.append({
                'ID_Ensembl': gene_info['Ensembl'],
                'Tecido_enriquecido': gene_info['Tissue'],
                'Categoria': gene_info['Category']
            })
        
        df_interactors_consolidated = pd.DataFrame(interactors_consolidated)
        interactors_file = os.path.join(membrane_interactors_dir, "all_interactors_proteins.csv")
        df_interactors_consolidated.to_csv(interactors_file, sep=';', index=False, encoding='utf-8')
        print(f"Arquivo consolidado de interatores salvo: {interactors_file}")
    
    return filtered_stats, tissue_distribution, membrane_interactors_dir

def generate_log_file(categories, categorized_info, all_enriched_genes, filtered_stats, tissue_distribution, output_dir):
    """
    Gera arquivo de log com estatísticas detalhadas
    """
    log_lines = []
    log_lines.append("="*80)
    log_lines.append("RELATÓRIO DE ANÁLISE DE ENRIQUECIMENTO DE EXPRESSÃO")
    log_lines.append("="*80)
    
    # Estatísticas gerais
    log_lines.append("\nRESUMO GERAL")
    log_lines.append("-"*40)
    total_genes = len(all_enriched_genes)
    log_lines.append(f"Total de genes analisados: {total_genes}")
    
    for category in categories:
        log_lines.append(f"Genes {category}: {len(categories[category])}")
    
    # Distribuição por fonte (arquivo de origem)
    log_lines.append("\n\nDISTRIBUIÇÃO POR FONTE (ARQUIVO DE ORIGEM)")
    log_lines.append("-"*50)
    
    source_distribution = {}
    for category in categories:
        for gene_info in categorized_info[category]:
            source = gene_info['Source']
            if source not in source_distribution:
                source_distribution[source] = {cat: 0 for cat in categories}
            source_distribution[source][category] += 1
    
    for source, counts in source_distribution.items():
        log_lines.append(f"\n{source}:")
        for category in categories:
            if counts[category] > 0:
                log_lines.append(f"  {category}: {counts[category]}")
    
    # Distribuição por tecido
    log_lines.append("\n\nDISTRIBUIÇÃO POR TECIDO ENRIQUECIDO (GERAL)")
    log_lines.append("-"*50)
    
    general_tissue_dist = {}
    for gene_info in all_enriched_genes:
        tissue = gene_info['Tissue']
        category = gene_info['Category']
        if tissue not in general_tissue_dist:
            general_tissue_dist[tissue] = {cat: 0 for cat in categories}
        general_tissue_dist[tissue][category] += 1
    
    # Ordenar tecidos pelo total
    sorted_tissues = sorted(general_tissue_dist.items(), 
                          key=lambda x: sum(x[1].values()), 
                          reverse=True)
    
    for tissue, counts in sorted_tissues:
        total = sum(counts.values())
        log_lines.append(f"\n{tissue} (Total: {total}):")
        for category in categories:
            if counts[category] > 0:
                log_lines.append(f"  {category}: {counts[category]}")
    
    # Estatísticas de membrana
    log_lines.append("\n\nRESUMO MEMBRANA")
    log_lines.append("-"*40)
    for category in categories:
        log_lines.append(f"Genes {category}: {filtered_stats['membrane'][category]}")
    
    if tissue_distribution['membrane']:
        log_lines.append("\nDistribuição por tecido (Membrana):")
        for tissue, counts in tissue_distribution['membrane'].items():
            total = sum(counts.values())
            if total > 0:
                log_lines.append(f"  {tissue}: {total} genes")
    
    # Estatísticas de interatores
    log_lines.append("\n\nRESUMO INTERACTORAS")
    log_lines.append("-"*40)
    for category in categories:
        log_lines.append(f"Genes {category}: {filtered_stats['interactors'][category]}")
    
    if tissue_distribution['interactors']:
        log_lines.append("\nDistribuição por tecido (Interatores):")
        for tissue, counts in tissue_distribution['interactors'].items():
            total = sum(counts.values())
            if total > 0:
                log_lines.append(f"  {tissue}: {total} genes")
    
    # Lista detalhada de genes por categoria
    log_lines.append("\n\nDETALHAMENTO POR CATEGORIA")
    log_lines.append("="*80)
    
    for category in categories:
        if categorized_info[category]:
            log_lines.append(f"\n{category} ({len(categorized_info[category])} genes):")
            log_lines.append("-"*60)
            
            # Agrupar por tecido
            tissue_groups = {}
            for gene_info in categorized_info[category]:
                tissue = gene_info['Tissue']
                if tissue not in tissue_groups:
                    tissue_groups[tissue] = []
                tissue_groups[tissue].append(gene_info)
            
            for tissue, genes in tissue_groups.items():
                log_lines.append(f"\n  Tecido: {tissue} ({len(genes)} genes)")
                for gene in genes[:10]:  # Mostrar apenas primeiros 10 por tecido
                    log_lines.append(f"    {gene['Ensembl']} ({gene['Gene']}) - Ratio: {gene['Ratio']:.2f}X, "
                                   f"Expressão: {gene['Max_Expression']:.2f}, Fonte: {gene['Source']}")
                
                if len(genes) > 10:
                    log_lines.append(f"    ... e mais {len(genes) - 10} genes")
    
    # Salvar arquivo de log
    log_file = os.path.join(output_dir, "enrichment_analysis_log.txt")
    with open(log_file, 'w', encoding='utf-8') as f:
        f.write('\n'.join(log_lines))
    
    print(f"\nArquivo de log salvo: {log_file}")
    
    # Também salvar estatísticas em CSV
    stats_data = []
    for gene_info in all_enriched_genes:
        stats_data.append({
            'Category': gene_info['Category'],
            'Ensembl': gene_info['Ensembl'],
            'Gene': gene_info['Gene'],
            'Tissue': gene_info['Tissue'],
            'Ratio': gene_info['Ratio'],
            'Max_Expression': gene_info['Max_Expression'],
            'Source': gene_info['Source'],
            'Uniprot': gene_info['Uniprot'],
            'Num_Tissues': gene_info['Num_Tissues']
        })
    
    if stats_data:
        stats_df = pd.DataFrame(stats_data)
        stats_file = os.path.join(output_dir, "detailed_statistics.csv")
        stats_df.to_csv(stats_file, index=False, encoding='utf-8')
        print(f"Estatísticas detalhadas salvas: {stats_file}")
    
    return log_file

def main():
    """
    Função principal
    """
    # Configuração dos arquivos de entrada
    input_files = [
        "Unique_Prostate.tsv",
        "Unique_Prostate_Cancer.tsv",
        "Shared_Prostate.tsv",
        "Shared_Prostate_Cancer.tsv"
    ]
    
    # Arquivos de filtro
    membrane_file = "ENSG_Membrane_Proteins.txt"
    interactors_file = "ENSG_Interactors_Proteins.txt"
    
    # Diretório de saída principal
    output_base_dir = "Enrichment_Analysis_Results"
    
    print("="*80)
    print("ANÁLISE DE ENRIQUECIMENTO DE EXPRESSÃO - HUMAN PROTEIN ATLAS")
    print("="*80)
    
    # 1. Criar diretório de saída principal
    os.makedirs(output_base_dir, exist_ok=True)
    
    # 2. Processar arquivos
    print("\n1. Processando arquivos de entrada...")
    all_genes = process_files(input_files)
    
    if not all_genes:
        print("Nenhum gene foi processado. Verifique os arquivos de entrada.")
        return
    
    # 3. Classificar genes
    print("\n2. Classificando genes por enriquecimento...")
    categories, categorized_info, all_enriched_genes = classify_genes(all_genes)
    
    # 4. Salvar resultados gerais (All_Proteins)
    print("\n3. Salvando resultados gerais em All_Proteins...")
    all_proteins_dir = save_all_proteins_results(categories, all_enriched_genes, output_base_dir)
    
    # 5. Carregar arquivos de filtro
    print("\n4. Carregando arquivos de filtro...")
    membrane_genes, interactors_genes = load_filter_files(membrane_file, interactors_file)
    
    # 6. Filtrar e salvar resultados filtrados (Membrane_Interactors_Proteins)
    print("\n5. Filtrando e salvando resultados em Membrane_Interactors_Proteins...")
    filtered_stats, tissue_distribution, membrane_interactors_dir = save_membrane_interactors_results(
        categories, membrane_genes, interactors_genes, output_base_dir
    )
    
    # 7. Gerar arquivo de log
    print("\n6. Gerando relatório completo...")
    generate_log_file(categories, categorized_info, all_enriched_genes, filtered_stats, 
                     tissue_distribution, output_base_dir)
    
    print("\n" + "="*80)
    print("PROCESSAMENTO CONCLUÍDO!")
    print("="*80)
    print(f"\nResultados salvos em: {output_base_dir}")
    print("\nEstrutura criada:")
    print(f"  {output_base_dir}/")
    print(f"  ├── All_Proteins/")
    print(f"  │   ├── all_enriched_proteins.csv (todos os genes enriquecidos)")
    print(f"  │   ├── Enriquecido10X/")
    print(f"  │   │   └── Enriquecido10X_All.csv")
    print(f"  │   ├── Enriquecido5X/")
    print(f"  │   │   └── Enriquecido5X_All.csv")
    print(f"  │   └── Enriquecido2X/")
    print(f"  │       └── Enriquecido2X_All.csv")
    print(f"  ├── Membrane_Interactors_Proteins/")
    print(f"  │   ├── all_membrane_proteins.csv")
    print(f"  │   ├── all_interactors_proteins.csv")
    print(f"  │   ├── Enriquecido10X/")
    print(f"  │   │   ├── Enriquecido10X_Membrane.csv")
    print(f"  │   │   └── Enriquecido10X_Interactors.csv")
    print(f"  │   ├── Enriquecido5X/")
    print(f"  │   │   ├── Enriquecido5X_Membrane.csv")
    print(f"  │   │   └── Enriquecido5X_Interactors.csv")
    print(f"  │   └── Enriquecido2X/")
    print(f"  │       ├── Enriquecido2X_Membrane.csv")
    print(f"  │       └── Enriquecido2X_Interactors.csv")
    print(f"  ├── enrichment_analysis_log.txt")
    print(f"  └── detailed_statistics.csv")

if __name__ == "__main__":
    main()