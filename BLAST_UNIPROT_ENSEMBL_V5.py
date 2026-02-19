#!/usr/bin/env python3
import subprocess
import os
import sys
import argparse
from datetime import datetime
from collections import defaultdict

def parse_ensembl_header(header):
    """Parseia o header do Ensembl e retorna (gene_id, transcript_id, protein_id)"""
    parts = header.split('|')
    gene_id = parts[0] if len(parts) > 0 else 'NA'
    transcript_id = parts[1] if len(parts) > 1 else 'NA'
    protein_id = parts[2].split('.')[0] if len(parts) > 2 else 'NA'  # Remove versão se existir
    return gene_id, transcript_id, protein_id

def parse_uniprot_header(header):
    """Parseia o header do UniProt e retorna uniprot_id"""
    parts = header.split()
    uniprot_id = parts[0] if len(parts) > 0 else 'NA'
    return uniprot_id

def create_blast_db(fasta_file, db_name):
    """Cria banco de dados BLAST"""
    cmd = f"makeblastdb -in {fasta_file} -dbtype prot -out {db_name}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Erro ao criar banco de dados BLAST: {result.stderr}")
        sys.exit(1)

def run_blast(query_file, db_file, output_file, threads):
    """Executa BLASTp e retorna resultados"""
    cmd = f"blastp -query {query_file} -db {db_file} -out {output_file} -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' -num_threads {threads}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Erro ao executar BLAST: {result.stderr}")
        sys.exit(1)

def calculate_coverage(qstart, qend, qlen, sstart, send, slen):
    """Calcula a porcentagem de coverage para query e subject"""
    query_coverage = (abs(qend - qstart) + 1) / qlen * 100
    subject_coverage = (abs(send - sstart) + 1) / slen * 100
    return query_coverage, subject_coverage

def main():
    # Configuração do parser de argumentos
    parser = argparse.ArgumentParser(description='Correlaciona IDs UniProt com sequências Ensembl usando BLAST')
    parser.add_argument('-t', '--threads', type=int, default=23, 
                       help='Número de threads para o BLAST (padrão: 1)')
    parser.add_argument('-e', '--ensembl', type=str, default='ENSEMBL.fasta',
                       help='Arquivo multifasta do Ensembl (padrão: EXEMPLO_ENSEMBL.txt)')
    parser.add_argument('-u', '--uniprot', type=str, default='3_UniProt_Sequences.fasta',
                       help='Arquivo multifasta do UniProt (padrão: EXEMPLO_UniProt_Sequences.txt)')
    
    args = parser.parse_args()
    
    # Configuração de arquivos
    ensembl_file = args.ensembl
    uniprot_file = args.uniprot
    blast_db = "ensembl_blast_db"
    blast_results = "blast_results.txt"
    threads = args.threads
    
    print(f"Usando {threads} thread(s) para o BLAST")
    
    # Verificar se os arquivos existem
    if not os.path.exists(ensembl_file):
        print(f"Erro: Arquivo {ensembl_file} não encontrado!")
        sys.exit(1)
    
    if not os.path.exists(uniprot_file):
        print(f"Erro: Arquivo {uniprot_file} não encontrado!")
        sys.exit(1)
    
    # Criar dicionário de mapeamento para headers do Ensembl
    ensembl_map = {}
    ensembl_lengths = {}
    print("Processando arquivo Ensembl...")
    with open(ensembl_file) as f:
        current_header = None
        sequence = ""
        for line in f:
            if line.startswith('>'):
                if current_header and sequence:
                    gene_id, transcript_id, protein_id = parse_ensembl_header(current_header)
                    ensembl_map[current_header] = (gene_id, transcript_id, protein_id)
                    ensembl_lengths[current_header] = len(sequence)
                current_header = line[1:].strip()
                sequence = ""
            else:
                sequence += line.strip()
        
        # Processar a última sequência
        if current_header and sequence:
            gene_id, transcript_id, protein_id = parse_ensembl_header(current_header)
            ensembl_map[current_header] = (gene_id, transcript_id, protein_id)
            ensembl_lengths[current_header] = len(sequence)

    # Processar arquivo UniProt e coletar IDs
    uniprot_ids = set()
    uniprot_lengths = {}
    print("Processando arquivo UniProt...")
    
    with open(uniprot_file) as f:
        current_header = None
        sequence = ""
        for line in f:
            if line.startswith('>'):
                if current_header:
                    uniprot_id = parse_uniprot_header(current_header)
                    uniprot_ids.add(uniprot_id)
                    uniprot_lengths[uniprot_id] = len(sequence)
                current_header = line[1:].strip()
                sequence = ""
            else:
                sequence += line.strip()
        
        # Processar a última sequência
        if current_header:
            uniprot_id = parse_uniprot_header(current_header)
            uniprot_ids.add(uniprot_id)
            uniprot_lengths[uniprot_id] = len(sequence)

    # VERIFICAR SE O ARQUIVO DE RESULTADOS DO BLAST JÁ EXISTE
    if os.path.exists(blast_results):
        print(f"Arquivo {blast_results} encontrado. Pulando etapa do BLAST...")
    else:
        # Criar banco de dados BLAST e executar
        print("Criando banco de dados BLAST...")
        create_blast_db(ensembl_file, blast_db)
        
        print(f"Executando BLAST com {threads} thread(s)...")
        run_blast(uniprot_file, blast_db, blast_results, threads)

    # Processar resultados do BLAST
    blast_matches = defaultdict(list)
    print("Processando resultados do BLAST...")
    
    with open(blast_results) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 14:
                uniprot_id = parts[0]
                ensembl_header = parts[1]
                identity = float(parts[2])
                alignment_length = int(parts[3])
                qstart = int(parts[6])
                qend = int(parts[7])
                sstart = int(parts[8])
                send = int(parts[9])
                bitscore = float(parts[11])
                qlen = int(parts[12])
                slen = int(parts[13])
                
                # Aceitar identidade de 99% ou 100%
                if identity >= 99.0:
                    if ensembl_header in ensembl_map:
                        # Calcular coverage e aplicar filtro
                        query_coverage, subject_coverage = calculate_coverage(
                            qstart, qend, qlen, sstart, send, slen
                        )
                        
                        # MODIFICAÇÃO: Aceitar coverage a partir de 90%
                        if query_coverage >= 90.0 and subject_coverage >= 90.0:
                            gene_id, transcript_id, protein_id = ensembl_map[ensembl_header]
                            blast_matches[uniprot_id].append({
                                'gene_id': gene_id,
                                'transcript_id': transcript_id,
                                'protein_id': protein_id,
                                'ensembl_header': ensembl_header,
                                'alignment_length': alignment_length,
                                'bitscore': bitscore,
                                'identity': identity,
                                'query_coverage': query_coverage,
                                'subject_coverage': subject_coverage
                            })

    # Processar correspondências
    correlation_lines = []
    multiple_matches_lines = []
    
    stats = {
        'total': len(uniprot_ids),
        'found_single': 0,
        'found_multiple': 0,
        'found_gene': set(),
        'multiple_matches': 0,
        'unique_genes': set()
    }

    print("Processando correlações...")
    
    # Para cada UniProt, encontrar o melhor gene
    for uniprot_id in uniprot_ids:
        if uniprot_id in blast_matches:
            matches = blast_matches[uniprot_id]
            
            # Agrupar por gene e manter apenas a melhor entrada por gene
            gene_best_matches = {}
            for match in matches:
                gene_id = match['gene_id']
                if gene_id not in gene_best_matches:
                    gene_best_matches[gene_id] = match
                else:
                    # Manter apenas a melhor entrada para este gene
                    current_best = gene_best_matches[gene_id]
                    # Critérios: identidade > alignment_length > bitscore
                    if (match['identity'] > current_best['identity'] or
                        (match['identity'] == current_best['identity'] and 
                         match['alignment_length'] > current_best['alignment_length']) or
                        (match['identity'] == current_best['identity'] and 
                         match['alignment_length'] == current_best['alignment_length'] and 
                         match['bitscore'] > current_best['bitscore'])):
                        gene_best_matches[gene_id] = match
            
            # Se há múltiplos genes únicos, considerar como múltipla correspondência
            unique_genes = list(gene_best_matches.keys())
            if len(unique_genes) > 1:
                stats['multiple_matches'] += 1
                stats['found_multiple'] += 1
                
                # Adicionar todas as correspondências ao arquivo de múltiplas correspondências
                for gene_id, match in gene_best_matches.items():
                    multiple_matches_lines.append([
                        uniprot_id,
                        match['gene_id'],
                        match['transcript_id'],
                        match['protein_id'],
                        "BLAST",
                        str(match['alignment_length']),
                        str(match['bitscore']),
                        f"{match['identity']:.1f}%",
                        f"{match['query_coverage']:.1f}%",
                        f"{match['subject_coverage']:.1f}%"
                    ])
                
                # No arquivo principal, colocar apenas a melhor correspondência
                best_gene_id = None
                best_match = None
                for gene_id, match in gene_best_matches.items():
                    if best_match is None:
                        best_gene_id = gene_id
                        best_match = match
                    else:
                        # Comparar com a hierarquia de critérios
                        if (match['identity'] > best_match['identity'] or
                            (match['identity'] == best_match['identity'] and 
                             match['alignment_length'] > best_match['alignment_length']) or
                            (match['identity'] == best_match['identity'] and 
                             match['alignment_length'] == best_match['alignment_length'] and 
                             match['bitscore'] > best_match['bitscore'])):
                            best_gene_id = gene_id
                            best_match = match
                
                correlation_lines.append([
                    uniprot_id,
                    best_match['gene_id'],
                    best_match['transcript_id'],
                    best_match['protein_id'],
                    "BLAST"
                ])
                stats['found_single'] += 1
                stats['found_gene'].add(best_gene_id)
                stats['unique_genes'].add(best_gene_id)
                
            else:
                # Apenas uma correspondência
                stats['found_single'] += 1
                gene_id = unique_genes[0]
                match = gene_best_matches[gene_id]
                correlation_lines.append([
                    uniprot_id,
                    match['gene_id'],
                    match['transcript_id'],
                    match['protein_id'],
                    "BLAST"
                ])
                stats['found_gene'].add(gene_id)
                stats['unique_genes'].add(gene_id)

    # MODIFICAÇÃO 1: Escrever arquivo de correlações principal como TSV com 5 colunas
    with open("uniprot_ensembl_correlations.tsv", "w") as f:
        f.write("UNIPROT\tENSEMBL(Gene)\tENSEMBL(Transcrito)\tENSEMBL(Proteina)\tOrigem\n")
        for line in correlation_lines:
            f.write("\t".join(line) + "\n")

    # MODIFICAÇÃO 1: Escrever arquivo de múltiplas correspondências como TSV com colunas adicionais
    with open("uniprot_ensembl_multiple_matches.tsv", "w") as f:
        f.write("UNIPROT\tENSEMBL(Gene)\tENSEMBL(Transcrito)\tENSEMBL(Proteina)\tOrigem\tAlignment_Length\tBitscore\tIdentity\tQuery_Coverage\tSubject_Coverage\n")
        for line in multiple_matches_lines:
            f.write("\t".join(line) + "\n")

    # Calcular estatísticas finais
    stats['found_any'] = stats['found_single'] + stats['found_multiple']
    
    # Gerar estatísticas
    with open("uniprot_ensembl_stats.txt", "w") as f:
        f.write("ESTATÍSTICAS DE CORRESPONDÊNCIA UNIPROT->ENSEMBL (BLAST)\n")
        f.write("=" * 55 + "\n")
        f.write(f"Data de processamento: {datetime.now().strftime('%c')}\n")
        f.write(f"Arquivo de entrada UniProt: {uniprot_file}\n")
        f.write(f"Arquivo de entrada Ensembl: {ensembl_file}\n")
        f.write(f"Threads utilizadas: {threads}\n")
        f.write(f"Total de IDs processados: {stats['total']}\n")
        f.write(f"IDs com correspondência única: {stats['found_single']}\n")
        f.write(f"IDs com múltiplas correspondências: {stats['found_multiple']}\n")
        f.write(f"IDs com correspondência encontrada: {stats['found_single']}\n")
        f.write(f"IDs sem correspondência: {stats['total'] - stats['found_any']}\n")
        f.write(f"Genes Ensembl únicos encontrados: {len(stats['unique_genes'])}\n")
        
        success_rate = (stats['found_single'] / stats['total'] * 100) if stats['total'] > 0 else 0
        f.write(f"Taxa de sucesso: {success_rate:.1f}%\n\n")
        
        f.write("CRITÉRIOS UTILIZADOS:\n")
        f.write("- Identidade mínima: 99%\n")
        f.write("- Query coverage mínimo: 90%\n")  # MODIFICAÇÃO 2
        f.write("- Subject coverage mínimo: 90%\n")  # MODIFICAÇÃO 2
        f.write("- Ordem de prioridade: Identidade > Alignment Length > Bitscore\n")
        f.write("- Múltiplas UniProts podem correlacionar com o mesmo gene Ensembl (splicing alternativo)\n\n")
        
        f.write("ARQUIVOS GERADOS:\n")
        f.write("- uniprot_ensembl_correlations.tsv: Correlações únicas (melhor match por UniProt)\n")
        f.write("- uniprot_ensembl_multiple_matches.tsv: UniProts com múltiplos genes Ensembl\n")
        f.write("- uniprot_ensembl_stats.txt: Estatísticas do processamento\n")

    # Limpeza (apenas se o banco de dados foi criado nesta execução)
    if not os.path.exists(blast_results):
        for ext in ['phr', 'pin', 'pog', 'psd', 'psi', 'psq']:
            file = f"{blast_db}.{ext}"
            if os.path.exists(file):
                os.remove(file)
    
    print("\nProcessamento concluído!")
    print(f"Total de IDs processados: {stats['total']}")
    print(f"IDs com correspondência única: {stats['found_single']}")
    print(f"IDs com múltiplas correspondências: {stats['found_multiple']}")
    print(f"IDs sem correspondência: {stats['total'] - stats['found_any']}")
    print(f"Genes Ensembl únicos encontrados: {len(stats['unique_genes'])}")
    print(f"Taxa de sucesso: {success_rate:.1f}%")
    print(f"\nCritérios utilizados:")
    print("- Identidade mínima: 99%")
    print("- Query coverage mínimo: 90%")  # MODIFICAÇÃO 2
    print("- Subject coverage mínimo: 90%")  # MODIFICAÇÃO 2
    print("- Múltiplas UniProts podem correlacionar com o mesmo gene Ensembl")
    print(f"\nArquivos gerados:")
    print("- uniprot_ensembl_correlations.tsv")
    print("- uniprot_ensembl_multiple_matches.tsv")
    print("- uniprot_ensembl_stats.txt")

if __name__ == "__main__":
    main()