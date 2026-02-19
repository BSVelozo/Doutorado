import glob
import time
from collections import defaultdict
from datetime import datetime
import subprocess
import os

# Ordem hierárquica das origens - incluindo Manualmente
HIERARQUIA = {'NCBI': 0, 'ENSEMBL': 1, 'UNIPROT': 2, 'Human Protein Atlas': 3, 'BLAST': 4, 'Manualmente': 5}

# Total fixo de proteínas UNIPROT para cálculo das porcentagens
TOTAL_PROTEINAS_UNIPROT = 79658

def processar_arquivos():
    # Lista específica de arquivos a processar - incluindo todos os 6 arquivos
    arquivos = [
        '1_uniprot_ensembl_found_NCBI.tsv',
        '2_uniprot_ensembl_found_BIOMART.tsv',
        '3_UniProt_ensembl_found_UNIPROT.tsv',
        '4_UniProt_ensembl_found_HPA.tsv',
        '5_uniprot_ensembl_correlations_BLAST.tsv',
        '6_NCBI_PDB_ensembl_found.tsv'
    ]
    
    # Mapear arquivos para seus tipos (quantidade de colunas)
    tipos_arquivos = {
        '1_uniprot_ensembl_found_NCBI.tsv': 5,  # UNIPROT, Gene, Transcrito, Proteína, Origem
        '2_uniprot_ensembl_found_BIOMART.tsv': 5,  # UNIPROT, Gene, Transcrito, Proteína, Origem  
        '3_UniProt_ensembl_found_UNIPROT.tsv': 5,  # UNIPROT, Gene, Transcrito, Proteína, Origem
        '4_UniProt_ensembl_found_HPA.tsv': 3,  # UNIPROT, Gene, Origem
        '5_uniprot_ensembl_correlations_BLAST.tsv': 5,  # UNIPROT, Gene, Transcrito, Proteína, Origem
        '6_NCBI_PDB_ensembl_found.tsv': 3  # ID_NCBI_PDB, Gene, Origem
    }
    
    # Verificar se todos os arquivos existem
    arquivos_encontrados = []
    contagem_por_arquivo = defaultdict(int)
    
    # Primeiro, contar quantas entradas cada arquivo tem originalmente
    for arquivo in arquivos:
        try:
            with open(arquivo, 'r') as f:
                cabecalho = next(f)  # Pular cabeçalho
                contador = 0
                for linha in f:
                    partes = linha.strip().split('\t')
                    if len(partes) >= 2 and partes[0]:  # Verificar se tem ID
                        contador += 1
                contagem_por_arquivo[arquivo] = contador
                arquivos_encontrados.append(arquivo)
        except FileNotFoundError:
            print(f"AVISO: Arquivo {arquivo} não encontrado. Continuando com os arquivos disponíveis.")
    
    # Dicionário para armazenar a melhor entrada de cada ID
    dados_combinados = {}
    contagem_origem_final = defaultdict(int)
    
    for arquivo in arquivos_encontrados:
        try:
            num_colunas = tipos_arquivos[arquivo]
            with open(arquivo, 'r') as f:
                cabecalho = next(f)  # Pular cabeçalho
                for linha in f:
                    partes = linha.strip().split('\t')
                    if len(partes) < 2 or not partes[0]:
                        continue
                        
                    id_entrada = partes[0]
                    
                    # Processar de acordo com o número de colunas
                    if num_colunas == 5:
                        # Formato: ID, Gene, Transcrito, Proteína, Origem
                        gene = partes[1] if len(partes) > 1 else ''
                        transcrito = partes[2] if len(partes) > 2 else ''
                        proteina = partes[3] if len(partes) > 3 else ''
                        origem = partes[4].strip() if len(partes) > 4 else 'Desconhecida'
                        entrada_formatada = [id_entrada, gene, transcrito, proteina, origem]
                        
                    elif num_colunas == 3:
                        # Formato: ID, Gene, Origem
                        gene = partes[1] if len(partes) > 1 else ''
                        origem = partes[2].strip() if len(partes) > 2 else 'Desconhecida'
                        # Para arquivos com 3 colunas, transcrito e proteína são vazios
                        entrada_formatada = [id_entrada, gene, '', '', origem]
                    
                    # Verificar se a origem está na hierarquia
                    if origem not in HIERARQUIA:
                        print(f"AVISO: Origem '{origem}' não encontrada na hierarquia. Pulando entrada.")
                        continue
                    
                    # Manter apenas a entrada com maior prioridade
                    if id_entrada not in dados_combinados or \
                       HIERARQUIA[origem] < HIERARQUIA[dados_combinados[id_entrada][4]]:
                        dados_combinados[id_entrada] = entrada_formatada
                        
        except FileNotFoundError:
            continue  # Pular arquivos não encontrados
    
    # Contar origens finais
    for entrada in dados_combinados.values():
        origem = entrada[4]
        contagem_origem_final[origem] += 1

    return dados_combinados, contagem_origem_final, arquivos_encontrados, contagem_por_arquivo

def gerar_estatisticas(dados_combinados):
    stats = {
        'gene': 0,
        'transcrito': 0,
        'proteina': 0,
        'com_correspondencia': 0,
        'total_uniprot': 0  # Contar apenas UNIPROTs
    }
    
    for entrada in dados_combinados.values():
        id_entrada = entrada[0]
        
        # Verificar valores vazios ou inválidos
        tem_gene = entrada[1] not in ['', 'NA', 'N/A', '-', '.', 'null']
        tem_transcrito = entrada[2] not in ['', 'NA', 'N/A', '-', '.', 'null']
        tem_proteina = entrada[3] not in ['', 'NA', 'N/A', '-', '.', 'null']
        
        if tem_gene:
            stats['gene'] += 1
        if tem_transcrito:
            stats['transcrito'] += 1
        if tem_proteina:
            stats['proteina'] += 1
        if any([tem_gene, tem_transcrito, tem_proteina]):
            stats['com_correspondencia'] += 1
        
        # Contar apenas UNIPROTs (IDs que seguem o padrão UNIPROT)
        if id_entrada.startswith(('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 
                                 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z')):
            stats['total_uniprot'] += 1

    return stats

def extrair_ensg_unicos_combinados(arquivo_combinado, arquivo_ensembl_previa):
    """Extrai IDs ENSG únicos combinando o arquivo combinado com o arquivo de ENSG prévio"""
    ensg_ids_set = set()
    
    # 1. Extrair do arquivo combinado usando awk
    try:
        # Usar awk para extrair a coluna 2 (ENSEMBL Gene) e remover duplicatas
        comando = f"awk -F'\t' 'NR>1 && $2 != \"\" && $2 != \"NA\" && $2 != \"-\" {{print $2}}' {arquivo_combinado} | sort | uniq"
        resultado = subprocess.run(comando, shell=True, capture_output=True, text=True)
        
        if resultado.returncode == 0:
            ensg_ids_comb = resultado.stdout.strip().split('\n')
            for id_ in ensg_ids_comb:
                if id_:
                    ensg_ids_set.add(id_)
        else:
            print(f"Erro ao executar awk: {resultado.stderr}")
    except Exception as e:
        print(f"Erro ao extrair ENSG do arquivo combinado: {e}")
    
    # 2. Ler arquivo de ENSG prévio (All_Unique_ensembl.txt)
    try:
        if os.path.exists(arquivo_ensembl_previa):
            with open(arquivo_ensembl_previa, 'r') as f:
                for linha in f:
                    ensg_id = linha.strip()
                    if ensg_id:
                        ensg_ids_set.add(ensg_id)
            print(f"Arquivo {arquivo_ensembl_previa} encontrado e processado.")
        else:
            print(f"AVISO: Arquivo {arquivo_ensembl_previa} não encontrado.")
    except Exception as e:
        print(f"Erro ao ler arquivo {arquivo_ensembl_previa}: {e}")
    
    # Converter para lista ordenada
    ensg_ids = sorted(list(ensg_ids_set))
    return len(ensg_ids), ensg_ids

def main():
    inicio = time.time()
    dados_combinados, contagem_origem_final, arquivos_encontrados, contagem_por_arquivo = processar_arquivos()
    stats = gerar_estatisticas(dados_combinados)
    
    # Mapear nomes de arquivos para origens
    arquivo_para_origem = {
        '1_uniprot_ensembl_found_NCBI.tsv': 'NCBI',
        '2_uniprot_ensembl_found_BIOMART.tsv': 'ENSEMBL',
        '3_UniProt_ensembl_found_UNIPROT.tsv': 'UNIPROT',
        '4_UniProt_ensembl_found_HPA.tsv': 'Human Protein Atlas',
        '5_uniprot_ensembl_correlations_BLAST.tsv': 'BLAST',
        '6_NCBI_PDB_ensembl_found.tsv': 'Manualmente'
    }
    
    # Escrever arquivo combinado (uma entrada por ID)
    arquivo_combinado = 'UNIPROT_ENSEMBL_ALL.tsv'
    with open(arquivo_combinado, 'w') as f:
        f.write("ID\tENSEMBL(Gene)\tENSEMBL(Transcrito)\tENSEMBL(Proteína)\tOrigem\n")
        for entrada in dados_combinados.values():
            linha_formatada = [entrada[0], entrada[1], entrada[2], entrada[3], entrada[4]]
            f.write('\t'.join(linha_formatada) + '\n')
    
    # Escrever arquivo com correspondências (todas têm correspondência)
    with open('uniprot_ensembl_found.txt', 'w') as f:
        f.write("ID\tENSEMBL(Gene)\tENSEMBL(Transcrito)\tENSEMBL(Proteína)\tOrigem\n")
        for entrada in dados_combinados.values():
            linha_formatada = [entrada[0], entrada[1], entrada[2], entrada[3], entrada[4]]
            f.write('\t'.join(linha_formatada) + '\n')
    
    # Extrair ENSG únicos combinando arquivo combinado com arquivo prévio
    print("\nExtraindo IDs ENSG únicos (combinados)...")
    arquivo_ensembl_previa = 'All_Unique_ensembl.txt'
    num_ensg_unicos, ensg_ids = extrair_ensg_unicos_combinados(arquivo_combinado, arquivo_ensembl_previa)
    
    # Salvar IDs ENSG únicos combinados em arquivo
    arquivo_ensg_combinados = 'ensg_ids_completos.txt'
    with open(arquivo_ensg_combinados, 'w') as f:
        for ensg_id in ensg_ids:
            f.write(f"{ensg_id}\n")
    
    # Calcular estatísticas de fontes
    num_do_combinado = len(dados_combinados)
    num_do_previa = 0
    if os.path.exists(arquivo_ensembl_previa):
        with open(arquivo_ensembl_previa, 'r') as f:
            num_do_previa = sum(1 for linha in f if linha.strip())
    
    # Gerar arquivo de estatísticas detalhado
    data_processamento = datetime.now().strftime("%a %d %b %Y %H:%M:%S -03")
    with open('estatisticas_detalhadas.txt', 'w') as estat:
        estat.write("ESTATÍSTICAS DETALHADAS - CORRESPONDÊNCIA ID->ENSEMBL\n")
        estat.write("=" * 60 + "\n")
        estat.write(f"Data de processamento: {data_processamento}\n")
        estat.write(f"Arquivos processados: {', '.join(arquivos_encontrados)}\n")
        estat.write(f"Total de proteínas UNIPROT (referência): {TOTAL_PROTEINAS_UNIPROT}\n")
        estat.write(f"IDs únicos processados: {len(dados_combinados)}\n")
        estat.write(f"IDs ENSG únicos combinados: {num_ensg_unicos}\n\n")
        
        estat.write("FONTES DOS IDs ENSG:\n")
        estat.write("-" * 25 + "\n")
        estat.write(f"IDs do arquivo combinado: {num_do_combinado}\n")
        estat.write(f"IDs do arquivo All_Unique_ensembl.txt: {num_do_previa}\n")
        estat.write(f"Total combinado (sem redundância): {num_ensg_unicos}\n\n")
        
        estat.write("CORRESPONDÊNCIAS ENSEMBL:\n")
        estat.write("-" * 30 + "\n")
        estat.write(f"Ensembl Gene: {stats['gene']} de {stats['total_uniprot']} ({stats['gene']/stats['total_uniprot']*100:.2f}%)\n")
        estat.write(f"Ensembl Transcript: {stats['transcrito']} de {stats['total_uniprot']} ({stats['transcrito']/stats['total_uniprot']*100:.2f}%)\n")
        estat.write(f"Ensembl Protein: {stats['proteina']} de {stats['total_uniprot']} ({stats['proteina']/stats['total_uniprot']*100:.2f}%)\n")
        estat.write(f"Pelo menos uma correspondência: {stats['com_correspondencia']} de {stats['total_uniprot']} ({stats['com_correspondencia']/stats['total_uniprot']*100:.2f}%)\n")
        estat.write(f"Nenhuma correspondência: {stats['total_uniprot'] - stats['com_correspondencia']} de {stats['total_uniprot']} ({(stats['total_uniprot'] - stats['com_correspondencia'])/stats['total_uniprot']*100:.2f}%)\n\n")
        
        estat.write("DISTRIBUIÇÃO POR ORIGEM (APÓS HIERARQUIA):\n")
        estat.write("-" * 45 + "\n")
        
        # Calcular totais por origem
        totais_por_origem = {}
        for arquivo in arquivos_encontrados:
            origem = arquivo_para_origem.get(arquivo, arquivo)
            totais_por_origem[origem] = contagem_por_arquivo[arquivo]
        
        # Escrever estatísticas por origem
        for origem, count in sorted(contagem_origem_final.items(), key=lambda x: HIERARQUIA[x[0]]):
            total_na_origem = totais_por_origem.get(origem, 0)
            if total_na_origem > 0:
                porcentagem_preservada = count / total_na_origem * 100
                porcentagem_total = count / stats['total_uniprot'] * 100
                estat.write(f"{origem}: {count} de {total_na_origem} ({porcentagem_preservada:.2f}%) - {porcentagem_total:.2f}% do total UNIPROT\n")
            else:
                porcentagem_total = count / stats['total_uniprot'] * 100
                estat.write(f"{origem}: {count} de 0 (N/A) - {porcentagem_total:.2f}% do total UNIPROT\n")
        
        # Estatísticas adicionais
        estat.write(f"\nRESUMO:\n")
        estat.write(f"- IDs com correspondência: {stats['com_correspondencia']} ({stats['com_correspondencia']/stats['total_uniprot']*100:.2f}%)\n")
        estat.write(f"- IDs sem correspondência: {stats['total_uniprot'] - stats['com_correspondencia']} ({(stats['total_uniprot'] - stats['com_correspondencia'])/stats['total_uniprot']*100:.2f}%)\n")
        estat.write(f"- IDs processados (únicos): {len(dados_combinados)}\n")
        estat.write(f"- IDs UNIPROT processados: {stats['total_uniprot']}\n")
        estat.write(f"- IDs não-UNIPROT processados: {len(dados_combinados) - stats['total_uniprot']}\n")
        estat.write(f"- IDs ENSG únicos combinados: {num_ensg_unicos}\n")
        
        # Mostrar totais por arquivo
        estat.write(f"\nTOTAIS POR ARQUIVO (ANTES DA HIERARQUIA):\n")
        for arquivo in arquivos_encontrados:
            origem = arquivo_para_origem.get(arquivo, arquivo)
            total = contagem_por_arquivo.get(arquivo, 0)
            estat.write(f"- {origem} ({arquivo}): {total} entradas\n")
    
    # Gerar log resumido
    with open('processamento.log', 'w') as log:
        log.write("ESTATÍSTICAS DE CORRESPONDÊNCIA ID->ENSEMBL (FINAL)\n")
        log.write("=" * 55 + "\n")
        log.write(f"Data de processamento: {data_processamento}\n")
        log.write(f"Arquivos de entrada: {', '.join(arquivos_encontrados)}\n")
        log.write(f"Total de proteínas UNIPROT (referência): {TOTAL_PROTEINAS_UNIPROT}\n")
        log.write(f"IDs UNIPROT processados: {stats['total_uniprot']}\n")
        log.write(f"IDs ENSG únicos combinados: {num_ensg_unicos}\n\n")
        
        log.write("FONTES DOS IDs ENSG:\n")
        log.write(f"- IDs do arquivo combinado: {num_do_combinado}\n")
        log.write(f"- IDs do arquivo All_Unique_ensembl.txt: {num_do_previa}\n")
        log.write(f"- Total combinado (sem redundância): {num_ensg_unicos}\n\n")
        
        log.write("CORRESPONDÊNCIAS ENCONTRADAS (UNIPROT):\n")
        log.write(f"Ensembl Gene: {stats['gene']} de {stats['total_uniprot']} ({stats['gene']/stats['total_uniprot']*100:.2f}%)\n")
        log.write(f"Ensembl Transcript: {stats['transcrito']} de {stats['total_uniprot']} ({stats['transcrito']/stats['total_uniprot']*100:.2f}%)\n")
        log.write(f"Ensembl Protein: {stats['proteina']} de {stats['total_uniprot']} ({stats['proteina']/stats['total_uniprot']*100:.2f}%)\n")
        log.write(f"Pelo menos uma correspondência: {stats['com_correspondencia']} de {stats['total_uniprot']} ({stats['com_correspondencia']/stats['total_uniprot']*100:.2f}%)\n\n")
        
        log.write("ARQUIVOS GERADOS:\n")
        log.write(f"- {arquivo_combinado}: Todos os IDs processados\n")
        log.write(f"- uniprot_ensembl_found.txt: Apenas IDs com correspondências\n")
        log.write(f"- {arquivo_ensg_combinados}: IDs ENSG únicos combinados\n")
        log.write(f"- estatisticas_detalhadas.txt: Estatísticas completas\n")
        
        log.write("\nORIGEM DAS CORRELAÇÕES (APÓS APLICAÇÃO DA HIERARQUIA):\n")
        for origem, count in sorted(contagem_origem_final.items(), key=lambda x: HIERARQUIA[x[0]]):
            total_na_origem = totais_por_origem.get(origem, 0)
            if total_na_origem > 0:
                porcentagem_preservada = count / total_na_origem * 100
                log.write(f"- {origem}: {count} de {total_na_origem} ({porcentagem_preservada:.2f}% preservados)\n")
            else:
                log.write(f"- {origem}: {count} de 0 (N/A)\n")

    tempo_execucao = time.time() - inicio
    
    print("\n" + "="*70)
    print("PROCESSAMENTO CONCLUÍDO!")
    print("="*70)
    print(f"Arquivos processados: {len(arquivos_encontrados)}")
    print(f"IDs únicos processados: {len(dados_combinados)}")
    print(f"IDs UNIPROT processados: {stats['total_uniprot']}")
    print(f"\nFONTES DOS IDs ENSG:")
    print(f"- IDs do arquivo combinado: {num_do_combinado}")
    print(f"- IDs do arquivo All_Unique_ensembl.txt: {num_do_previa}")
    print(f"- IDs ENSG únicos combinados: {num_ensg_unicos}")
    print(f"\nCORRESPONDÊNCIAS ENSEMBL (UNIPROT):")
    print(f"- Ensembl Gene: {stats['gene']}/{stats['total_uniprot']} ({stats['gene']/stats['total_uniprot']*100:.2f}%)")
    print(f"- Ensembl Transcript: {stats['transcrito']}/{stats['total_uniprot']} ({stats['transcrito']/stats['total_uniprot']*100:.2f}%)")
    print(f"- Ensembl Protein: {stats['proteina']}/{stats['total_uniprot']} ({stats['proteina']/stats['total_uniprot']*100:.2f}%)")
    print(f"- Pelo menos uma correspondência: {stats['com_correspondencia']}/{stats['total_uniprot']} ({stats['com_correspondencia']/stats['total_uniprot']*100:.2f}%)")
    
    print(f"\nArquivos gerados:")
    print(f"- {arquivo_combinado}")
    print(f"- uniprot_ensembl_found.txt")
    print(f"- {arquivo_ensg_combinados}")
    print(f"- processamento.log")
    print(f"- estatisticas_detalhadas.txt")
    print(f"\nTempo de execução: {tempo_execucao:.2f} segundos")

if __name__ == "__main__":
    main()