import time
from datetime import datetime
from collections import defaultdict
import re

def extrair_ids_do_tsv(arquivo_tsv):
    """Extrai todos os IDs do arquivo TSV que têm correlação ENSEMBL"""
    ids_com_correlacao = set()
    
    with open(arquivo_tsv, 'r') as f:
        cabecalho = next(f)  # Pular cabeçalho
        for linha in f:
            partes = linha.strip().split('\t')
            if len(partes) >= 2 and partes[0] and partes[1]:  # ID e ENSEMBL(Gene) preenchidos
                ids_com_correlacao.add(partes[0])
    
    return ids_com_correlacao

def extrair_ids_do_fasta(arquivo_fasta):
    """Extrai todos os IDs do arquivo FASTA, identificando ENSG em qualquer posição do header e anotações"""
    ids_ensg = set()  # IDs que já contêm ENSG em qualquer parte do header
    todos_os_bancos = set()  # Todos os bancos de dados encontrados
    ids_por_banco = defaultdict(list)  # IDs agrupados por banco de dados
    anotacoes_por_id = {}  # Anotações por ID
    total_headers = 0  # Contador total de headers
    
    # Padrão regex para identificar códigos ENSG
    padrao_ensg = re.compile(r'ENSG\d{11}')
    
    # Lista de anotações conhecidas
    anotacoes_conhecidas = ['TM', 'SP+TM', 'GPI-Anchor']
    
    with open(arquivo_fasta, 'r') as f:
        for linha in f:
            if linha.startswith('>'):
                total_headers += 1  # Contar cada header
                
                # Pegar o header completo
                header = linha[1:].strip()
                
                # Extrair o ID inicial (primeira palavra)
                id_proteina = header.split()[0]
                
                # Identificar o banco de dados (segunda palavra)
                partes = header.split()
                banco_dados = "Desconhecido"
                if len(partes) > 1:
                    banco_dados = partes[1]
                
                # Adicionar à lista de todos os bancos
                todos_os_bancos.add(banco_dados)
                
                # Identificar anotações (última palavra, se for uma anotação conhecida)
                anotacao = None
                if len(partes) > 2:
                    ultima_palavra = partes[-1]
                    if ultima_palavra in anotacoes_conhecidas:
                        anotacao = ultima_palavra
                
                # Armazenar anotação para este ID
                anotacoes_por_id[id_proteina] = anotacao
                
                # Verificar se há algum código ENSG em qualquer parte do header
                if padrao_ensg.search(header):
                    ids_ensg.add(id_proteina)
                    # Mesmo que tenha ENSG, ainda assim registramos o banco
                    ids_por_banco[banco_dados].append((id_proteina, True))
                else:
                    # Se não tem ENSG, processar normalmente
                    ids_por_banco[banco_dados].append((id_proteina, False))
    
    return ids_ensg, ids_por_banco, todos_os_bancos, anotacoes_por_id, total_headers

def analisar_anotacoes_sem_correlacao(ids_sem_correlacao, anotacoes_por_id):
    """Analisa as anotações das proteínas sem correlação"""
    contagem_anotacoes = defaultdict(int)
    
    for id_proteina in ids_sem_correlacao:
        anotacao = anotacoes_por_id.get(id_proteina, None)
        if anotacao is None:
            contagem_anotacoes["Sem anotação"] += 1
        else:
            contagem_anotacoes[anotacao] += 1
    
    return contagem_anotacoes

def main():
    inicio = time.time()
    
    # Arquivos de entrada
    arquivo_tsv = "Correlation_All.tsv"
    arquivo_fasta = "All_Unique.fasta"
    
    # Extrair IDs
    print("Extraindo IDs do arquivo TSV...")
    ids_com_correlacao = extrair_ids_do_tsv(arquivo_tsv)
    print(f"Encontrados {len(ids_com_correlacao)} IDs com correlação ENSEMBL")
    
    print("Extraindo IDs do arquivo FASTA...")
    ids_ensg_fasta, ids_por_banco, todos_os_bancos, anotacoes_por_id, total_headers = extrair_ids_do_fasta(arquivo_fasta)
    print(f"Encontrados {len(ids_ensg_fasta)} IDs ENSG no arquivo FASTA")
    print(f"Total de headers no arquivo FASTA: {total_headers}")
    
    # Calcular estatísticas por banco de dados
    estatisticas_banco = {}
    
    for banco in todos_os_bancos:
        # Inicializar estatísticas para este banco
        total_banco = 0
        total_com_ensg = 0
        total_sem_ensg = 0
        ids_sem_correlacao = []
        
        # Processar todas as proteínas deste banco
        if banco in ids_por_banco:
            for id_proteina, tem_ensg in ids_por_banco[banco]:
                total_banco += 1
                if tem_ensg:
                    total_com_ensg += 1
                else:
                    total_sem_ensg += 1
                    # Verificar se não tem correlação no TSV
                    if id_proteina not in ids_com_correlacao:
                        ids_sem_correlacao.append(id_proteina)
        
        # Calcular porcentagens
        if total_sem_ensg > 0:
            porcentagem_sem_correlacao = (len(ids_sem_correlacao) / total_sem_ensg) * 100
        else:
            porcentagem_sem_correlacao = 0.0
        
        # Analisar anotações das proteínas sem correlação
        anotacoes_sem_correlacao = analisar_anotacoes_sem_correlacao(ids_sem_correlacao, anotacoes_por_id)
            
        estatisticas_banco[banco] = {
            'total': total_banco,
            'com_ensg': total_com_ensg,
            'sem_ensg': total_sem_ensg,
            'sem_correlacao': len(ids_sem_correlacao),
            'porcentagem_sem_correlacao': porcentagem_sem_correlacao,
            'ids_sem_correlacao': ids_sem_correlacao,
            'anotacoes_sem_correlacao': anotacoes_sem_correlacao
        }
    
    # Calcular totais gerais
    total_nao_ensg = sum([estat['sem_ensg'] for estat in estatisticas_banco.values()])
    total_sem_correlacao_geral = sum([estat['sem_correlacao'] for estat in estatisticas_banco.values()])
    total_fasta = total_headers  # Usar o contador direto de headers
    
    if total_nao_ensg > 0:
        porcentagem_sem_correlacao_geral = (total_sem_correlacao_geral / total_nao_ensg) * 100
    else:
        porcentagem_sem_correlacao_geral = 0
    
    # Verificar consistência das contagens
    total_por_bancos = sum([estat['total'] for estat in estatisticas_banco.values()])
    print(f"\nVerificação de consistência:")
    print(f"Total de headers no FASTA: {total_fasta}")
    print(f"Soma de proteínas por banco: {total_por_bancos}")
    print(f"Proteínas com ENSG no header: {len(ids_ensg_fasta)}")
    print(f"Proteínas sem ENSG no header: {total_nao_ensg}")
    
    if total_fasta != total_por_bancos:
        print(f"AVISO: Inconsistência na contagem! Total de headers ({total_fasta}) diferente da soma por bancos ({total_por_bancos})")
    
    # Gerar relatório
    data_processamento = datetime.now().strftime("%a %d %b %Y %H:%M:%S -03")
    
    with open('relatorio_correlacao_por_banco.txt', 'w') as relatorio:
        relatorio.write("RELATÓRIO DE CORRELAÇÃO ENSEMBL POR BANCO DE DADOS\n")
        relatorio.write("=" * 60 + "\n")
        relatorio.write(f"Data de processamento: {data_processamento}\n")
        relatorio.write(f"Arquivo de correlação: {arquivo_tsv}\n")
        relatorio.write(f"Arquivo FASTA: {arquivo_fasta}\n\n")
        
        relatorio.write("ESTATÍSTICAS GERAIS:\n")
        relatorio.write("-" * 25 + "\n")
        relatorio.write(f"Total de proteínas no FASTA: {total_fasta}\n")
        relatorio.write(f"Proteínas que já contêm ENSG: {len(ids_ensg_fasta)}\n")
        relatorio.write(f"Proteínas sem ENSG no header: {total_nao_ensg}\n")
        relatorio.write(f"Proteínas sem ENSG no header e SEM correlação: {total_sem_correlacao_geral} ({porcentagem_sem_correlacao_geral:.2f}%)\n\n")
        
        relatorio.write("ESTATÍSTICAS POR BANCO DE DADOS:\n")
        relatorio.write("-" * 35 + "\n")
        
        # Ordenar bancos por porcentagem de não correlação (decrescente)
        bancos_ordenados = sorted(estatisticas_banco.items(), 
                                 key=lambda x: x[1]['porcentagem_sem_correlacao'], 
                                 reverse=True)
        
        for banco, estat in bancos_ordenados:
            relatorio.write(f"{banco}:\n")
            relatorio.write(f"  Total: {estat['total']}\n")
            relatorio.write(f"  Já têm ENSG: {estat['com_ensg']}\n")
            relatorio.write(f"  Precisam de correlação: {estat['sem_ensg']}\n")
            relatorio.write(f"  Sem correlação: {estat['sem_correlacao']} ({estat['porcentagem_sem_correlacao']:.2f}% das que precisam)\n")
            
            # Adicionar informações de anotações para proteínas sem correlação
            if estat['sem_correlacao'] > 0:
                relatorio.write(f"  Distribuição de anotações (sem correlação):\n")
                for anotacao, count in estat['anotacoes_sem_correlacao'].items():
                    porcentagem_anotacao = (count / estat['sem_correlacao']) * 100
                    relatorio.write(f"    {anotacao}: {count} ({porcentagem_anotacao:.2f}%)\n")
            relatorio.write("\n")
        
        relatorio.write("DETALHES DOS IDs SEM CORRELAÇÃO POR BANCO:\n")
        relatorio.write("-" * 45 + "\n")
        
        for banco, estat in bancos_ordenados:
            if estat['sem_correlacao'] > 0:
                relatorio.write(f"\n{banco} ({estat['sem_correlacao']} IDs sem correlação):\n")
                for id_proteina in estat['ids_sem_correlacao']:
                    anotacao = anotacoes_por_id.get(id_proteina, "Sem anotação")
                    relatorio.write(f"  {id_proteina} - Anotação: {anotacao}\n")
    
    # Gerar arquivos separados por banco de dados
    for banco, estat in estatisticas_banco.items():
        if estat['sem_correlacao'] > 0:
            nome_arquivo = f"proteinas_sem_correlacao_{banco.replace(' ', '_')}.txt"
            with open(nome_arquivo, 'w') as f:
                for id_proteina in estat['ids_sem_correlacao']:
                    anotacao = anotacoes_por_id.get(id_proteina, "Sem anotação")
                    f.write(f"{id_proteina}\t{anotacao}\n")
    
    # Gerar relatório consolidado de anotações
    with open('relatorio_anotacoes_sem_correlacao.txt', 'w') as relatorio_anotacoes:
        relatorio_anotacoes.write("RELATÓRIO DE ANOTAÇÕES - PROTEÍNAS SEM CORRELAÇÃO ENSEMBL\n")
        relatorio_anotacoes.write("=" * 65 + "\n")
        relatorio_anotacoes.write(f"Data de processamento: {data_processamento}\n\n")
        
        # Consolidar anotações de todos os bancos
        anotacoes_consolidadas = defaultdict(int)
        for banco, estat in estatisticas_banco.items():
            for anotacao, count in estat['anotacoes_sem_correlacao'].items():
                anotacoes_consolidadas[anotacao] += count
        
        relatorio_anotacoes.write("DISTRIBUIÇÃO GERAL DE ANOTAÇÕES:\n")
        relatorio_anotacoes.write("-" * 35 + "\n")
        for anotacao, count in sorted(anotacoes_consolidadas.items(), key=lambda x: x[1], reverse=True):
            porcentagem = (count / total_sem_correlacao_geral) * 100
            relatorio_anotacoes.write(f"{anotacao}: {count} ({porcentagem:.2f}%)\n")
        
        relatorio_anotacoes.write("\nDISTRIBUIÇÃO POR BANCO DE DADOS:\n")
        relatorio_anotacoes.write("-" * 30 + "\n")
        for banco, estat in bancos_ordenados:
            if estat['sem_correlacao'] > 0:
                relatorio_anotacoes.write(f"\n{banco}:\n")
                for anotacao, count in estat['anotacoes_sem_correlacao'].items():
                    porcentagem = (count / estat['sem_correlacao']) * 100
                    relatorio_anotacoes.write(f"  {anotacao}: {count} ({porcentagem:.2f}%)\n")
    
    tempo_execucao = time.time() - inicio
    
    print("\n" + "="*70)
    print("RELATÓRIO DE CORRELAÇÃO ENSEMBL POR BANCO DE DADOS")
    print("="*70)
    print(f"Total de proteínas no arquivo FASTA: {total_fasta}")
    print(f"Proteínas que já contêm ENSG: {len(ids_ensg_fasta)}")
    print(f"Proteínas sem ENSG no header: {total_nao_ensg}")
    print(f"Proteínas sem ENSG no header e SEM correlação: {total_sem_correlacao_geral} ({porcentagem_sem_correlacao_geral:.2f}%)")
    print(f"\nESTATÍSTICAS POR BANCO DE DADOS:")
    print("-" * 35)
    
    # Ordenar bancos por porcentagem de não correlação (decrescente)
    bancos_ordenados = sorted(estatisticas_banco.items(), 
                             key=lambda x: x[1]['porcentagem_sem_correlacao'], 
                             reverse=True)
    
    for banco, estat in bancos_ordenados:
        print(f"{banco}:")
        print(f"  Total: {estat['total']}")
        print(f"  Já têm ENSG: {estat['com_ensg']}")
        print(f"  Precisam de correlação: {estat['sem_ensg']}")
        print(f"  Sem correlação: {estat['sem_correlacao']} ({estat['porcentagem_sem_correlacao']:.2f}% das que precisam)")
        
        # Adicionar informações de anotações para proteínas sem correlação
        if estat['sem_correlacao'] > 0:
            print(f"  Distribuição de anotações (sem correlação):")
            for anotacao, count in estat['anotacoes_sem_correlacao'].items():
                porcentagem_anotacao = (count / estat['sem_correlacao']) * 100
                print(f"    {anotacao}: {count} ({porcentagem_anotacao:.2f}%)")
        print()
    
    print(f"\nArquivos gerados:")
    print(f"- relatorio_correlacao_por_banco.txt (relatório completo)")
    print(f"- relatorio_anotacoes_sem_correlacao.txt (relatório de anotações)")
    for banco, estat in estatisticas_banco.items():
        if estat['sem_correlacao'] > 0:
            nome_arquivo = f"proteinas_sem_correlacao_{banco.replace(' ', '_')}.txt"
            print(f"- {nome_arquivo}")
    print(f"\nTempo de execução: {tempo_execucao:.2f} segundos")

if __name__ == "__main__":
    main()