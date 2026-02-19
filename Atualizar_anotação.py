import sys

def carregar_anotacoes(arquivo_anotacoes):
    anotacoes = {}
    with open(arquivo_anotacoes, 'r') as f:
        for linha in f:
            if linha.startswith('>'):
                partes = linha.strip().split('|')
                if len(partes) == 2:
                    id_proteina = partes[0].strip().lstrip('>')
                    anotacao = partes[1].strip()
                    anotacoes[id_proteina] = anotacao
    return anotacoes

def anotar_multifasta(arquivo_multifasta, anotacoes, arquivo_saida):
    with open(arquivo_multifasta, 'r') as entrada, open(arquivo_saida, 'w') as saida:
        for linha in entrada:
            if linha.startswith('>'):
                id_proteina = linha.strip().lstrip('>')
                if id_proteina in anotacoes:
                    saida.write(f">{id_proteina} {anotacoes[id_proteina]}\n")
                else:
                    saida.write(f">{id_proteina}\n")
            else:
                saida.write(linha)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Uso: python script.py <arquivo_multifasta> <arquivo_anotacoes> <arquivo_saida>")
        sys.exit(1)

    arquivo_multifasta = sys.argv[1]
    arquivo_anotacoes = sys.argv[2]
    arquivo_saida = sys.argv[3]

    anotacoes = carregar_anotacoes(arquivo_anotacoes)
    anotar_multifasta(arquivo_multifasta, anotacoes, arquivo_saida)
