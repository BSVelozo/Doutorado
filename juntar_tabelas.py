#Função: Unir as tabelas de interações proteicas de interesse em uma só

#Descrição: A partir dos arquivos tabulares de input, o script cria uma tabela única contendo a junção das informações sem redundância.

#Uso: python juntar_tabelas.py arquivo1.txt arquivo2.txt [arquivo3.txt ...]


import pandas as pd
import sys

# ==== VERIFICAÇÃO DO INPUT ====
if len(sys.argv) < 3:
    print("Uso: python juntar_tabelas.py arquivo1.txt arquivo2.txt [arquivo3.txt ...]")
    sys.exit(1)

arquivos = sys.argv[1:]  # caminhos passados na linha de comando
sep = "\t"  # separador TAB
saida = "interacoes_unificadas.txt"

# ==== Dicionário para guardar todas as interações ====
interacoes = {}

# ==== LER E COMBINAR ====
for arquivo in arquivos:
    df = pd.read_csv(arquivo, sep=sep, dtype=str).fillna("")
    
    for _, linha in df.iterrows():
        prot_input = linha.iloc[0]
        interatores = [x for x in linha.iloc[1:] if x != ""]
        
        if prot_input not in interacoes:
            interacoes[prot_input] = set()
        
        interacoes[prot_input].update(interatores)

# ==== CRIAR DATAFRAME FINAL ====
max_interatores = max(len(v) for v in interacoes.values())

dados_finais = []
for prot_input, inters in sorted(interacoes.items()):
    linha = [prot_input] + sorted(inters)
    linha += [""] * (max_interatores - len(inters))
    dados_finais.append(linha)

colunas = ["Protein_Input"] + [f"Interactor_{i+1}" for i in range(max_interatores)]
df_final = pd.DataFrame(dados_finais, columns=colunas)

# ==== SALVAR ====
df_final.to_csv(saida, sep=sep, index=False)
print(f"✅ Arquivo final salvo em: {saida}")
