#!/usr/bin/env python3

import sys
from Bio import SeqIO

if len(sys.argv) != 4:
    print("Uso: python anotar_headers.py <arquivo1.fasta> <arquivo2_anotado.fasta> <saida.fasta>")
    sys.exit(1)

arquivo1 = sys.argv[1]
arquivo2 = sys.argv[2]
saida = sys.argv[3]

# Carrega os codigos com anotacao TM ou SP+TM
anotacoes_validas = {}
for registro in SeqIO.parse(arquivo2, "fasta"):
    partes = registro.description.strip().split()
    if len(partes) >= 2:
        codigo = partes[0]
        anotacao = partes[1]
        if anotacao in ["TM", "SP+TM"]:
            anotacoes_validas[codigo] = anotacao

# Processa o arquivo 1 e adiciona as anotacoes somente se forem TM ou SP+TM
registros_anotados = []
for registro in SeqIO.parse(arquivo1, "fasta"):
    codigo = registro.id
    if codigo in anotacoes_validas:
        registro.id = f"{codigo} {anotacoes_validas[codigo]}"
        registro.description = ""
    else:
        registro.id = codigo
        registro.description = ""
    registros_anotados.append(registro)

# Escreve o arquivo de saida
SeqIO.write(registros_anotados, saida, "fasta")
print(f"[?] Arquivo anotado salvo como: {saida}")
