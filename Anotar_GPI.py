import sys

def carregar_ids(lista_path):
    with open(lista_path) as f:
        return set(line.strip() for line in f if line.strip())

def processar_fasta(fasta_path, ids_gpi, output_path):
    with open(fasta_path) as fasta_in, open(output_path, 'w') as fasta_out:
        for line in fasta_in:
            if line.startswith('>'):
                header = line.strip()
                # Extrai o ID (assume que esta logo apos o '>' ate o primeiro espaco ou fim de linha)
                id_base = header[1:].split()[0]
                # Se o ID esta na lista e nao contem TM ou SP+TM, adiciona GPI-Anchor
                if id_base in ids_gpi and not header.endswith('TM') and not header.endswith('SP+TM'):
                    header += ' GPI-Anchor'
                fasta_out.write(header + '\n')
            else:
                fasta_out.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Uso: python anotar_gpi.py <fasta_entrada> <lista_ids.txt> <saida.fasta>")
        sys.exit(1)

    fasta_path = sys.argv[1]
    lista_path = sys.argv[2]
    output_path = sys.argv[3]

    ids_gpi = carregar_ids(lista_path)
    processar_fasta(fasta_path, ids_gpi, output_path)
