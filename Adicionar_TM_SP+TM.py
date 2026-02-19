import sys
import re

def carregar_sufixos(f3line_path):
    sufixos = {}
    with open(f3line_path) as f:
        for line in f:
            if line.startswith('>'):
                match = re.match(r'^>(\S+)\s*\|\s*(TM|SP\+TM)', line)
                if match:
                    uniprot_id, sufixo = match.groups()
                    sufixos[uniprot_id] = sufixo
    return sufixos

def processar_multifasta(fasta_path, sufixos, saida_path):
    with open(fasta_path) as f_in, open(saida_path, 'w') as f_out:
        escrever = False
        for line in f_in:
            if line.startswith('>'):
                seq_id = line[1:].strip().split()[0]
                if seq_id in sufixos:
                    f_out.write(f">{seq_id} {sufixos[seq_id]}\n")
                else:
                    f_out.write(line)
            else:
                f_out.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Uso: python adiciona_sufixos.py <arquivo.3line> <arquivo_multifasta> <saida.fasta>")
        sys.exit(1)

    f3line = sys.argv[1]
    multifasta = sys.argv[2]
    saida = sys.argv[3]

    sufixos_dict = carregar_sufixos(f3line)
    processar_multifasta(multifasta, sufixos_dict, saida)
