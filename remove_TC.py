import sys

def carregar_classes_removidas(arquivo_classes):
    with open(arquivo_classes, 'r') as f:
        return set(linha.strip() for linha in f if linha.strip())

def remover_sequencias_por_classe(fasta_in, classes_remover, fasta_out):
    with open(fasta_in, 'r') as fin, open(fasta_out, 'w') as fout:
        escrever = False
        header = ""

        for linha in fin:
            if linha.startswith(">"):
                header = linha.strip()
                # Extrai o numero da classe (por ex: 1.A.9.5.14)
                partes = header.split()
                if len(partes) > 1 and any(partes[1].startswith(cl) for cl in classes_remover):
                    escrever = False
                else:
                    escrever = True
                    fout.write(header + "\n")
            else:
                if escrever:
                    fout.write(linha)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Uso: python {sys.argv[0]} entrada.fasta classes_remover.txt saida.fasta")
        sys.exit(1)

    fasta_in = sys.argv[1]
    arquivo_classes = sys.argv[2]
    fasta_out = sys.argv[3]

    classes_remover = carregar_classes_removidas(arquivo_classes)
    remover_sequencias_por_classe(fasta_in, classes_remover, fasta_out)

    print(f"[INFO] Arquivo gerado: {fasta_out}")
