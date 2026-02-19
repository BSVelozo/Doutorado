import sys

def parse_fasta(file_path):
    fasta = {}
    with open(file_path) as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header and seq_lines:
                    seq = "".join(seq_lines)
                    fasta[header] = seq
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if header and seq_lines:
            seq = "".join(seq_lines)
            fasta[header] = seq
    return fasta

def remove_duplicates(fasta_dict):
    seen = set()
    unique = {}
    for header, seq in fasta_dict.items():
        if seq not in seen:
            seen.add(seq)
            unique[header] = seq
    return unique

def write_fasta(fasta_dict, output_file):
    with open(output_file, 'w') as f:
        for header, seq in fasta_dict.items():
            f.write(f"{header}\n")
            # Divide em linhas de no maximo 60 caracteres (padrao FASTA)
            for i in range(0, len(seq), 60):
                f.write(f"{seq[i:i+60]}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Uso: python remover_sequencias_duplicadas.py <entrada.fasta> <saida.fasta>")
        sys.exit(1)

    entrada = sys.argv[1]
    saida = sys.argv[2]

    fasta = parse_fasta(entrada)
    unicos = remove_duplicates(fasta)
    write_fasta(unicos, saida)
