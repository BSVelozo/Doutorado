import sys

def parse_fasta(file_path):
    sequences = {}
    with open(file_path, 'r') as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header and seq_lines:
                    seq = ''.join(seq_lines)
                    sequences[header] = seq
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if header and seq_lines:
            seq = ''.join(seq_lines)
            sequences[header] = seq
    return sequences

def normalize_seq(seq):
    return len(seq), ''.join(sorted(set(seq)))  # tamanho e diversidade

def main(fasta1_path, fasta2_path, output_fasta, report_path):
    fasta1 = parse_fasta(fasta1_path)
    fasta2 = parse_fasta(fasta2_path)

    # Indexa as sequencias do segundo FASTA com base em (tamanho, diversidade)
    seq_map = {}
    for header2, seq2 in fasta2.items():
        key = normalize_seq(seq2)
        seq_map[key] = header2  # Assume que nao ha duplicatas com mesmo perfil

    with open(output_fasta, 'w') as out_fasta, open(report_path, 'w') as out_report:
        out_report.write("UniProt_ID\tOriginal_ENSP_Header\n")
        for header1, seq1 in fasta1.items():
            key = normalize_seq(seq1)
            if key in seq_map:
                uniprot_id = seq_map[key].lstrip('>')
                new_header = f">{uniprot_id} {header1.lstrip('>')}"
                out_fasta.write(f"{new_header}\n")
                out_report.write(f"{uniprot_id}\t{header1.lstrip('>')}\n")
            else:
                out_fasta.write(f"{header1}\n")
            # Escreve a sequencia em blocos de 60 caracteres
            for i in range(0, len(seq1), 60):
                out_fasta.write(f"{seq1[i:i+60]}\n")

    print(f"[?] FASTA anotado salvo em: {output_fasta}")
    print(f"[?] Relatorio salvo em: {report_path}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Uso: python anotar_com_uniprot.py <fasta1> <fasta2> <saida_fasta> <saida_relatorio>")
        sys.exit(1)
    
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
