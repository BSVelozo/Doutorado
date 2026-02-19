#!/usr/bin/env python3
import sys
import re

def parse_ids(ids_file):
    with open(ids_file) as f:
        return set(line.strip() for line in f if line.strip())

def extract_fasta_sequences(fasta_file, ids_set, output_file):
    with open(fasta_file) as fasta, open(output_file, 'w') as out:
        write_flag = False
        header_id = None
        for line in fasta:
            line = line.strip()
            if line.startswith(">"):
                match = re.match(r'^>(ENSP[0-9]+\.[0-9]+)', line)
                if match:
                    header_id = match.group(1)
                    if header_id in ids_set:
                        write_flag = True
                        out.write(f">{header_id}\n")
                    else:
                        write_flag = False
                else:
                    write_flag = False
            elif write_flag:
                out.write(line + '\n')

def main():
    if len(sys.argv) != 4:
        print("Uso: python extrai_sequencias.py <arquivo_ids.txt> <arquivo_multifasta.fasta> <saida.fasta>")
        sys.exit(1)

    ids_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_file = sys.argv[3]

    ids_set = parse_ids(ids_file)
    extract_fasta_sequences(fasta_file, ids_set, output_file)
    print(f"[?] Sequencias extraidas para '{output_file}'")

if __name__ == "__main__":
    main()
