import sys

def load_gpi_ids(gpi_file):
    with open(gpi_file) as f:
        return set(line.strip() for line in f if line.strip())

def modify_fasta_headers(fasta_file, gpi_ids, output_file):
    with open(fasta_file) as fin, open(output_file, 'w') as fout:
        for line in fin:
            if line.startswith(">"):
                parts = line.strip().split()
                protein_id = parts[0][1:]  # remove o ">"
                annotation = " ".join(parts[1:]) if len(parts) > 1 else ""

                if (protein_id in gpi_ids and 
                    "SP+TM" not in annotation and 
                    "TM" not in annotation and 
                    "GPI-Anchor" not in annotation):
                    fout.write(f">{protein_id} GPI-Anchor\n")
                else:
                    fout.write(line)
            else:
                fout.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Uso: python script.py <multifasta> <lista_ids_gpi.txt> <saida.fasta>")
        sys.exit(1)

    fasta_path = sys.argv[1]
    gpi_list_path = sys.argv[2]
    output_path = sys.argv[3]

    gpi_ids = load_gpi_ids(gpi_list_path)
    modify_fasta_headers(fasta_path, gpi_ids, output_path)
