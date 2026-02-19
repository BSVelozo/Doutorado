import sys

def carregar_ids_predicted_membrane(tsv_file):
    ids = set()
    with open(tsv_file, 'r', encoding='utf-8') as f:
        header = f.readline()
        header_cols = [col.strip().strip('"').lower() for col in header.strip().split('\t')]
        try:
            idx_uniprot = header_cols.index("uniprot")
            idx_protein_class = header_cols.index("protein class")
        except ValueError:
            print("[ERRO] Colunas 'Uniprot' ou 'Protein class' nao encontradas no cabecalho")
            return ids

        for line in f:
            cols = line.strip().split('\t')
            if len(cols) <= max(idx_uniprot, idx_protein_class):
                continue
            protein_class = cols[idx_protein_class]
            if "Predicted membrane proteins" in protein_class:
                uniprot_id = cols[idx_uniprot].strip()
                if uniprot_id:
                    ids.add(uniprot_id)
    print(f"[INFO] Total de IDs UniProt com 'Predicted membrane proteins': {len(ids)}")
    return ids

def filtrar_fasta_por_ids(fasta_in, fasta_out, ids_uniprot):
    with open(fasta_in, 'r', encoding='utf-8') as fin, open(fasta_out, 'w', encoding='utf-8') as fout:
        grava = False
        total_seq = 0
        seq_gravadas = 0

        for line in fin:
            if line.startswith(">"):
                total_seq += 1
                header_id = line[1:].strip().split()[0]
                grava = header_id in ids_uniprot
                if grava:
                    seq_gravadas += 1
                    fout.write(line)
            else:
                if grava:
                    fout.write(line)

        print(f"[INFO] Total de sequencias no multifasta: {total_seq}")
        print(f"[INFO] Sequencias gravadas no arquivo de saida: {seq_gravadas}")

def main():
    if len(sys.argv) != 4:
        print(f"Uso: python {sys.argv[0]} arquivo.tsv arquivo_entrada.fasta arquivo_saida.fasta")
        sys.exit(1)

    tsv_file = sys.argv[1]
    fasta_in = sys.argv[2]
    fasta_out = sys.argv[3]

    ids_uniprot = carregar_ids_predicted_membrane(tsv_file)
    if not ids_uniprot:
        print("[ERRO] Nenhum ID UniProt carregado, encerrando.")
        sys.exit(1)

    filtrar_fasta_por_ids(fasta_in, fasta_out, ids_uniprot)
    print("[INFO] Processamento concluido.")

if __name__ == "__main__":
    main()
