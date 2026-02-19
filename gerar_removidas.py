import re
import sys

def processar_tcdb(entrada, saida_tabela, saida_ids):
    removidas = []

    with open(entrada, 'r', encoding='utf-8') as infile:
        for linha in infile:
            linha = linha.strip()
            if "(RETIRAR)" in linha:
                match = re.match(r"\d+\s+([0-9A-Z\.\-]+)\s*-\s*(.+?)\s*\(RETIRAR\)", linha)
                if match:
                    tc_number = match.group(1).strip()
                    title = match.group(2).strip()
                    removidas.append((tc_number, title))

    with open(saida_tabela, 'w', encoding='utf-8') as tsv_out:
        tsv_out.write("TC Number\tTitle\tRemoved\n")
        for tc, title in removidas:
            tsv_out.write(f"{tc}\t{title}\tYes\n")

    with open(saida_ids, 'w', encoding='utf-8') as txt_out:
        for tc, _ in removidas:
            txt_out.write(f"{tc}\n")

    print(f"[INFO] Total de classes removidas: {len(removidas)}")
    print(f"[INFO] Arquivo de tabela salvo como: {saida_tabela}")
    print(f"[INFO] Lista de TC Numbers salva como: {saida_ids}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Uso: python {sys.argv[0]} entrada.txt removidas.tsv removidas.txt")
        sys.exit(1)

    entrada = sys.argv[1]
    saida_tabela = sys.argv[2]
    saida_ids = sys.argv[3]

    processar_tcdb(entrada, saida_tabela, saida_ids)
