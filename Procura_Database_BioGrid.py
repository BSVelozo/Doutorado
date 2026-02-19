"""Este script integra informações de um arquivo tabular de interações com dois arquivos multifasta para recuperar sequências de proteínas interagentes.

Fluxo geral:
1) Lê dois arquivos multifasta (um com todas as proteínas do Uniprot e outro com as proteínas “alvo”) e armazena cada sequência em um dicionário {header_completo: sequencia}.
2) Extrai de cada header FASTA o identificador principal (primeiro “token” logo após '>') e cria um índice {UNIPROT_ID} para permitir busca rápida de sequências pelo código.
3) Lê um arquivo tabular (separado por TAB) e, nas colunas 3 e 4 (índices 2 e 3), procura padrões do tipo 'uniprot/swiss-prot:CODIGO'(típico do BioGRID). Para cada par encontrado, monta um mapeamento bidirecional de interações {CODIGO: {CODIGOS_INTERAGENTES}}.
4) Para cada proteína do primeiro multifasta (input), busca seus interatores no mapeamento e:
- (Opcional) remove interatores listados em um arquivo de exclusão (um ID por linha).
- Registra quais interatores foram removidos (gerando removidos.tsv).
- Escreve um arquivo FASTA por proteína de entrada (<ID>.fasta) contendo as sequências dos interatores encontrados no segundo multifasta.
5) Gera uma tabela 'tabela_interacoes.tsv' com uma linha por proteína de entrada e colunas Interactor_1..Interactor_N (com preenchimento por vazio para manter a mesma largura).
6) Combina todos os FASTAs individuais gerados em um único arquivo de resumo (nome definido pelo usuário),salvo dentro da pasta de saída.

Observações:
- IDs são normalizados para maiúsculas para comparação.
- O script imprime avisos quando um ID de entrada não tem correspondência na tabela ou quando interatores não são encontrados no segundo multifasta.
- Uso: python script.py <primeiro_multifasta> <arquivo_tabular> <segundo_multifasta> <pasta_saida> <arquivo_resumo> [arquivo_ids_excluir.txt]"""

import re
from pathlib import Path

def parse_fasta(file_path):
    fasta = {}
    with open(file_path, 'r') as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    fasta[header] = "".join(seq_lines)
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if header:
            fasta[header] = "".join(seq_lines)
    return fasta

def extract_uniprot_ids(fasta_dict):
    uniprot_ids = {}
    for header in fasta_dict:
        match = re.match(r'>\s*([\w\.-]+)', header, re.IGNORECASE)
        if match:
            code = match.group(1).upper()
            uniprot_ids[code] = header
    return uniprot_ids

def build_mapping_from_tabular(tab_file):
    mapping = {}
    with open(tab_file, 'r') as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            col3 = parts[2]
            col4 = parts[3]

            match3 = re.search(r'uniprot/swiss-prot:([\w\.-]+)', col3, re.IGNORECASE)
            match4 = re.search(r'uniprot/swiss-prot:([\w\.-]+)', col4, re.IGNORECASE)

            if match3 and match4:
                code3 = match3.group(1).upper()
                code4 = match4.group(1).upper()
                mapping.setdefault(code3, set()).add(code4)
                mapping.setdefault(code4, set()).add(code3)
    return mapping

def carregar_lista_exclusao(caminho):
    with open(caminho, 'r') as f:
        return set(line.strip().upper() for line in f if line.strip())

def write_sequence(output_dir, query_code, match_codes, second_fasta_ids, second_fasta):
    written = False
    with open(output_dir / f"{query_code}.fasta", "a") as out:
        for code in match_codes:
            header = second_fasta_ids.get(code)
            if header:
                sequence = second_fasta[header]
                out.write(f"{header}\n{sequence}\n")
                written = True
    return written

def main(first_fasta_path, tabular_path, second_fasta_path, output_dir, summary_file, excluir_ids_path=None):
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    lista_exclusao = set()
    if excluir_ids_path:
        lista_exclusao = carregar_lista_exclusao(excluir_ids_path)

    first_fasta = parse_fasta(first_fasta_path)
    second_fasta = parse_fasta(second_fasta_path)

    first_ids = extract_uniprot_ids(first_fasta)
    second_fasta_ids = extract_uniprot_ids(second_fasta)

    mapping = build_mapping_from_tabular(tabular_path)

    interaction_table_lines = []
    removidos = []

    # Cabe�alho para a tabela de intera��es
    max_interactors = 0
    interacao_tmp = []

    for query_code in first_ids:
        if query_code in mapping:
            match_codes = mapping[query_code]
            match_codes_filtrados = {code for code in match_codes if code not in lista_exclusao}

            # Registra os removidos individualmente
            for excluido in match_codes:
                if excluido in lista_exclusao:
                    removidos.append([query_code, excluido])

            if not match_codes_filtrados:
                continue

            found = write_sequence(output_dir, query_code, match_codes_filtrados, second_fasta_ids, second_fasta)

            if not found:
                print(f"[!] Codigo(s) {match_codes_filtrados} nao encontrados no segundo multifasta.")

            interacao_tmp.append([query_code] + sorted(match_codes_filtrados))
            max_interactors = max(max_interactors, len(match_codes_filtrados))
        else:
            print(f"[!] Codigo {query_code} nao tem correspondencia na tabela.")

    # Escreve a tabela de intera��es com cabe�alho
    with open(output_dir / "tabela_interacoes.tsv", "w") as tsv:
        header = ["Input"] + [f"Interactor_{i+1}" for i in range(max_interactors)]
        tsv.write("\t".join(header) + "\n")
        for linha in interacao_tmp:
            tsv.write("\t".join(linha + [""] * (max_interactors - len(linha[1:]))) + "\n")

    # Combina os FASTAs individuais (evita incluir o próprio resumo)
    summary_path = output_dir / summary_file
    with open(summary_path, "w") as summary:
        for fasta_file in sorted(output_dir.glob("*.fasta")):
            if fasta_file.name == summary_file:
                continue  # Ignora o arquivo de resumo
            with open(fasta_file, "r") as infile:
                summary.write(infile.read())

    # Escreve relat�rio de exclus�es com input + identificadores exclu�dos
    if removidos:
        with open(output_dir / "removidos.tsv", "w") as out:
            out.write("Input\tProteina_excluida\n")
            for linha in removidos:
                out.write("\t".join(linha) + "\n")

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 6 or len(sys.argv) > 7:
        print("Uso: python script.py <primeiro_multifasta> <arquivo_tabular> <segundo_multifasta> <pasta_saida> <arquivo_resumo> [arquivo_ids_excluir.txt]")
        sys.exit(1)

    excluir_path = sys.argv[6] if len(sys.argv) == 7 else None
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], excluir_path)

