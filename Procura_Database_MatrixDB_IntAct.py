"""
Este script integra informações de um arquivo tabular de interações proteína-proteína com dois arquivos multifasta
para recuperar e organizar sequências de proteínas interagentes.

O script lê dois multifastas: um contendo todas as proteínas do Uniprot (queries) e outro contendo apenas as proteínas “alvo”.
Cada arquivo FASTA é convertido em um dicionário {header_completo: sequência}. Em seguida, é extraído de cada
header o identificador principal (primeiro token após '>'), normalizado para maiúsculas, criando índices do tipo
{UNIPROT_ID: header_completo} para permitir busca rápida de sequências por código.

O arquivo tabular (separado por TAB) é lido linha a linha e, nas colunas 1 e 2 (índices 0 e 1), são procurados
padrões do tipo 'uniprotkb:CODIGO'. Para cada par válido encontrado, é construído um mapeamento
bidirecional de interações no formato {CODIGO: {CODIGOS_INTERAGENTES}}.

Para cada proteína do primeiro multifasta, o script busca seus interatores nesse mapeamento e:
(1) opcionalmente remove interatores presentes em um arquivo de exclusão (um ID por linha),
(2) registra os pares removidos em um arquivo 'removidos.tsv',
(3) escreve um arquivo FASTA por proteína de entrada (<ID>.fasta) contendo as sequências dos interatores
    encontrados no segundo multifasta,
(4) acumula os interatores para gerar uma tabela de interações.

Ao final, o script gera:
- 'tabela_interacoes.tsv': uma tabela com uma linha por proteína de entrada e colunas Interactor_1..Interactor_N,
  preenchendo com campos vazios para manter largura uniforme;
- um arquivo FASTA de resumo (nome definido pelo usuário) contendo a concatenação de todos os FASTAs individuais;
- 'removidos.tsv' (se aplicável), listando os interatores excluídos por proteína de entrada.

Os IDs são tratados de forma case-insensitive (convertidos para maiúsculas), e avisos são impressos quando
uma proteína de entrada não possui correspondência na tabela de interações ou quando interatores não são
encontrados no segundo multifasta.

Uso:
python script.py <primeiro_multifasta> <arquivo_tabular> <segundo_multifasta> <pasta_saida> <arquivo_resumo> [arquivo_ids_excluir.txt]
"""

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
            col1 = parts[0]
            col2 = parts[1]

            match1 = re.search(r'uniprotkb:([\w\.-]+)', col1, re.IGNORECASE)
            match2 = re.search(r'uniprotkb:([\w\.-]+)', col2, re.IGNORECASE)

            if match1 and match2:
                code1 = match1.group(1).upper()
                code2 = match2.group(1).upper()
                mapping.setdefault(code1, set()).add(code2)
                mapping.setdefault(code2, set()).add(code1)
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
