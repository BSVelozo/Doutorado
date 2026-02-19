"""Este script combina múltiplos arquivos FASTA em um único FASTA, deduplicando proteínas pelo
identificador UniProt e agregando automaticamente informações de fonte e de origem.

O programa recebe uma lista arbitrária de FASTAs de entrada e gera um único FASTA de saída.
Cada arquivo deve estar em um dos dois formatos suportados:
(1) >ID
(2) >ID | source=fonteA,fonteB |

Para cada registro FASTA, o script extrai o UniProt ID, preserva apenas a primeira sequência
encontrada para cada ID e agrega:
- as fontes vindas dos campos source= dos headers (quando existirem);
- os nomes dos arquivos FASTA de origem (sem extensão).

O script alterna automaticamente entre dois formatos de saída:
- Modo simples (nenhum header contém source=):
  >UNIPROT_ID | source=arquivoA,arquivoB |
- Modo estendido (pelo menos um header contém source=):
  >UNIPROT_ID | source=fonteA,fonteB | origin=arquivoA,arquivoB |

A saída contém uma única entrada por UniProt ID, com sequências deduplicadas e headers
consistentes.

Uso:
python unificar_IdentificarOrigem_ProteinasInteractoras.py --i file1.fasta file2.fasta ... --o unico.fasta"""

import argparse
import os
import re
import sys

def strip_fasta_ext(filename: str) -> str:
    base = os.path.basename(filename)
    name, _ext = os.path.splitext(base)
    return name


def extract_uniprot_id(header: str) -> str:
    h = header[1:].strip()
    h = h.split("|", 1)[0].strip()
    h = h.split(None, 1)[0].strip()
    return h


def _normalize_tag(x: str) -> str:
    """
    Normaliza um item de source/origin:
    - remove espaços nas bordas
    - colapsa múltiplos espaços internos para um só
    """
    x = x.strip()
    x = re.sub(r"\s+", " ", x)
    return x


def extract_source_list(header: str):
    m = re.search(r"source=([^|]+)", header, flags=re.IGNORECASE)
    if not m:
        return []
    raw_items = m.group(1).split(",")
    items = []
    for it in raw_items:
        it = _normalize_tag(it)
        if it:
            items.append(it)
    return items


def parse_fasta(filepath: str):
    records = []
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        header = None
        seq_lines = []

        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_lines)))
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)

        if header is not None:
            records.append((header, "".join(seq_lines)))

    return records


def main():
    parser = argparse.ArgumentParser(
        description="Combina múltiplos FASTAs deduplicando por UniProt ID e agregando source/origin automaticamente."
    )
    parser.add_argument(
        "-i", "--inputs", nargs="+", required=True,
        help="Arquivos FASTA de entrada (todos no mesmo formato)."
    )
    parser.add_argument(
        "-o", "--output", required=True,
        help="Arquivo FASTA de saída."
    )

    args = parser.parse_args()

    proteins = {}
    sources = {}
    file_origins = {}
    extended_mode = False

    for fasta in args.inputs:
        file_tag = strip_fasta_ext(fasta)
        for header, sequence in parse_fasta(fasta):
            uniprot = extract_uniprot_id(header)
            src_list = extract_source_list(header)

            if src_list:
                extended_mode = True

            if uniprot not in proteins:
                proteins[uniprot] = sequence
                sources[uniprot] = set()
                file_origins[uniprot] = set()
            else:
                if proteins[uniprot] != sequence:
                    print(
                        f"[!] Aviso: ID {uniprot} aparece com sequências diferentes. Mantendo a primeira ocorrência.",
                        file=sys.stderr,
                    )

            # deduplicação robusta (TCDB vs "TCDB ")
            if src_list:
                sources[uniprot].update(src_list)

            # origem do arquivo (já vem sem .fasta pelo strip_fasta_ext)
            file_origins[uniprot].add(_normalize_tag(file_tag))

    with open(args.output, "w", encoding="utf-8") as out:
        for uniprot in sorted(proteins):
            seq = proteins[uniprot]

            if not extended_mode:
                origin_str = ",".join(sorted(file_origins[uniprot]))
                out.write(f">{uniprot} | source={origin_str} |\n{seq}\n")
            else:
                # também normaliza saída: sem espaços extras após vírgula
                src_str = ",".join(sorted(sources[uniprot]))
                origin_str = ",".join(sorted(file_origins[uniprot]))
                out.write(f">{uniprot} | source={src_str} | origin={origin_str} |\n{seq}\n")


if __name__ == "__main__":
    main()
