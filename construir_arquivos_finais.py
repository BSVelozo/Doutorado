#Função: Cria os arquivos finais das proteínas extracelulares interactoras

#Descrição: O script lê os IDs presentes no arquivo CSV gerado pelo DeepLoc e comparacom os IDs do arquivo fasta com as proteínas interactoras para salvar
# apenas as sequências das proteínas interactoras extracelulares. Enquanto isso, ele também monta uma tabela que mostra as proteínas de membrana e suas interagentes.

#Uso:python construir_arquivos_finais.py --d [deeploc.csv] --p [tabela_interacoes.tsv] --f [proteins.fasta] --o [Nome_arquivo_saida]

import csv, argparse
from collections import defaultdict

p = argparse.ArgumentParser()
p.add_argument("--d", required=True)
p.add_argument("--p", required=True)
p.add_argument("--f", required=True)
p.add_argument("--o", default="filtered")
p.add_argument("--sep", default="\t")
a = p.parse_args()

# DeepLoc: proteínas extracelulares
extra = {
    r["Protein_ID"]
    for r in csv.DictReader(open(a.d))
    if float(r["Extracellular"]) > 0
}

# PPIs
ppi, keep = defaultdict(list), set()
with open(a.p) as f:
    next(f)
    for l in f:
        c = l.rstrip().split(a.sep)
        v = [i for i in c[1:] if i in extra]
        if v:
            ppi[c[0]] = v
            keep.update(v)

# Tabela PPI filtrada
m = max(map(len, ppi.values()))
with open(f"{a.o}_ppi.tsv", "w") as o:
    o.write(a.sep.join(["Input"] + [f"I{i+1}" for i in range(m)]) + "\n")
    for k, v in ppi.items():
        o.write(a.sep.join([k] + v + [""] * (m - len(v))) + "\n")

# FASTA filtrado (>ID | ...)
with open(a.f) as f, open(f"{a.o}.fasta", "w") as o:
    w = False
    for l in f:
        if l[0] == ">":
            w = l[1:].split("|", 1)[0].strip() in keep
            if w:
                o.write(l)
        elif w:
            o.write(l)

