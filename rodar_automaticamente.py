#Função: Rodar scripts automáticamente de forma remota

#Descrição: Esse script permite com que o usuário rode outros scripts automáticamente de acordo com a sequência que ele colocar na hora de chamá-lo. Infelizmente,
# esta versão aceita apenas 3 scripts de uma vez e você precisa colocar todos os elementos de execução mencionando todo o caminho um de cada vez.

#Uso: Coloque o caminho de cada elemento de execução dos scripts desejados nas listas.

import subprocess
import sys

def test(lista):
    for i, caminho in enumerate(lista, start=1):
        print(f'Execução {i} iniciada')
        cmd = ["python3"] + caminho
        try:
            # check=True -> levanta exceção se o comando retornar erro
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"[ERRO] Script {caminho[0]} retornou erro (código {e.returncode}).")
            print("Saída de erro:")
            if e.stderr:
                print(e.stderr.decode(errors="ignore"))
            continue
        except FileNotFoundError:
            print(f"[ERRO] Arquivo '{caminho[0]}' não encontrado.")
            continue
        except KeyboardInterrupt:
            print("\n[INTERRUPÇÃO] Execução cancelada pelo usuário.")
            sys.exit(1)
        except Exception as e:
            print(f"[ERRO DESCONHECIDO] {e}")
            break
        else:
            print(f"--- Execução {i} terminada com sucesso ---\n")

if __name__ == "__main__":
    caminho1 = [
        "/home/single1/Doutorado_Bernardo/3.3_Identificaçao_proteinas_extracelulares_interactoras_com_proteinas_transmembrana/scripts/Main_script/Procura_Database_V6.py",
        "/home/single1/Doutorado_Bernardo/3.2_Identificaçao_proteinas_membrana/4_Results_Final/TCDB/TCDB_Final.fasta",
        "/home/single1/Doutorado_Bernardo/3.3_Identificaçao_proteinas_extracelulares_interactoras_com_proteinas_transmembrana/Bancos_de_interacao/IntAct/psimitab/Final_IntAct.txt",
	"/home/single1/Doutorado_Bernardo/All_Proteins_UNIPROT/UNIPROT_All_Clean.fasta",
	"TCDBxIntAct",
	"resumo_TCDB_IntAct.fasta",
	"/home/single1/Doutorado_Bernardo/3.3_Identificaçao_proteinas_extracelulares_interactoras_com_proteinas_transmembrana/Headers_UNIQUE.txt"
    ]

    caminho2 = [
        "/home/single1/Doutorado_Bernardo/3.3_Identificaçao_proteinas_extracelulares_interactoras_com_proteinas_transmembrana/scripts/Main_script/Procura_Database_V6_1.py",
        "/home/single1/Doutorado_Bernardo/3.2_Identificaçao_proteinas_membrana/4_Results_Final/TCDB/TCDB_Final.fasta",
        "/home/single1/Doutorado_Bernardo/3.3_Identificaçao_proteinas_extracelulares_interactoras_com_proteinas_transmembrana/Bancos_de_interacao/BioGRID/BIOGRID-ALL-4.4.243.mitab/Final_BIOGRID.txt",
        "/home/single1/Doutorado_Bernardo/All_Proteins_UNIPROT/UNIPROT_All_Clean.fasta",
        "TCDBxBioGRID",
        "resumo_TCDB_BioGRID.fasta",
        "/home/single1/Doutorado_Bernardo/3.3_Identificaçao_proteinas_extracelulares_interactoras_com_proteinas_transmembrana/Headers_UNIQUE.txt"
    ]

    caminho3 = [
        "/home/single1/Doutorado_Bernardo/3.3_Identificaçao_proteinas_extracelulares_interactoras_com_proteinas_transmembrana/scripts/Main_script/Procura_Database_V6.py",
        "/home/single1/Doutorado_Bernardo/3.2_Identificaçao_proteinas_membrana/4_Results_Final/TCDB/TCDB_Final.fasta",
        "/home/single1/Doutorado_Bernardo/3.3_Identificaçao_proteinas_extracelulares_interactoras_com_proteinas_transmembrana/Bancos_de_interacao/MatrixDB/matrixdb_all/Final_MatrixDB.txt",
        "/home/single1/Doutorado_Bernardo/All_Proteins_UNIPROT/UNIPROT_All_Clean.fasta",
        "TCDBxMatrixDB",
        "resumo_TCDB_MatrixDB.fasta",
        "/home/single1/Doutorado_Bernardo/3.3_Identificaçao_proteinas_extracelulares_interactoras_com_proteinas_transmembrana/Headers_UNIQUE.txt"
    ]

    paths_list = [caminho1, caminho2, caminho3]

    test(paths_list)
