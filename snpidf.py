from Bio.Align import PairwiseAligner  # Importação correta do PairwiseAligner

# Definindo exceções personalizadas
class NotATCG(Exception):
    pass

class Notdigit(Exception):
    pass

# Função para encontrar a posição de uma sequência por similaridade
def find_sequence_by_similarity(seq, subseq):
    # Cria a instância do PairwiseAligner
    aligner = PairwiseAligner()
    
    # Realiza alinhamento local (Smith-Waterman)
    alignments = aligner.align(seq, subseq)

    if alignments:
        # Seleciona o melhor alinhamento
        best_alignment = alignments[0]
        start = best_alignment.start + 1  # Ajuste para 1-indexação
        end = best_alignment.end
        print("Melhor alinhamento encontrado:")
        print(best_alignment)  # Exibe o alinhamento encontrado
        print(f"Subsequência encontrada na posição: {start} a {end}")
        return start
    else:
        print("Nenhuma correspondência encontrada.")
        return None

# Função principal para processar sequências e identificar SNPs
def process_sequences():
    try:
        # Solicita informações do usuário
        file01 = input("Digite o nome do primeiro arquivo (extensões válidas: .fasta, .fa, .nt): ").strip()
        file02 = input("Digite o nome do segundo arquivo (extensões válidas: .fasta, .fa, .nt): ").strip()

        # Verificação das extensões dos arquivos
        if not (file01.endswith(('.fasta', '.fa', '.nt')) and file02.endswith(('.fasta', '.fa', '.nt'))):
            raise ValueError("Os arquivos devem ter uma das seguintes extensões: .fasta, .fa, .nt.")

        # Lendo os arquivos
        def read_fasta(file):
            name = ""
            sequence = ""
            with open(file, "r") as f:
                for line in f:
                    line = line.strip()
                    if line.startswith(">"):
                        name = line
                    else:
                        sequence += line.upper()
            return name, sequence

        name01, seq01 = read_fasta(file01)
        name02, seq02 = read_fasta(file02)

        # Verificando se as sequências contêm apenas caracteres válidos
        valid_chars = set("ATCGN")
        if not set(seq01).issubset(valid_chars) or not set(seq02).issubset(valid_chars):
            raise NotATCG("Erro: Sequências contêm caracteres inválidos. Use apenas A, T, C, G, N.")

        # Encontrando posição inicial das sequências por similaridade
        print("Buscando posição inicial da sequência 02 na sequência 01...")
        start_pos = find_sequence_by_similarity(seq01, seq02)
        if start_pos is None:
            print("Não foi possível determinar a posição inicial.")
            return

        # Ajustando as sequências para o alinhamento
        posi01 = 1  # A sequência 01 começa no início
        posi02 = start_pos
        add01 = "x" * (posi01 - 1)
        add02 = "x" * (posi02 - 1)
        seq01 = add01 + seq01
        seq02 = add02 + seq02

        # Deixando as sequências com o mesmo comprimento
        len_diff = abs(len(seq01) - len(seq02))
        if len(seq01) < len(seq02):
            seq01 += "x" * len_diff
        else:
            seq02 += "x" * len_diff

        # Buscando SNPs
        print("Comparando sequências...")
        snps = []
        inter = 5  # Intervalo de nucleotídeos a verificar antes e depois do SNP
        for i in range(len(seq01)):
            if seq01[i] != "x" and seq02[i] != "x" and seq01[i] != seq02[i]:
                if (
                    i >= inter and i + inter < len(seq01) and
                    seq01[i - inter:i] == seq02[i - inter:i] and
                    seq01[i + 1:i + 1 + inter] == seq02[i + 1:i + 1 + inter]
                ):
                    snps.append((i + 1, seq01[i], seq02[i]))

        # Exibindo os resultados
        if snps:
            print("SNPs identificados:")
            for pos, base1, base2 in snps:
                print(f"Posição: {pos}, {base1} -> {base2}")
            print(f"Total de {len(snps)} SNPs identificados.")
        else:
            print("Nenhum SNP identificado.")

    except ValueError as ve:
        print(f"Erro: {ve}")
    except FileNotFoundError as fnfe:
        print(f"Erro: Arquivo não encontrado - {fnfe}")
    except NotATCG as natcg:
        print(natcg)
    except Exception as e:
        print(f"Erro inesperado: {e}")

# Executa a função principal
if __name__ == "__main__":
    process_sequences()
