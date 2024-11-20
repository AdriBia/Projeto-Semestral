from Bio import SeqIO
from Bio.Align import PairwiseAligner

# Função para ler a sequência de um arquivo FASTA
def read_fasta(file):
    with open(file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            return str(record.seq)

# Função para realizar o alinhamento das sequências
def align_sequences(seq1, seq2):
    aligner = PairwiseAligner()
    aligner.mode = 'local'  # Alinhamento local
    aligner.match = 2       # Pontuação para correspondência
    aligner.mismatch = -1   # Penalidade por divergência
    aligner.open_gap_score = -0.5  # Penalidade para abertura de gap
    aligner.extend_gap_score = -0.1  # Penalidade para extensão de gap

    alignment = aligner.align(seq1, seq2)  # Alinhamento
    return alignment

# Função principal para processar as sequências e salvar os resultados
def process_sequences():
    try:
        # Solicita os arquivos de entrada
        file1 = input("Digite o nome do primeiro arquivo (extensões válidas: .fasta, .fa, .nt): ").strip()
        file2 = input("Digite o nome do segundo arquivo (extensões válidas: .fasta, .fa, .nt): ").strip()

        # Lê as sequências
        seq1 = read_fasta(file1)
        seq2 = read_fasta(file2)

        # Realiza o alinhamento
        print("Buscando posição inicial da sequência 02 na sequência 01...")
        alignment = align_sequences(seq1, seq2)

        if not alignment:
            print("Nenhum alinhamento encontrado.")
            return

        # Salva o alinhamento em um arquivo
        with open("alignment_output.txt", "w") as output_file:
            output_file.write("Melhor alinhamento encontrado:\n")
            for aln in alignment:
                output_file.write(str(aln) + "\n")

        print("Resultado do alinhamento salvo em alignment_output.txt")

    except ValueError as ve:
        print(f"Erro: {ve}")
    except FileNotFoundError as fnfe:
        print(f"Erro: Arquivo não encontrado - {fnfe}")
    except Exception as e:
        print(f"Erro inesperado: {e}")

# Executa a função principal
if __name__ == "__main__":
    process_sequences()
