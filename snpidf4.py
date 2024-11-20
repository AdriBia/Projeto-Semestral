from Bio import Align
from Bio import SeqIO

# Solicita os nomes dos arquivos
seq_file1 = input("Digite o nome do primeiro arquivo (extensões válidas: .fasta, .fa, .nt): ")
seq_file2 = input("Digite o nome do segundo arquivo (extensões válidas: .fasta, .fa, .nt): ")

# Solicita o valor mínimo de correspondência
min_length = int(input("Digite o valor mínimo de correspondência (em pares de bases): "))

# Lê as sequências dos arquivos
seq1 = str(SeqIO.read(seq_file1, "fasta").seq)
seq2 = str(SeqIO.read(seq_file2, "fasta").seq)

# Cria um objeto PairwiseAligner
aligner = Align.PairwiseAligner()

# Alinha as sequências
alignments = aligner.align(seq1, seq2)

# Filtra os alinhamentos que têm pelo menos o comprimento mínimo de correspondência
valid_alignments = []
for alignment in alignments:
    # Verifica o comprimento da correspondência (número de bases iguais)
    match_length = sum(1 for a, b in zip(alignment.target, alignment.query) if a == b)
    if match_length >= min_length:
        valid_alignments.append(alignment)

# Exibe os alinhamentos válidos
if valid_alignments:
    print(f"{len(valid_alignments)} alinhamentos válidos encontrados.")
    for alignment in valid_alignments:
        print(f"Alinhamento: \n{alignment}")
else:
    print("Nenhum alinhamento válido encontrado.")
