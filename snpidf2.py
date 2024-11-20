from Bio.Align import PairwiseAligner

# Função para encontrar a posição de uma sequência por similaridade
def find_sequence_by_similarity(seq, subseq):
    # Cria a instância do PairwiseAligner
    aligner = PairwiseAligner()
    
    # Ajustando os parâmetros de pontuação para evitar alinhamentos excessivos
    aligner.mode = 'local'  # Alinhamento local
    aligner.match = 2  # Pontuação para correspondência
    aligner.mismatch = -1  # Penalidade por divergência
    aligner.open_gap_score = -0.5  # Penalidade para abertura de gap
    aligner.extend_gap_score = -0.1  # Penalidade para extensão de gap

    # Realiza o alinhamento
    alignments = aligner.align(seq, subseq)

    if alignments:
        # Seleciona o melhor alinhamento
        best_alignment = alignments[0]
        aligned_seq1 = best_alignment.aligned[0]  # Alinhamento para a sequência 1
        aligned_seq2 = best_alignment.aligned[1]  # Alinhamento para a sequência 2
        
        # A posição de início no alinhamento pode ser extraída do índice de correspondência
        start_pos = aligned_seq1[0] + 1  # Ajuste para 1-indexação
        end_pos = aligned_seq1[-1] + 1  # Ajuste para 1-indexação

        return best_alignment, start_pos, end_pos
    else:
        print("Nenhuma correspondência encontrada.")
        return None, None, None

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
            raise ValueError("Erro: Sequências contêm caracteres inválidos. Use apenas A, T, C, G, N.")

        # Encontrando posição inicial das sequências por similaridade
        print("Buscando posição inicial da sequência 02 na sequência 01...")
        best_alignment, start_pos, end_pos = find_sequence_by_similarity(seq01, seq02)
        if best_alignment is None:
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

        # Solicitando o nome do arquivo de saída
        output_file = input("Digite o nome do arquivo para salvar o resultado (ex: resultado.txt): ").strip()

        # Salvando o alinhamento em um arquivo de saída
        alignment_file = "alignment_output.txt"
        with open(alignment_file, "w") as f:
            f.write("Melhor Alinhamento encontrado:\n")
            f.write(str(best_alignment))
            f.write(f"\nSubsequência encontrada na posição: {start_pos} a {end_pos}\n")

        # Salvando os resultados no arquivo de SNPs
        with open(output_file, "w") as f:
            if snps:
                f.write("SNPs identificados:\n")
                for pos, base1, base2 in snps:
                    f.write(f"Posição: {pos}, {base1} -> {base2}\n")
                f.write(f"Total de {len(snps)} SNPs identificados.\n")
            else:
                f.write("Nenhum SNP identificado.\n")

        print(f"Resultado do alinhamento salvo em {alignment_file}")
        print(f"Resultado dos SNPs salvo em {output_file}")

    except ValueError as ve:
        print(f"Erro: {ve}")
    except FileNotFoundError as fnfe:
        print(f"Erro: Arquivo não encontrado - {fnfe}")
    except Exception as e:
        print(f"Erro inesperado: {e}")

# Executa a função principal
if __name__ == "__main__":
    process_sequences()
