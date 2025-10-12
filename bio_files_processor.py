def convert_multiline_fasta_to_oneline(input_fasta, output_fasta):
    with open(input_fasta, 'r') as raw_fastq, open(output_fasta, 'a'):
        read_line = input_fasta.readline().strip()
        if not read_line.startswith('>'):
            return 'У прочтениий нет названий'
        while input_fasta.readline().strip() != '':
            if line.startwith('>'):
                output_fasta.write(f'{read_line}\n')
                oneline = ''
                sequence_line = input_fasta.readline().strip()
                while not sequence_line.startswith('>'):
                    oneline += str(sequence_line)
                output_fasta.write(f'{oneline}\n')

def parse_blast_output():
    return


def select_genes_from_gbk_to_fasta():
    return