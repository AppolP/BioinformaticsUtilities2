def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str) -> None:
    """
    Converts multiline sequences of DNA/RNA/protein given as a path to file.fasta into output.fasta file with glued lines into one

    Arguments:
    input_fasta: str, pathway to fasta file to convert
    output_fasta: str, pathway to fasta file to save converted sequences in

    Raises an error in case of unnamed sequences. Notice: unnamed sequences (except the first one) will be converted into one line
    by the given name no matter if there are two different sequences in a row.
    """

    with open(input_fasta, "r") as raw_fasta, open(
        output_fasta, "w"
    ) as multiline_fasta:
        read_line = raw_fasta.readline().strip()
        if not read_line.startswith(">"):
            print("Reads have no name")
            return None

        while read_line:
            if read_line.startswith(">"):
                multiline_fasta.write(f"{read_line}\n")
                multiline_sequences = []
                sequence_line = raw_fasta.readline().strip()
                flag = True
                while flag:
                    multiline_sequences.append(sequence_line)
                    sequence_line = raw_fasta.readline().strip()
                    if sequence_line.startswith(">") or not sequence_line:
                        flag = False
                multiline_fasta.write(f'{"".join(multiline_sequences)}\n')
                read_line = sequence_line
            else:
                read_line = raw_fasta.readline().strip()


def parse_blast_output(input_file: str, output_file: str) -> None:
    """
    Reads a file within the pathway user typed and saves sorted list of Description first lines for each Query

    Arguments:
    input_file: str, pathway to the txt file
    output_file:  str, pathway to save parsed file
    """

    with open(input_file, "r") as infile, open(output_file, "w") as res_file:
        line = infile.readline()
        proteins = []
        while line:
            if line.startswith("Description"):
                line = infile.readline().rstrip()
                proteins.append(line[:66])
            line = infile.readline()
        res_file.write("\n".join(sorted(proteins)))


def select_genes_from_gbk_to_fasta():
    return
