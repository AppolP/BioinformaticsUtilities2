with open('example_blast_results.txt', 'r') as infile, open('result.txt', 'w') as res_file:
    line = infile.readline()
    proteins = []
    while line:
        if line.startswith('Description'):
            line = infile.readline().rstrip()
            proteins.append(line[:66])
        line = infile.readline()
    res_file.write('\n'.join(sorted(proteins)))        
