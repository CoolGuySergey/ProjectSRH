GENE_BLOCK_START = "BEGIN SETS;"
BLOCK_END = "end;"

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--input_nexus', required=True, type=str)
    parser.add_argument('--output_partitioning', required=True, type=str)
    parser.add_argument('--tool', required=True, type=str, choices=('raxml', 'iq'))
    parser.add_argument('--level', required=False, type=str, choices=('gene', 'strand', 'pos'), default='strand')
    parser.add_argument('--reverse_strand', required=False, type=str, nargs='*')
    params = parser.parse_args()

    gene2pos = {}
    with open(params.input_nexus, 'r') as f:
        started = False
        for line in f.readlines():
            if GENE_BLOCK_START in line.upper():
                started = True
                continue
            if started and BLOCK_END in line.lower():
                break
            if not started:
                continue
            line = line.replace('charset', '').strip(' \t\n;')
            gene, positions = line.split('=')
            gene = gene.strip(' \t')
            start, end = [int(_) for _ in positions.strip(' \t').split('-')]
            gene2pos[gene] = start, end

    model = 'GTR+G6+FO+IO' if 'raxml' == params.tool else 'DNA'

    with open(params.output_partitioning, 'w+') as f:
        if 'gene' == params.level:
            for gene, (start_pos, stop_pos) in gene2pos.items():
                f.write("{model}, {gene}_12 = {start_1}-{stop}\\3, {start_2}-{stop}\\3\n"
                        "{model}, {gene}_3 = {start_3}-{stop}\\3\n"
                        .format(gene=gene, start_1=start_pos, start_2=start_pos + 1,
                                start_3=start_pos + 2, stop=stop_pos, model=model))
        elif 'strand' == params.level and params.reverse_strand:
            forward_str_12 = []
            forward_str_3 = []
            reverse_str_12 = []
            reverse_str_3 = []
            print(params.reverse_strand)
            print(gene2pos.keys())
            for gene, (start_pos, stop_pos) in gene2pos.items():
                if gene.lower() in params.reverse_strand:
                    reverse_str_12.append('{start_1}-{stop}\\3, {start_2}-{stop}\\3'.format(
                        start_1=start_pos, start_2=start_pos + 1, stop=stop_pos))
                    reverse_str_3.append('{start_3}-{stop}\\3'.format(
                        start_3=start_pos + 2, stop=stop_pos))
                else:
                    forward_str_12.append('{start_1}-{stop}\\3, {start_2}-{stop}\\3'.format(
                        start_1=start_pos, start_2=start_pos + 1, stop=stop_pos))
                    forward_str_3.append('{start_3}-{stop}\\3'.format(
                        start_3=start_pos + 2, stop=stop_pos))
            f.write("{model}, forward_12 = {forward_12}\n"
                    "{model}, forward_3 = {forward_3}\n"
                    "{model}, reverse_12 = {reverse_12}\n"
                    "{model}, reverse_3 = {reverse_3}\n"
                    .format(model=model,
                            forward_12=', '.join(forward_str_12), forward_3=', '.join(forward_str_3),
                            reverse_12=', '.join(reverse_str_12), reverse_3=', '.join(reverse_str_3)))
        else:
            str_12 = []
            str_3 = []
            for gene, (start_pos, stop_pos) in gene2pos.items():
                str_12.append('{start_1}-{stop}\\3, {start_2}-{stop}\\3'.format(
                    start_1=start_pos, start_2=start_pos + 1, stop=stop_pos))
                str_3.append('{start_3}-{stop}\\3'.format(
                    start_3=start_pos + 2, stop=stop_pos))
            f.write("{model}, pos_12 = {s_12}\n"
                    "{model}, pos_3 = {s_3}\n"
                    .format(model=model, s_12=', '.join(str_12), s_3=', '.join(str_3)))
