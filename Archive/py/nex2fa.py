import re


MATRIX_START = 'MATRIX'
MATRIX_STOP = ';'
if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input_nexus', required=True, type=str)
    parser.add_argument('-o', '--output_fasta', required=True, type=str)
    params = parser.parse_args()

    with open(params.input_nexus, 'r') as f:
        with open(params.output_fasta, 'w+') as fa:
            started = False
            for line in f.readlines():
                if MATRIX_START in line:
                    started = True
                    continue
                if started and MATRIX_STOP in line:
                    break
                if not started:
                    continue
                line = line.strip(' \t\n')
                if line:
                    name_seq = re.compile('\s+').split(line)
                    fa.write('>{}\n{}\n'.format(name_seq[0].replace('|', '_').replace(':', '-'), name_seq[-1]))
