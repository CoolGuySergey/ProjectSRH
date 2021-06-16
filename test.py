# RenameFastaToSource

from Bio import SeqIO
import argparse

def run(args):

    PathToFasta = args.fasta
    PathToGB = args.gb
    
    GenBank, CommonName = [], []
    for seq_record in SeqIO.parse(PathToGB, "genbank"):
        GenBank.append(">" + seq_record.name)
        CommonName.append(seq_record.annotations["source"][14:])

    NameDict = dict(zip(GenBank, CommonName))
    NameDict['>NC_042922'] = "Himaloaesalus gaoligongshanus"
    NameDict['>NC_042614'] = "Sinodendron rugosum"

    fasta = open(PathToFasta)
    newfasta= open(f"{PathToFasta}_renamed", "a")

    for line in fasta:
        if line.startswith('>'):
            line = line.strip('\n')
            newname= NameDict[line]
            newfasta.write('>' + newname + '\n')
        else:
            newfasta.write(line)

    fasta.close()
    newfasta.close()

def main():
    parser.add_argument("-fasta", help="Relative path of input fasta", required=True, dest="fasta", type=str)
    parser.add_argument("-gb", help="Relative path of GB reference", required=True, dest="gb", type=str)
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
