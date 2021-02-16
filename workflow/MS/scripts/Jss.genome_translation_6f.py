#!/usr/bin/env python3
import argparse, sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def frame_translate_6(r_id, dna, aa_cutoff, writeFile):
    orf_numer = 0
    
    ### sense strand
    for frame in range(3):
        newrna = dna[frame:]
        newrna = newrna[:int(len(newrna)/3)*3]
        end = frame

        cds = newrna.translate()
        block_list = cds.split('*')
        for block in block_list[:-1]:
            start = end
            end += len(block)*3+3
            Pm = block.find('M')
            start += Pm*3
            block = block[Pm:]
            if len(block) >= aa_cutoff: 
                orf_numer += 1
                proid = r_id + ":" + str(start) + "-" + str(end) + ":+" + ":dnaORF_%s"%orf_numer 
                proseq = block
                orf = SeqRecord(proseq, id=proid, description="")                 
                SeqIO.write(orf, writeFile, "fasta")

    ### antisense strand
    dna = dna.reverse_complement()
    dna_len = len(dna)
    for frame in range(3):
        newrna = dna[frame:]
        newrna = newrna[:int(len(newrna)/3)*3]
        end = frame

        
        cds = newrna.translate()
        
        block_list = cds.split('*')
        for block in block_list[:-1]:
            start = end
            end += len(block)*3+3
            Pm = block.find('M')
            start += Pm*3
            block = block[Pm:]
            if len(block) >= aa_cutoff: 
                orf_numer += 1
                proid = r_id + ":" + str(dna_len-end) + "-" + str(dna_len-start) + ":-" + ":dnaORF_%s"%orf_numer 
                proseq = block
                orf = SeqRecord(proseq, id=proid, description="")                 
                SeqIO.write(orf, writeFile, "fasta")


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("-i", "--input", type=str, dest="input", metavar="input.fa", required=True)
    base_group.add_argument("-f", "--outFa", type=str, dest="outFa", metavar="output.fa", required=True)
    base_group.add_argument("-b", "--outBed", type=str, dest="outBed", metavar="output.bed", required=True)
    base_group.add_argument("-c", "--cutoff", type=int, dest="cutoff", metavar="aa_cutoff=20", required=False,
                            default = 20)

    return(parser.parse_args(args))


def main():
    args = sys.argv[1:]

    args = parse_args(args)
    f_in = args.input
    out_fa = args.outFa
    out_bed = args.outBed
    aa_cutoff = args.cutoff

    ### make fasta
    handle_fa = open(out_fa, "w")
    orf_numer = 0
    with open(f_in, "r") as f:      
        for r in SeqIO.parse(f, "fasta"):
            dna = Seq(str(r.seq))      
            frame_translate_6(r.id, dna, aa_cutoff, handle_fa)
            
    handle_fa.close()

    ### convert fasta to bed12
    handle_bed = open(out_bed, "w")
    with open(out_fa, "r") as f_i:
        for r in SeqIO.parse(f_i, "fasta"):
            orf_s, orf_e = map(int, r.id.split(":")[1].split("-"))
            strand = r.id.split(":")[2]

            bed_info = ["MN996528", str(orf_s), str(orf_e), r.id, "0", strand, 
                str(orf_s), str(orf_e), "0,0,0", "1", str(orf_e-orf_s), "0"]
            handle_bed.write("\t".join(bed_info) + "\n")
    handle_bed.close()


if __name__ == "__main__":
    main()
