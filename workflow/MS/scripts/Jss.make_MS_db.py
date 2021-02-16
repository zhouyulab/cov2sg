#!/usr/bin/env python3

import os, sys, argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from interval import Interval


def frame_translate_3(rna, aa_cutoff, orf):
    newORFs = list()

    for frame in range(3):
        newrna = rna[frame:]
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
                proid = orf + ":" + str(start) + "-" + str(end) 
                proseq = block
                orf_new = (proid, proseq, start, end)    
                newORFs.append(orf_new)       

    return newORFs    


def filter_orfs(newORFs, junctions):
    js_s, js_e = map(int, junctions.split("-"))

    ORFs_with_juction = list()
    for newo in newORFs:
        neworf_s, neworf_e = newo[2], newo[3]    
        Interval_orf = Interval(neworf_s, neworf_e)

        if js_s in Interval_orf:
            ORFs_with_juction.append(newo)
    
    return ORFs_with_juction
        


def parse_args(args):
    parser = argparse.ArgumentParser()
    base_group = parser.add_argument_group("Base")

    base_group.add_argument("-p", "--path", type=str, dest="path", metavar="fasta.path", required=True)
    base_group.add_argument("-c", "--cutoff", type=int, dest="cutoff", metavar="aa_cutoff=20", required=False, default = 20)
    base_group.add_argument("-o", "--output", type=str, dest="output", metavar="output.fa", required=True)

    return(parser.parse_args(args))


def main():
    args = sys.argv[1:]
    args = parse_args(args)

    path= args.path
    aa_cutoff = args.cutoff
    f_output = args.output

    ### 
    handle = open(f_output, "w")

    for orf in ["ORF1ab", "S2N"]:
        f_fasta = os.path.join(path, "WIV04.nanopore.%s.iso.fa"%orf)

        ### 3-frame traslateion
        with open(f_fasta, "r") as f:      
            for r in SeqIO.parse(f, "fasta"):
                cdna = Seq(str(r.seq))
                newORFs = frame_translate_3(cdna, aa_cutoff, orf)

                ### only reserve orfs with junction 
                ORFs_with_juction = filter_orfs(newORFs, r.id)
                
                ### output fasta
                
                for newO in ORFs_with_juction:
                    newO_id = newO[0] + ":js_" + r.id.replace("-", "_")
                    newO_seq = newO[1]

                    record = SeqRecord(newO_seq, id=newO_id, description="")
                    SeqIO.write(record, handle, "fasta")
    
    handle.close()

if __name__ == "__main__":
    main()
