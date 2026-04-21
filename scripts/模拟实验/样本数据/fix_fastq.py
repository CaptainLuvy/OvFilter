# -*- coding: utf-8 -*-
import sys

in_fq = sys.argv[1]
out_fq = sys.argv[2]

with open(in_fq, 'r') as f, open(out_fq, 'w') as out:
    while True:
        header = f.readline()
        if not header: break
        if not header.startswith('@'): continue
        
        header = header.strip()
        seq = ""
        
        # Read sequence until '+' is encountered
        while True:
            # We need to remember where we are in the file to handle cases where 
            # we read the next header by mistake
            pos = f.tell()
            line = f.readline()
            if not line: 
                break
            if line.startswith('+'): 
                break
            if line.startswith('@'):
                # We accidentally read the next header! This means this record is broken.
                # Seek back so the next iteration can process this header properly
                f.seek(pos)
                break
            seq += line.strip()
            
        qual = ""
        # Read qualities until length matches sequence length
        while len(qual) < len(seq):
            pos = f.tell()
            line = f.readline()
            if not line: 
                break
            if line.startswith('@'):
                # We accidentally read the next header! 
                f.seek(pos)
                break
            qual += line.strip()
            
        if len(qual) > len(seq):
            qual = qual[:len(seq)]
        elif len(qual) < len(seq):
            qual += '?' * (len(seq) - len(qual))
            
        # Discard sequences with length 0
        if len(seq) == 0:
            continue
            
        out.write("{}\n{}\n+\n{}\n".format(header, seq, qual))
print("Fix completed!")