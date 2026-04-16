#DNA工具测试
from DNAToolkit import *
from utilities import colored
import random

#输入一个随机序列测试
rndDNAStr=''.join([random.choice(Nucleotides)
                   for nuc in range(50)])

DNAStr= validateSeq(rndDNAStr)

# print(f'\nSequence:{colored(DNAStr)}\n')
# print(f'[1]+Sequence Length:{len(DNAStr)}\n')
# print(colored(f'[2]+Nucleotide Frequency:{countFrequency(DNAStr)}\n'))
# print(f'[3]+DNA/RNA Transcription:{colored(transcription(DNAStr))}\n')
# print(f'[4]+Reverse Complememt:{colored(reverse_complement(DNAStr))}\n')
#
# print(f"[5]+DNA String + Reverse Complement:\n5'{colored(DNAStr)}3'")
# print(f"  {''.join(['|' for c in range(len(DNAStr))])}")
# print(f"3'{colored(reverse_complement(DNAStr))}5'\n")
# print(f'[7]+GC Content in Subsection k=5:{gc_content_subsec(DNAStr,k=5)}\n')
#
# print(f'[8]+Aminoacids Sequence from DNA:{translate_seq(DNAStr,0)}\n')
# print(f'[9]+Codon frequency (L):{codon_usage(DNAStr,"L")}\n')

# print(f'[10]+Reading_frames:')
# for frame in gen_reading_frames(DNAStr):
#     print(frame)

print('\n[11]+ All prots in 6 open reading frames:')
for prot in  all_proteins_from_orfs(DNAStr,0,0,True):
    print(f'{prot}')
