#!/usr/bin/env python3
import sys

samples_file = sys.argv[1]
geno_file = sys.argv[2]
out_file = sys.argv[3]

with open(samples_file) as f:
    samples = [x.strip() for x in f if x.strip()]

seqs = {s: [] for s in samples}

iupac = {
    frozenset(['A', 'G']): 'R',
    frozenset(['C', 'T']): 'Y',
    frozenset(['G', 'C']): 'S',
    frozenset(['A', 'T']): 'W',
    frozenset(['G', 'T']): 'K',
    frozenset(['A', 'C']): 'M'
}

def base(ref, alt, gt):
    gt = gt.replace('|', '/').strip()
    if gt in ('.', './.'):
        return 'N'
    if gt in ('0/0', '0'):
        return ref if ref in 'ACGT' else 'N'
    if gt in ('1/1', '1'):
        return alt if alt in 'ACGT' else 'N'
    if gt in ('0/1', '1/0') and ref in 'ACGT' and alt in 'ACGT':
        return iupac.get(frozenset([ref, alt]), 'N')
    return 'N'

with open(geno_file) as f:
    for line in f:
        if not line.strip():
            continue
        ref, alt, *gts = line.rstrip('\n').split('\t')
        if len(gts) != len(samples):
            continue
        for s, gt in zip(samples, gts):
            seqs[s].append(base(ref, alt, gt))

L = len(next(iter(seqs.values())))
assert all(len(x) == L for x in seqs.values())

with open(out_file, 'w') as out:
    out.write('#NEXUS\nBegin data;\n')
    out.write(f'  Dimensions ntax={len(samples)} nchar={L};\n')
    out.write('  Format datatype=dna missing=N gap=- interleave=no;\n')
    out.write('  Matrix\n')
    for s in samples:
        out.write(f'  {s} {"".join(seqs[s])}\n')
    out.write('  ;\nEnd;\n')
