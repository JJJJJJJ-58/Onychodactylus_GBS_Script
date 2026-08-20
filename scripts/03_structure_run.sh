#!/usr/bin/env bash
set -euo pipefail

IN=${1:-strc1.for_structure.txt}
KMIN=1
KMAX=10
REPS=10
THREADS=32
BURNIN=100000
MCMC=2000000

STAMP=$(date +%Y%m%d_%H%M%S)
RESDIR="results_${STAMP}"
mkdir -p "$RESDIR"

N=$(wc -l < "$IN")
L=$(awk '{print (NF-2)/2; exit}' "$IN")

cat > mainparams.auto <<EOF2
#define MAXPOPS       20
#define BURNIN        $BURNIN
#define NUMREPS       $MCMC
#define INFILE        $(readlink -f "$IN")
#define OUTFILE       $(readlink -f "$RESDIR")/test
#define NUMINDS       $N
#define NUMLOCI       $L
#define PLOIDY        2
#define MISSING       0
#define ONEROWPERIND  1
#define LABEL         1
#define POPDATA       1
#define POPFLAG       0
#define LOCDATA       0
#define PHENOTYPE     0
#define EXTRACOLS     0
#define MARKERNAMES   0
#define MAPDISTANCES  0
EOF2

cat > extraparams <<'EOF2'
#define NOADMIX        0
#define FREQSCORR      1
#define LINKAGE        0
#define COMPUTEPROB    1
#define PRINTQHAT      1
EOF2

STRUCTURE=${STRUCTURE:-structure}
MP=$(readlink -f mainparams.auto)
EP=$(readlink -f extraparams)

parallel -j "$THREADS" --halt soon,fail=1 \
  --joblog "$RESDIR/parallel.log" \
  'out='"$RESDIR"'/K{1}_rep{2}; \
   '"$STRUCTURE"' -K {1} -i '"$IN"' -o "$out" \
     -m '"$MP"' -e '"$EP"' -D $((800000 + {1}*100 + {2}))' \
  ::: $(seq "$KMIN" "$KMAX") ::: $(seq 1 "$REPS")
