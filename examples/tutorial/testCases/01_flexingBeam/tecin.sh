#!/bin/sh

cat > TECIN.DAT.tmpl << EOF
Benchmark for the Rayleigh Taylor instability
       10777	   12345           1
    2    1    1    1    0    0    0    0    0    0    0
    0
.so tecin.dat.partf.nps
.so tecin.dat.partf.elm
.so tecin.dat.bcs
    0
    0
  sec
  0.5
    0    1    1    2    1    0    0    0    0    0
    0
           0           0           0           0           0           0           0           0
           1        7.5E10           0.3  5.000000e+31           1.0           1.0           1.0
end material props
           0.0           0.0           0.0
EOF

# Count number of nodes/elements
NUMNP=`wc tecin.dat.nps | awk '{print $1-1}'`
NUMEL=`wc tecin.dat.elm | awk '{print $1-1}'`
#NUMFN=`awk '{if ($1 != "end") print $0}' tecin.dat.spl | wc | awk '{print $1}'`
#NUMFN  = "0"
#NUMAT =`awk '{if ($1 != "end") print $0}' tecin.dat.elm | sort -nuk2 | wc | awk '{print $1}'`
#NUMAT  = 2
#NUMWNK=`awk '{if ($1 != "end") print $0}' tecin.dat.wnk | wc | awk '{print $1}'`
#NUMWNK ="0"

echo "Number of nodes:       $NUMNP" 1>&2
echo "Number of elements:    $NUMEL" 1>&2
echo "Number of materials:   $NUMAT" 1>&2
echo "Number of split nodes: $NUMFN" 1>&2
echo "Number of Winklers:    $NUMWNK" 1>&2

# Make TECIN.DAT
awk '{
    if (NR==2) {printf("%12i%12i%12i\n"),'"$NUMNP"','"$NUMEL"','"1"'}
    else if (NR==14) {printf("%12i%12i%12i%12i%12i%12i%12i%12i\n"),$1,$2,$3,$4,$5,$6,$7,$8}
    else {print $0}
}' TECIN.DAT.tmpl > TECIN.DAT


rm TECIN.DAT.tmpl

exit

