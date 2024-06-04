file=ShapeFunctionsTri.dat
TYPE=dNdx; /* options are: N, dNdx, dNdy */
case $TYPE

er
gscan -s $file; NBLOCKS=${1}
if case N then
  bscan -s 1 $file; fmn=$fmin; fmx=$fmax
  bscan -s 4 $file; fmn=^min($fmn,$fmin)^; fmx=^max($fmx,$fmax)^
  bscan -s 7 $file; fmn=^min($fmn,$fmin)^; fmx=^max($fmx,$fmax)^
  if [ ${NBLOCKS} > 9 ] then
 bscan -s 10 $file; fmn=^min($fmn,$fmin)^; fmx=^max($fmx,$fmax)^
  endif
else if case dNdx then
  bscan -s 2 $file; fmn=$fmin; fmx=$fmax
  bscan -s 5 $file; fmn=^min($fmn,$fmin)^; fmx=^max($fmx,$fmax)^
  bscan -s 8 $file; fmn=^min($fmn,$fmin)^; fmx=^max($fmx,$fmax)^
  if [ ${NBLOCKS} > 9 ] then
 bscan -s 11 $file; fmn=^min($fmn,$fmin)^; fmx=^max($fmx,$fmax)^
  endif
else if case dNdy then
  bscan -s 3 $file; fmn=$fmin; fmx=$fmax
  bscan -s 6 $file; fmn=^min($fmn,$fmin)^; fmx=^max($fmx,$fmax)^
  bscan -s 9 $file; fmn=^min($fmn,$fmin)^; fmx=^max($fmx,$fmax)^
  if [ ${NBLOCKS} > 9 ] then
 bscan -s 12 $file; fmn=^min($fmn,$fmin)^; fmx=^max($fmx,$fmax)^
  endif
endif

pf $file
nd 1 1
s=^($xmax-$xmin)*0^; xmn=^$xmin-$s^; xmx=^$xmax+$s^
s=^($ymax-$ymin)*0^; ymn=^$ymin-$s^; ymx=^$ymax+$s^
xa $xmn 0.5 $xmx
ya $ymn 0.5 $ymx
set ff
xl 10; yl $xl

* Panel 1
sor 0 ^$yl*1.2^
cp re
cb noce
if case N then
  cont tr 32 10
  cb $fmn $fmx 10
  tt @N sub 1@
  pd 1
else if case dNdx then
  cont tr 10 32
  cb $fmn $fmx 10
  pd 2
  tt @d N sub 1 / d x@
else if case dNdy then
  cont tr 10 32
  cb $fmn $fmx 10
  pd 3
  tt @d N sub 1 / d y@
endif
xt s
yt t
si .4 .4
fr it ac 1111 Tc 0100 lc 0100

* panel 2
sor ^$xl*1.1^ 0
cp re
cb noce
if case N then
  cont tr 32 10
  cb $fmn $fmx 10
  pd 4
  tt @N sub 2@
else if case dNdx then
  cont tr 10 32
  cb $fmn $fmx 10
  pd 5
  tt @d N sub 2 / d x@
else if case dNdy then
  cont tr 10 32
  cb $fmn $fmx 10
  pd 6
  tt @d N sub 2 / d y@
endif
fr it Tc 0000 lc 0000

* panel 3
sor ^-$xl*1.1^ ^-$yl*1.2^
cp re
cb noce
if case N then
  cont tr 32 10
  cb $fmn $fmx 10
  pd 7
  tt @N sub 3@
else if case dNdx then
  cont tr 10 32
  cb $fmn $fmx 10
  pd 8
  tt @d N sub 3 / d x@
else if case dNdy then
  cont tr 10 32
  cb $fmn $fmx 10
  pd 9
  tt @d N sub 3 / d y@
endif
fr it Tc 1100 lc 1100

* panel 4
sor ^$xl*1.1^ 0
if [ $NBLOCKS > 9 ] then
  cp re
  cb noce
  if case N then
 cont tr 32 10
 cb $fmn $fmx 10
 pd 10
 tt @N sub 4@
  else if case dNdx then
 cont tr 10 32
 cb $fmn $fmx 10
 pd 11
 tt @d N sub 4 / d x@
  else if case dNdy then
 cont tr 10 32
 cb $fmn $fmx 10
 pd 12
 tt @d N sub 4 / d y@
  endif
  fr it Tc 1000 lc 1000
endif

sy 1 .4 1
cf -u Hleg ^-$xl*0.05^ ^$yl*1.15^ $fmn $fmx ^$xl/3^ ^$yl/15^
