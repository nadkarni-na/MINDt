#!/bin/bash

rawimdir=$(readlink -e $1)
NIfTIdir=$(readlink -e $2)
brain=$(readlink -e $3)
atlas=$(readlink -e $4)
mask=$(readlink -e $5)
head=$(readlink -e $6)
headweight=$(readlink -e $7)
basetype=$8
tmpdir=$9
registerfunctional=${10}
subpipeline=${11}
Urad=${12}
brainvol=${13}
scale=${14}

echo "1 rawimdir" $rawimdir
echo "2 NIfTIdir" $NIfTIdir
echo "3 brain" $brain
echo "4 atlas" $atlas
echo "5 mask" $mask
echo "6 head" $head
echo "7 headweight" $headweight
echo "8 basetype" $basetype
echo "9 tmpdir" $tmpdir
echo "10 registerfunctional" $registerfunctional
echo "11 subpipeline" $subpipeline
echo "12 Urad" $Urad
echo "13 brainvol" $brainvol
echo "14 scale" $scale

python PVEnDCMtoNIfTI.py /usr/bin/dcmdump $rawimdir $NIfTIdir yes
IDsequencetypes.bash $NIfTIdir
anattotemplate.bash $NIfTIdir yes 3dUnifize $brainvol $scale $brain $atlas $mask $brain $headweight $basetype $Urad
perfFAIREPI.bash $NIfTIdir
perfFAIREPI_spatnorm.bash $NIfTIdir $subpipeline $brain $tmpdir $registerfunctional $scale $brainvol
rs.bash $NIfTIdir $subpipeline $brain $tmpdir $registerfunctional yes $brainvol
