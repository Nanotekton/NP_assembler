#!/bin/bash

a=$1
b=$2
use_center=$3
app=$4
c=$(($a+$b))

output=structure$c$app.pdb

python  NP_assembler.py --shell ligand_data/shell-Au.xyz --names ligand_data/BYP.names ligand_data/TMA.names --monomers ligand_data/ByPyCu.xyz ligand_data/TMA.xyz --monomers_resnames BYP TMA --monomers_numbers $a $b --anchor_idxes 1 23 --out $output $use_center

cp $output tmp_struct
cp ligand_data/COR.pdb $output
grep -v COR tmp_struct >> $output
rm tmp_struct
python renum.py $output
