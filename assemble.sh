#!/bin/bash

a=$1
b=$2
use_center=$3
app=$4
c=$(($a+$b))

python  NP_assembler.py --shell ligand_data/shell-Au.xyz --names ligand_data/BYP.names ligand_data/TMA.names --monomers ligand_data/ByPyCu.xyz ligand_data/TMA.xyz --monomers_resnames BYP TMA --monomers_numbers $a $b --anchor_idxes 1 23 --out structure$c$app.pdb  $use_center
