#!/bin/bash
plnam=$1
for i in $(seq 1 1 $(($(wc fq3list | awk '{print $1}')))); do echo "TR ""$i"; done > pltodo.dat; for i in $(awk '{print $2}' pltodo.dat); do cat $plnam | awk -v select=$i '{if ($10==select){printf "%i %.4f %4f %s %s %i %i %.4f %4f %s %s %s\n", $1, $5, $6, $8, $9, $11/1000, $12/1000, $14, $13, $15, $16, $17}}' > 'TR'"$i".dat; done; mv 'TR'"$i".dat ref0.dat
tac pltodo.dat | sed "2,$ ! d" | tac > temp; mv temp pltodo.dat; echo "ref 0" >> pltodo.dat
