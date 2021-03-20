#!/bin/bash
cat resall_1.dat | sed "s/^60/F/g" | sed "s/^50/E/g" | sed "s/^40/D/g" | sed "s/^30/C/g" | sed "s/^20/B/g" | sed "s/^10/A/g" | awk '{{printf "%s ", $1}; for (i=2;i<NF;i++) {printf "%s ", $i-$NF}; printf "\n"}' > resall_1r.dat
cat res_1.dat | sed "s/^60/F/g" | sed "s/^50/E/g" | sed "s/^40/D/g" | sed "s/^30/C/g" | sed "s/^20/B/g" | sed "s/^10/A/g" > res_1r.dat
