#!/bin/bash
cat $1.dat | awk -v pat="None" '{if ($8 ~ pat) {if ($9 ~ pat) {print $0} else {printf "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n", $1, $2, $3, $4, $5, $6, $7, $9, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20}} else {print $0}}' | awk -v pat="None" '{if ($9 ~ pat) {if ($8 ~ pat) {print $0} else {printf "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n", $1, $2, $3, $4, $5, $6, $7, $8, $8, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20}} else {print $0}}' | awk -v pat="\[" '{if ($8 ~ pat) {if ($9 ~ pat) {print $0} else {printf "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n", $1, $2, $3, $4, $5, $6, $7, $9, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20}} else {print $0}}' | awk -v pat="\[" '{if ($9 ~ pat) {if ($8 ~ pat) {print $0} else {printf "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n", $1, $2, $3, $4, $5, $6, $7, $8, $8, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20}} else {print $0}}' > $1.csv
