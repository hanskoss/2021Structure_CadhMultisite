#!/bin/bash
echo 'need fq1list, relaxanalyse.sh'
echo 'need peak list dat file which is spelled out in procname'
echo 'need setlist.csv, space sep, with datasetname and comment#2 experiment type note'
SLhard=$1 #2.5
SLsoft=$2 #0.871
Nfield=$3 #81.09876
datasetname2=$4
datasetname=$datasetname2
lendataset=$(($(wc fq1list | awk '{print $1}')+1))
procname=$5
if [ -f $procname.new ]
then
procnam3="$procname".new
else
procnam3="$procname".dat
fi
echo $procnam3
for x in $(seq 1 1 $lendataset); do echo $datasetname $x 1; done > $procname.schedulem
cp $procname.schedulem $procname.schedulem2
cat $procnam3 | sed "s/\t/ /g" | awk -v pat="None" '{if ($8 ~ pat) {if ($9 ~ pat) {print $0} else {printf "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n", $1, $2, $3, $4, $5, $6, $7, $9, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20}} else {print $0}}' | awk -v pat="None" '{if ($9 ~ pat) {if ($8 ~ pat) {print $0} else {printf "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n", $1, $2, $3, $4, $5, $6, $7, $8, $8, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20}} else {print $0}}' | awk -v pat="\[" '{if ($8 ~ pat) {if ($9 ~ pat) {print $0} else {printf "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n", $1, $2, $3, $4, $5, $6, $7, $9, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20}} else {print $0}}' | awk -v pat="\[" '{if ($9 ~ pat) {if ($8 ~ pat) {print $0} else {printf "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n", $1, $2, $3, $4, $5, $6, $7, $8, $8, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20}} else {print $0}}' | awk -v dsn=$datasetname2 '{printf "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n", $1, dsn, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20}'> $procname.csv
#cat $procname.dat | awk -v pat="\[" '{if ($8 ~ pat) {if ($9 ~ pat) {print $0} else {printf "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n", $1, $2, $3, $4, $5, $6, $7, $9, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20}} else {print $0}}' > $procname.csv
echo $procname $SLhard $SLsoft $Nfield
./relaxanalyse2.sh "$procname" $SLhard $SLsoft $Nfield
./res2namadd.sh
