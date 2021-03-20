#!/bin/bash
printf "" > reslist.csv
./improvepeakl4.sh $1
data="$1"'.csv'
for b in $(pre=0; 
for z in $(cat $data | sed "2,$ ! d" | awk '{printf "%s\n", $8}' | sed "s/[^0-9]*$//g" | sed "s/\[//g" | sort -k 1); 
do 
if [[ $z != $pre ]]
then pre=$z; echo $z; fi; done); 
do cat $data | awk '{printf "%s %s %s\n", $7, $11, $8}' | awk '{printf "%i %s %s\n", $1, $2, $3}' | sed "s/[^0-9]*$//g"  | sed "s/\[//g" | grep ' '"$b"$ > tempc; 
ab=1; 
printf "" > "$b".tmp; 
echo $b >> reslist.csv; 
for aa in $(seq 1 1 $((1*$(wc fq3list | awk '{print $1}')-1))); 
do ab=$((ab+1));
grep ^"$aa"' ' tempc | awk -v freqx=$(cat fq3list | sed "$ab"' ! d') -v ref=$(sed "$ ! d" tempc | awk '{printf "%i", $2/1000}') '{printf "%s,%i\n", freqx, $2/1000-ref}' >> "$b".tmp;
done; 

cp "$b".tmp tmpa
printf "" > tmpb
for i in $(tac cestrefs)
do
sed "$i ! d" "$b".tmp | sed "s/,/ /g" | awk '{printf "95,%s\n",$2}' >> tmpb
sed "$i d" tmpa > tmpc
cp tmpc tmpa
done
cat tmpb > "$b".csv
if [[ $3 > 0 && $2 > 1 ]]
then
tac tmpa | sed "1,""$3"" d" | tac | sed "1,""$2"" d" >> "$b".csv
elif [[ $3 > 0 ]]
then
tac tmpa | sed "1,""$3"" d" | tac >> "$b".csv
elif [[ $2 > 0 ]]
then
cat tmpa | sed "1,""$2"" d" >> "$b".csv
else
cat tmpa  >> "$b".csv
fi
done

