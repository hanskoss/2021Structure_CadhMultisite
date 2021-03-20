#!/bin/bash
cp pltodo.dat peaklists.info
for i in $(seq 1 1 $(wc peaklists.info | awk '{print $1}'))
do 
name="$(sed "$i ! d" peaklists.info | awk '{print $1}')""$(sed "$i ! d" peaklists.info | awk '{print $2}')"
./improvepeakl3alt.sh $name
sed "1,$ ! d" $name.csv | awk '{print $5}'
done | sort | uniq | sed "s/^\[/ \[/g" | sed "s/\[/\\\[/g" | sed "s/\]//g" > reslist_x.csv

printf "" > reslist.csv
for i in $(cat reslist_x.csv)
do 
	rewritten=$(echo $i | sed "s/.\[//g" | sed "s/]//g" | sed "s/^\([A-Z]*[0-9]*\)[a-zA-Z]*/\1/g")
	echo $rewritten >> reslist.csv
	for j in $(seq 1 1 $(wc peaklists.info | awk '{print $1}'))
	do
	name="$(sed "$j ! d" peaklists.info | awk '{print $1}')""$(sed "$j ! d" peaklists.info | awk '{print $2}')"
	#$(sed "$j ! d" peaklists.info | awk '{print $1}')
	state=$(echo $rewritten |  sed "s/^\([A-Z]*\).*/\1/g")
        if [[ $state == '' ]]
        then
        state='X'
        fi
#	echo 't
	grep "$i" "$name".csv | awk -v name=$(sed "$j ! d" peaklists.info | awk '{print $1}') -v st=$state -v rn=$rewritten '{printf "%s %s %.4f %.4f %i %s\n", rn, name, $2, $3, $6, st}'
	done
	done > shiftlistb.csv
for i in $(seq 1 1 $(wc peaklists.info | awk '{print $1}'))
do
name="$(sed "$i ! d" peaklists.info | awk '{print $1}')""$(sed "$i ! d" peaklists.info | awk '{print $2}')"
#rm $name.dat
#rm $name.csv
done

#NEED TO COPY reslist.x and shiftlist to _mod files
