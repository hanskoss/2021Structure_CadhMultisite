#!/bin/bash

datafile=$1.csv
schedulefile=$1.schedulem

function spincheck {
	xa=$1
#	if [[ $xa -eq 68 ]] ####
#	then ####
        grep ' '"$xa"'$' temp > temp2
        for z in $(st=0; 
        	for x in $(cat $schedulefile | sort -k 3 -n | awk '{print $3}')
                do 
                if [[ $st -eq 0 ]]
                	then y=$x; st=1; echo $x
                        else 
                                if [[ $y -ne $x ]]
                                then y=$x; echo $x; 
                                fi; 
                        fi
                        done)
                do  #for each group of peaks
		numberofhits=0
		found=""
		printf "" > temp3
			echo $(for da in $(grep ' '"$z"'$' $schedulefile | awk '{printf "%sXSPACEX%s\n", $1, $2}')
                        do
			
	                        field1=$(echo $da | sed "s/XSPACEX/ /g" | awk '{print $1}')
        	                field2=$(echo $da | sed "s/XSPACEX/ /g" | awk '{print $2}')
				grep "$xa"'$' temp4 | grep '^'"$field1"' '"$field2"' ' |  awk '{printf "%.3f %.3f\n", $3, $4}'
			done) | awk 'BEGIN {eleme=0; sum=0; sum2=0} { eleme+=1; sum+=$1; sum2+=$2} END { printf "%.3f %.3f ", sum/eleme, sum2/eleme }' >> res_"$z".dat
			for da in $(grep ' '"$z"'$' $schedulefile | awk '{printf "%sXSPACEX%s\n", $1, $2}')
                        do 
                        field1=$(echo $da | sed "s/XSPACEX/ /g" | awk '{print $1}')
                        field2=$(echo $da | sed "s/XSPACEX/ /g" | awk '{print $2}')
                        found=$(grep ^"$field1" temp2 | awk '{print $4}' | grep '^'"$field2"'$' | sed "s/ /_/g")
#			echo $field1 $field2 'marker1'	
			if [[ $found != "" ]]; then 
				numberofhits=$((numberofhits+1))
				grep "$xa"'$' temp4 | grep '^'"$field1"' '"$field2"' ' | awk '{printf "%.2f ", $5}' >> res_"$z".dat
			else
				echo $field1 $field2 >> temp3 
			fi
			found=""
                        done
		if [[ $numberofhits -gt 0 ]]
			then
			if [[ $(grep ' '"$z"'$' $schedulefile | awk '{printf "%sXSPACEX%s\n", $1, $2}' | wc | awk '{print $1*2}') -ne $numberofhits ]]
#			then echo 'All '"$numberofhits"' peaks found.'
#			else
                        then
			echo "missing peaks for group  ""$z"" and spin system ""$xa"':'
			cat temp3
			echo '.'
			fi
		fi
        done;
#	fi
}


#cat $datafile | awk '{printf "%s %s %s %i %i %s\n", $2, $5, $6, $7, $11, $9}' | sed "s/}.*//g" | sed "s/{//g" | sed "s/^[^ ]*://g" | sed "2,$ ! d" | sort -k1,1 -k6,6n -k4,4n >> temp

cat $datafile  | awk '{printf "%s %s %s %i %i %s\n", $2, $5, $6, $7, $11/1000, $8}' | sed "s/}.*//g" | sed "s/{//g" | sed "s/^[^ ]*://g" | sed "2,$ ! d" | sed "s/ H\[/ \[/g" | sed -E "s/ ([^ 0-9])([0-9]*)[^ 0-9]*$/ \2 \1/g" | sed -E "s/([0-9])$/\1 0/g"  | sed "s/ A$/ 1000/g" | sed "s/ B$/ 2000/g" | sed "s/ C$/ 3000/g" | sed "s/ D$/ 4000/g" | sed "s/ E$/ 5000/g" | sed "s/ F$/ 6000/g" | awk '{printf "%s %s %s %s %s %i\n", $1, $2, $3, $4, $5, $6+$7}'  | sort -k1,1 -k6,6n -k4,4n > temp
#cat $datafile  | awk '{printf "%s %s %s %i %i %s\n", $2, $5, $6, $7, $11, $8}' | sed "s/}.*//g" | sed "s/{//g" | sed "s/^[^ ]*://g" | sed "2,$ ! d" | sed -E "s/ ([^ 0-9])([0-9]*)[^ 0-9]*$/ \2 \1/g" | sed "s/ H\[/ \[/g" | sed -E "s/([0-9])$/\1 0/g"  | sed "s/ A$/ 1000/g" | sed "s/ B$/ 2000/g" | sed "s/ C$/ 3000/g" | sed "s/ D$/ 4000/g" | awk '{printf "%s %s %s %s %s %i\n", $1, $2, $3, $4, $5, $6+$7}'  | sort -k1,1 -k6,6n -k4,4n > temp
#cat $datafile | awk '{printf "%s %s %s %i %i %s\n", $2, $5, $6, $7, $11, $8}' | sed "s/}.*//g" | sed "s/{//g" | sed "s/^[^ ]*://g" | sed "2,$ ! d" | sed -E "s/ [^0-9]([0-9]*)[^ 0-9]*/ \1/g" | sort -k1,1 -k6,6n -k4,4n >> temp
cp temp temp.init

for x in $(cat $schedulefile | awk '{printf "%sxxSPACExx%s\n", $1, $2}'); do grep ^"$(echo $x | sed "s/xxSPACExx/ /g" | awk '{print $1}')"' ' temp | grep ' '"$(echo $x | sed "s/xxSPACExx/ /g" | awk '{print $2}')"' '; done  | sort -k1,1 -k6,6n -k4,4n | awk 'BEGIN{x=0}{if (x!=$0) {print $0}; x=$0}' | awk 'BEGIN{x=0; y=0; z=0}{if (x!=$1||y!=$4||z!=$6) {print $0}else{printf "%s DUPLICATE\n", $0};x=$1;y=$4;z=$6}' > tempb; 
if [[ $(grep "DUPLICATE" tempb) != "" ]]; then 
grep "DUPLICATE" tempb; 
else echo 'no problem' 
fi


echo "x" > temp
cat tempb >> temp; sleep 0.1

#for x in $(seq 2 2 $(wc temp | awk '{printf "%i", ($1 - 1)}')); do echo $(cat temp | awk '{printf "%s\n", $5}' | sed "$x ! d") $(cat temp | awk '{printf "%s %s %s %s %s %s\n", $1, $2, $3, $4, $6, $5}' | sed "$((x+1)) ! d"); done | awk '{printf "%s %s %s %s %s %i %i\n", $2, $3, $4, $5, $6, $1, $7}' | awk '{printf "%s %s %s %s %.3f %i %s\n", $1, $4, $2, $3, $7/$6, $7, $5}' > temp4

#sleep 10

for x in $(seq 2 2 $(wc temp | awk '{printf "%i", ($1 - 1)}')); do echo $(cat temp | awk '{printf "%s\n", $5}' | sed "$x ! d") $(cat temp | awk '{printf "%s %s %s %s %s %s\n", $1, $2, $3, $4, $6, $5}' | sed "$((x+1)) ! d"); done | awk '{printf "%s %s %s %s %s %i %i\n", $2, $3, $4, $5, $6, $1, $7}' | awk '{printf "%s %s %s %s %.3f %i %i %s\n", $1, $4, $2, $3, $7/($6*2-$7), $6, $7, $5}' > temp4

printf ""  > temp4b

zlist=""
for z in $(st=0;
	for x in $(cat $schedulefile | sort -k 3 -n | awk '{print $3}')
                do 
                if [[ $st -eq 0 ]]
                        then y=$x; st=1; echo $x 
                        else 
                                if [[ $y -ne $x ]]
                                then y=$x; echo $x; 
                                fi; 
                        fi
                        done)
do
echo $z
zlist="$zlist""$z"' '
printf "" > resx_"$z".dat;
printf "" > resy_"$z".dat
printf "" > resall_"$z".dat
printf "SS CSH unique_name " > res_"$z".dat; for xg in $(sed "s/ /@/g" $schedulefile | grep '@'"$z"'$'); do printf "$xg "; done >> res_"$z".dat; printf "\n" >> res_"$z".dat
done
echo $zlist 'zlistecho'
sta=0
for xa in $(cat temp | sort -k 6 -n | awk '{print $6}') #for each spin system
do 
	if [[ $sta -eq 0 ]]
		then ya=$xa; sta=1;
		for z in $zlist
		do
		echo "$xa" 'xx'
		printf "$xa"' ' >> res_"$z".dat
		if [[ $(cat temp | grep ^"$(grep ' '"$z"'$' $schedulefile | awk '{printf "%s\n", $1}' | sed "1 ! d")"' ' | grep ' '"$xa"$) != "" ]]
		then
		cat temp | grep ^"$(grep ' '"$z"'$' $schedulefile | awk '{printf "%s\n", $1}' | sed "1 ! d")"' ' | grep ' '"$xa"$ | awk -v st=$xa 'BEGIN{printf "%s ", st}{printf "%s ", $5}END{printf "\n"}' >> resall_"$z".dat
#                cat temp4 | grep ^"$(grep ' '"$z"'$' $schedulefile | awk '{printf "%s\n", $1}' | sed "1 ! d")"' ' | grep ' '"$xa"$ | awk -v xa=$(cat temp4 | grep ^"$(grep ' '"$z"'$' $schedulefile | awk '{printf "%s\n", $1}')"' ' | grep ' '"$xa"$ | awk 'BEGIN{x=0;y=0;a=0;b=0;c=0;d=0}{x+=$6;y+=$7;a+=($6*$6);b+=($7*$7);c+=1}END{printf "%i\n", sqrt(a/c)}') -v xb=$(cat temp4 | grep ^"$(grep ' '"$z"'$' $schedulefile | awk '{printf "%s\n", $1}')"' ' | grep ' '"$xa"$ | awk 'BEGIN{x=0;y=0;a=0;b=0;c=0;d=0}{x+=$6;y+=$7;a+=($6*$6);b+=($7*$7);c+=1}END{printf "%i\n", sqrt(b/c)}') 'BEGIN{x=0;y=0;a=0;b=0;c=0;d=0}{a+=($6-xa)^2;b+=($7-xb)^2;c+=1}END{printf "%s %s %i %i %i %i %i\n", $1, $8, xa, sqrt(a/c), xb, sqrt(b/c), c}' >> resx_"$z".dat
		fi
		done
		spincheck "$xa"

		
		for z in $zlist
		do
		printf "\n" >> res_"$z".dat
		done
	else 
		if [[ $ya -ne $xa ]]
		then ya=$xa
			for z in $zlist
			do
			printf "$xa"' ' >> res_"$z".dat
			if [[ $(cat temp | grep ^"$(grep ' '"$z"'$' $schedulefile | awk '{printf "%s\n", $1}' | sed "1 ! d")"' ' | grep ' '"$xa"$) != "" ]]
	                then
			echo $z, 'znow', $xa
			cat temp | grep ^"$(grep ' '"$z"'$' $schedulefile | awk '{printf "%s\n", $1}' | sed "1 ! d")"' ' | grep ' '"$xa"$ | awk -v st=$xa 'BEGIN{printf "%s ", st}{printf "%s ", $5}END{printf "\n"}' >> resall_"$z".dat
##			cat temp4 | grep ^"$(grep ' '"$z"'$' $schedulefile | awk '{printf "%s\n", $1}')"' ' | grep ' '"$xa"$ | awk -v xa=$(cat temp4 | grep ^"$(grep ' '"$z"'$' $schedulefile | awk '{printf "%s\n", $1}')"' ' | grep ' '"$xa"$ | awk 'BEGIN{x=0;y=0;a=0;b=0;c=0;d=0}{x+=$6;y+=$7;a+=($6*$6);b+=($7*$7);c+=1}END{printf "%i\n", sqrt(a/c)}') -v xb=$(cat temp4 | grep ^"$(grep ' '"$z"'$' $schedulefile | awk '{printf "%s\n", $1}')"' ' | grep ' '"$xa"$ | awk 'BEGIN{x=0;y=0;a=0;b=0;c=0;d=0}{x+=$6;y+=$7;a+=($6*$6);b+=($7*$7);c+=1}END{printf "%i\n", sqrt(b/c)}') 'BEGIN{x=0;y=0;a=0;b=0;c=0;d=0}{a+=($6-xa)^2;b+=($7-xb)^2;c+=1}END{printf "%s %s %i %i %i %i %i\n", $1, $8, xa, sqrt(a/c), xb, sqrt(b/c), c}' >> resx_"$z".dat
##file=resx_"$z".dat; one=$(grep ^"$(cat $file | sed "1 ! d" | awk '{printf "%i", $7-1}')"' ' ttest_10.dat | awk '{print $2}'); two=$(grep ^"$(cat $file | sed "1 ! d" | awk '{print $7*2-1}')"' ' ttest_10.dat | awk '{print $3}'); echo $one $two; cat $file | awk '{printf "%s %s %s %s %s %s %s %.2f %.2f %.2f\n", $1, $2, $3, $4, $5, $6, $7, $3/($4/sqrt($7)), $5/($6/sqrt($7)), sqrt(($3-$5)^2)/(sqrt($4^2/$7+sqrt($6^2/$7)))}' | awk -v one=$one -v two=$two '{if ($8 > one) {xa=1}else {xa=0}; if ($9 > one) {xb=1}else {xb=0}; if ($10 > two) {xc=1}else {xc=0}; printf "%s %s %s %s %s %s %s %s %s %s %s %s %s\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, xa, xb, xc}' >> resy_"$z".dat
			fi
			done
			echo $xa
			spincheck "$xa"
			for z in $zlist
			do printf "\n" >> res_"$z".dat
			done
		fi
	fi;
done
for z in $zlist
do
grep " [^ ]" res_"$z".dat | grep -v "0.000 0.000" > temp5
grep -v "^.[^ ]* $" resall_"$z".dat > temp55
mv temp55 resall_"$z".dat
mv temp5 res_"$z".dat
done
printf "" > procsettings.dat; for y in $(seq 1 1 $(for x in $(sed "$ ! d" $1.schedulem | awk '{print $3}'); do echo $x; done)); do echo 'res_'"$y" >> procsettings.dat; done
if [[ -f fq1list ]]
then
for file in $(cat procsettings.dat); do echo "unique_name given_name offset/ppm spinlock/kHz duration/us" > $file.dat.headerdat; for stringcheck in $(sed "1 ! d" $file.dat | awk '{for (i=4; i<=NF; i++) {print $i}}'); do givenname='placeholder'; x1=$(echo $stringcheck | sed "s/@/ /g" | awk '{print $1}'); x2=$(echo $stringcheck | sed "s/@/ /g" | awk '{print $2}'); awk '{printf "%s %s %s %s %s %s\n", $22, $13, $6, $12, $2, $20}' setlist.csv | grep ' '"$x2"$ | sed "1 ! d"; done; done | awk '{print $5}' > tempnum.dat; i=2; for x in $(cat tempnum.dat | sort -u); do echo $x $(sed "$i ! d" fq1list); i=$((i+1)); done > tempnum2.dat; for x in $(cat tempnum.dat); do grep ^"$x"' ' tempnum2.dat | awk '{print $2}'; done > tempnum3.dat; rm tempnum.dat tempnum2.dat
fi ###### this re-sorts the fq1list offsets
echo $file "thatitis"
for file in $(cat procsettings.dat); do echo "unique_name given_name offset/ppm spinlock/kHz duration/us" > $file.dat.headerdat; for stringcheck in $(sed "1 ! d" $file.dat | awk '{for (i=4; i<=NF; i++) {print $i}}'); do givenname='placeholder'; echo $givenname; x1=$(echo $stringcheck | sed "s/@/ /g" | awk '{print $1}'); x2=$(echo $stringcheck | sed "s/@/ /g" | awk '{print $2}'); awk '{printf "%s %s %s %s %s %s\n", $22, $13, $6, $12, $18, $20}' setlist.csv | grep ' '"$x2"$ | sed "1 ! d" | awk -v x1=$x1 -v x2=$x2 -v x=$stringcheck -v givenname=$givenname -v slps=$3 -v slph=$2 -v mhz=$4 '{printf "%s %s %s %s %s %s %s %s\n", x, givenname, $4, slph, $2, $5, mhz, slps}' >> $file.dat.headerdat; done; done

echo 'abc1'
if [[ -f fq1list ]]
then
sed "1 ! d" $file.dat.headerdat > $file.dat.headerdat2
for x in $(seq 2 1 $(wc $file.dat.headerdat | awk '{print $1}'))
do
sed "$x ! d" $file.dat.headerdat | sed "s/placeholder/""$(sed "$((x-1)) ! d" tempnum3.dat)""/g" | awk -v substhh=$(sed "$(($x - 1))""! d" r1rtypelist.dat) '{printf "%s %s %s %s %s %s %s %s\n", $1, $2, $3, $4, $5, substhh, $7, $8}'
done >> $file.dat.headerdat2
fi
