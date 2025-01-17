#! /bin/zsh
TREE_file=$1;
num=$2;
subname1="$( echo $TREE_file |  sed -e s/.trees// )";
TREE_file_true="$( echo $TREE_file |  grep -Po "[^\/]+$" )";
subname="$( echo $TREE_file |  sed -e s/.trees// | grep -Po "[^\/]+$" )";
LOGFILE="${subname1}.log"
ID_files="${subname}_IDS.txt";
subTREES="${subname}_sub_selected.trees";

echo "Start of $TREE_file formatage ..."
date


line_start="$(grep "STATE_0" -Pno $TREE_file | cut -d ":" -f1)"
((line_start=line_start-1))
line_stop="$( wc -l $TREE_file  | cut -f1 -d " " )"

if [ $num != 0 ]; then

	echo "Start of no burn treefile generation ..."

	re='^[0-9]+$'

	if ! [[ $num =~ $re ]] ; then

		if [ -f $LOGFILE ];then 
			BURNIN="$(loganalyser $LOGFILE | grep -Po "burnIn   <= (\d+)" | grep -Po "\d+")";
		else BURNIN_on="$( tail -n1 $TREE_file  | cut -d " " -f2 | grep -Po "\d+" )";
			echo $BURNIN_on
			BURNIN=$((BURNIN_on / 10))
		fi;

		num=${BURNIN}

	fi;

	num=${num}

	nums="$(grep "STATE_\d+" -Po $TREE_file | sed -e s/"STATE_"// && echo $num )";

	echo "Burnin: $num"

	VALVAL=0

	for u in `echo $nums| tr "\ " '\n' |sort -V | grep -PC1 "$num$"`;do
	
		if [ $u != $num ];then
			VALVAL=$u;
		fi
	
	done

	echo "Close to Burnin: $VALVAL"

	line="$(grep "STATE_$VALVAL " -Pno $TREE_file | cut -d ":" -f1)";

	echo "Line of cut : $line"


else
	
	echo "no BURN-IN ..."
	
	line=$line_start

fi

sed -n "1,$line_start p" $TREE_file > $subTREES

line=${line}
line_stop=${line_stop}

shuf -i $line"-"$line_stop -n 500 > file_ran_values.txt

while read laline ; do 
	
	#echo "$laline"
	sed -n "$laline"p $TREE_file >> $subTREES
		
done < file_ran_values.txt

rm file_ran_values.txt
echo "Finished"
