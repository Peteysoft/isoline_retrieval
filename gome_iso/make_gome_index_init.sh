year=$1

outfile=$2

coord_path=data/coords/
overlap=75

last_year=$(($year-1))
next_year=$(($year-1))

auxbase=gomecoordfiles

cd ${coord_path}
ls ${last_year} > ${auxbase}${last_year}.txt
ls ${next_year} > ${auxbase}${next_year}.txt

tail -n ${overlap} ${auxbase}${last_year}.txt > ${outfile}
ls ${
head -n ${overlap} ${auxbase}${next_year}.txt > ${outfile}


