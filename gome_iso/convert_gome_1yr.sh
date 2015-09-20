#!/bin/bash

year=$1

base_path=/data/GOME1/V4
path=${base_path}/${year}
convert_command1="gdp01_ex -b nynnnnnnnn"
#convert_command1="${base_path}/gdp01_ex_lx -b nynyyynnnn"
#convert_command2=./convert_gome_pmd
convert_command2="convert_gome -b ynnn"

outpath=/tmp/int_data/gome/coords/${year}
outbase=coords
outfile1="/tmp/int_data/temp"
outfile1a="/dev/stdout"

if [ -n ${outpath} ]
then
  mkdir ${outpath}
fi

dirlist1=$(ls ${path} )

echo ${dirlist1}

#ln -s /dev/stdout ${outfile1}.el1

for dir1 in $dirlist1
do
  dirlist2=$(ls ${path}/${dir1})
  for dir2 in $dirlist2
  do
    flist=$(ls ${path}/${dir1}/${dir2})
    for fname in $flist
    do
      outfile2=${outpath}/${outbase}${fname}.dat
      echo cp ${path}/${dir1}/${dir2}/${fname} ${outfile1}
      echo "${convert_command1} ${outfile1} ${outfile1a} | ${convert_command2} ${outfile2}"
      #echo ${convert_command1} ${outfile1} ${outfile1}
      #echo ${convert_command2} ${outfile1}.el1 ${outfile2}
      cp ${path}/${dir1}/${dir2}/${fname} ${outfile1}
      #${convert_command1} ${path}/${dir1}/${dir2}/${fname} ${outfile1}
      ${convert_command1} ${outfile1} ${outfile1}
      ${convert_command2} ${outfile2} < ${outfile1}.el1
      #${convert_command1} ${outfile1} ${outfile1} | ${convert_command2} ${outfile2}
      rm -f ${outfile1}.el1
    done
  done
done

touch ${outpath}

