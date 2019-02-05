#!/bin/bash

##-------------------------------------------------------------##
##-- Script to renumber atom types in a lammps input file ---- ##
##-- 
##-- Usage :
##--       ./bead_list_to_lmpsinput.sh N beadlist infile outfile
##--        where N is the number of files to convert 
##--              beadlist is the list of atom types
##--              infile is the first of the lammps input files
##--              outfile is the first of the lammps output file
##--
##--        File names are in the format
##--              path/to/file.999
##-------------------------------------------------------------##

# check arguments
if [ $# -ne 4 ]; then
    echo ""
    echo " Usage :"
    echo "       ./bead_list_to_lmpsinput.sh N beadlist infile outfile"
    echo "        where N is the number of files to convert "
    echo "              beadlist is the list of atom types"
    echo "              infile is the first of the lammps input files"
    echo "              outfile is the first of the lammps output file"
    echo "        File names are in the format"
    echo "              path/to/file.999"
    exit 1
fi

Nfiles=$1
typesfile=$2
infile=$3
outfile=$4

# check number of files
if ! [[ "$Nfiles" =~ ^[0-9]+$ ]]; then
    echo "ERROR :: Must specify an integer number of files."
    exit 1
fi

# check files
if ! [ -f $typesfile ]; then
    echo "ERROR :: File $typesfile not found."
    exit 1
fi

for (( i=1 ; i<=$Nfiles; i++ )); do
    if [ -e ${outfile%.*}.$i ]; then
	echo "ERROR :: File ${outfile%.*}.$i already exists."
	exit 1
    fi
    if ! [  -f ${infile%.*}.$i ]; then
	echo "ERROR :: File ${infile%.*}.$i not found."
	exit 1
    fi
done


# load the types files
#declare -A type
Ndna=0
while read -r line
do
    a=($line)
    if ! [ ${a[0]} == "#" ]; then
	type[${a[0]}]=${a[1]}
	Ndna=$(( $Ndna + 1 ))
    fi
done < $typesfile

# set protein types
Ndnatypes=$( awk 'BEGIN{max=0;}{if ($1!="#"&&$2>max) {max=$2}}END{print max}' $typesfile )
Natoms=$( awk '{if ($2=="atoms") {print $1}}' ${infile%.*}.1 )
i=$(( $Ndna + 1 ))
while [ $i -le $Natoms ]; do
    type[$i]=$(( $Ndnatypes + 1 ))
    i=$(( $i + 1 ))
    type[$i]=$(( $Ndnatypes + 2 ))
    i=$(( $i + 1 ))
    type[$i]=$(( $Ndnatypes + 2 ))
    i=$(( $i + 1 ))
done
if [ -f tmp_types ]; then
    rm tmp_types
fi
for (( i=1; i<=$Natoms ; i++ )); do
    echo "$i ${type[$i]}" >> tmp_types
done

echo "Reading types file $typesfile"
echo "Editing lammps inputs with $Ndna polymer atoms"
echo "                           $Natoms total atoms"
echo "                           $Ndnatypes polymer atom types"
echo "                           $(( $Ndnatypes + 2 )) total atom types"

# edit files
for (( file=1; file<=$Nfiles; file++ )); do

fn=${infile%.*}.$file
ofn=${outfile%.*}.$file

echo -ne "Converting file $fn to $ofn \r"

# split originalfile into three parts
awk 'BEGIN{A=0;}{
if (A==0) {print;}
if ($1=="Atoms") {A=1;}
}' $fn > tmp_top
awk 'BEGIN{A=0;}{
if (A==1 && $1=="Velocities") {A=0;}
if (A==1&&$0!="") {print;}
if ($1=="Atoms")  {A=1;}
}' $fn > tmp_middle
awk 'BEGIN{B=0;}{
if ($1=="Velocities") {B=1;}
if (B==1) {print;}
}' $fn > tmp_bottom

# sort atoms
sort -g tmp_middle > tmp_mid_srt

# change bead types
paste tmp_types tmp_mid_srt | awk '{$5=$2;$1="";$2="";print}' > tmp_middle

# edit top of file
awk -v at=$(( $Ndnatypes + 2 )) '{if ($2=="atom"&&$3=="types") {$1=at}; print}' tmp_top > tmp_top2
awk -v dt=$(( $Ndnatypes )) '{ if ($0=="Masses") {print "Masses\n"; for (i=1;i<=dt;i++) {print i,1}; print i,"0.5\n" i+1,"0.25\n"; print "Atoms\n"; exit 1;}; print}' tmp_top2 > tmp_top

# stick them all back together again
cat tmp_top tmp_middle > $ofn
echo "" >> $ofn
cat tmp_bottom >> $ofn

rm tmp_top tmp_top2 tmp_middle tmp_bottom tmp_mid_srt

done
 
rm tmp_types
echo ""
echo "Done."