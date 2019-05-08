#    TaxidDetermination.sh - part of RIEMS - Reliable Information Extraction from Metagenomic Sequence datasets
#    Copyright (C) 2009-2016  Ariane Belka, Maria Jenckel, Matthias Scheuch, Dirk Hoeper
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


#!/bin/bash

function get_species() {
taxid=$tid                                                                                                  # assign tid to taxid (basically the same)
isSpecies=FALSE                                                                                             # set isSpecies to FALSE to start with

while [[ ${isSpecies} != TRUE ]]                                                                            # as long as $isSpecies is false, do ...
    do
        taxTemp=`grep "^\<${taxid}\>" ${taxdir}/nodes.dmp`                                                  # retrieve the line with information on the current TAX-ID from the nodes.dmp database
        if (( `echo ${taxTemp} | grep -wc "species"` > 0 ))                                                 # if taxonomy level of the current TAX-ID is species, then ...
            then                                                                                          
                isSpecies=TRUE                                                                              # set isSpecies TRUE to terminate the while-loop
        elif [[ $taxid == 1 ]]                                                                              # else if tax-id is 1 (reached root), then ...
            then                            
                taxid=$tid                                                                                  # set taxid back to tid
                isSpecies=TRUE                                                                              # and leave if an while loop
        elif [[ -z $taxTemp ]]                                                                              # is taxTemp is empty, then ...
            then
                break                                                                                       # leave loop
        else                            
            taxid=`echo ${taxTemp} | cut -d "|" -f2 | tr -d [:blank:]`                                      # then retrieve parent-TAX-ID (parent in column 2 -> cut -f2)
        fi                                                                                              
    done                        
                    
parentTax=$taxid    
tid=$parentTax                                                                                              # assign taxid to parentTax
}

function remove_ambiguities() {

csplit -s ${input} /\>/ {*} -f "csplit" -n 6 -z -q                                                          # split the input fasta-formatted sequence dataset into individual files for each individual entry
for split in csplit*                                                                                        # loop to execute the following code for sequence renaming and splitting into overlapping sub-sequences on all available sequence files
do
    grep ">" ${split} > ${split}_temp1                                                                      # extract the header line from the current file to temp1
    grep -v ">" ${split} | tr -d "[KMRYSWBVHDNkmryswbvhdn]" | sed '/^$/d' >> ${split}_temp1                 # extract the sequence from the current file, delete "N" and "n" characters from the sequence and write to file temp2
done
cat *_temp1 > ${input}.changed                                                                              # cat all changed files
rm *_temp1 ; rm csplit*                                                                                     # and remove subset and csplit files

}

function TaxidDetermination()                                                                               # for determination of tax-ids for given gis                                        
{
  case "${line}" in                                                                                         # line habors given gi
    1*) case "${line}" in                                                                                   # if 1st character of line (gi) begins with 1, continue with case in 1*)
          1)   grep -m1 "^\<$line\>" ${taxonomieiddb}-1.dmp >> gi-tid.txt ;;                                # if gi=1, grep gi and corresponding taxid in ${taxonomieiddb}-1.dmp
          10*) grep -m1 "^\<$line\>" ${taxonomieiddb}-10.dmp >> gi-tid.txt ;;                               # check 2nd character in gi and grep in corresponding file 
          11*) grep -m1 "^\<$line\>" ${taxonomieiddb}-11.dmp >> gi-tid.txt ;;                               # | 
          12*) grep -m1 "^\<$line\>" ${taxonomieiddb}-12.dmp >> gi-tid.txt ;;                               # |  
          13*) grep -m1 "^\<$line\>" ${taxonomieiddb}-13.dmp >> gi-tid.txt ;;                               # | 
          14*) grep -m1 "^\<$line\>" ${taxonomieiddb}-14.dmp >> gi-tid.txt ;;                               # |
          15*) grep -m1 "^\<$line\>" ${taxonomieiddb}-15.dmp >> gi-tid.txt ;;                               # |
          16*) grep -m1 "^\<$line\>" ${taxonomieiddb}-16.dmp >> gi-tid.txt ;;                               # |
          17*) grep -m1 "^\<$line\>" ${taxonomieiddb}-17.dmp >> gi-tid.txt ;;                               # |
          18*) grep -m1 "^\<$line\>" ${taxonomieiddb}-18.dmp >> gi-tid.txt ;;                               # |
          19*) grep -m1 "^\<$line\>" ${taxonomieiddb}-19.dmp >> gi-tid.txt ;;                               # |
        esac ;;                                                                                             # |
    2*) case "${line}" in                                                                                   # | 
          2)  grep -m1 "^\<$line\>" ${taxonomieiddb}-2.dmp >> gi-tid.txt ;;                                 # |
          20*) grep -m1 "^\<$line\>" ${taxonomieiddb}-20.dmp >> gi-tid.txt ;;                               # |
          21*) grep -m1 "^\<$line\>" ${taxonomieiddb}-21.dmp >> gi-tid.txt ;;                               # |
          22*) grep -m1 "^\<$line\>" ${taxonomieiddb}-22.dmp >> gi-tid.txt ;;                               # |
          23*) grep -m1 "^\<$line\>" ${taxonomieiddb}-23.dmp >> gi-tid.txt ;;                               # |
          24*) grep -m1 "^\<$line\>" ${taxonomieiddb}-24.dmp >> gi-tid.txt ;;                               # |
          25*) grep -m1 "^\<$line\>" ${taxonomieiddb}-25.dmp >> gi-tid.txt ;;                               # |
          26*) grep -m1 "^\<$line\>" ${taxonomieiddb}-26.dmp >> gi-tid.txt ;;                               # |
          27*) grep -m1 "^\<$line\>" ${taxonomieiddb}-27.dmp >> gi-tid.txt ;;                               # |
          28*) grep -m1 "^\<$line\>" ${taxonomieiddb}-28.dmp >> gi-tid.txt ;;                               # |
          29*) grep -m1 "^\<$line\>" ${taxonomieiddb}-29.dmp >> gi-tid.txt ;;                               # |
        esac ;;                                                                                             # |   
    3*) case "${line}" in                                                                                   # |               
          3)  grep -m1 "^\<$line\>" ${taxonomieiddb}-3.dmp >> gi-tid.txt ;;                                 # |
          30*) grep -m1 "^\<$line\>" ${taxonomieiddb}-30.dmp >> gi-tid.txt ;;                               # |
          31*) grep -m1 "^\<$line\>" ${taxonomieiddb}-31.dmp >> gi-tid.txt ;;                               # |
          32*) grep -m1 "^\<$line\>" ${taxonomieiddb}-32.dmp >> gi-tid.txt ;;                               # |
          33*) grep -m1 "^\<$line\>" ${taxonomieiddb}-33.dmp >> gi-tid.txt ;;                               # |
          34*) grep -m1 "^\<$line\>" ${taxonomieiddb}-34.dmp >> gi-tid.txt ;;                               # |
          35*) grep -m1 "^\<$line\>" ${taxonomieiddb}-35.dmp >> gi-tid.txt ;;                               # |
          36*) grep -m1 "^\<$line\>" ${taxonomieiddb}-36.dmp >> gi-tid.txt ;;                               # |
          37*) grep -m1 "^\<$line\>" ${taxonomieiddb}-37.dmp >> gi-tid.txt ;;                               # |
          38*) grep -m1 "^\<$line\>" ${taxonomieiddb}-38.dmp >> gi-tid.txt ;;                               # |
          39*) grep -m1 "^\<$line\>" ${taxonomieiddb}-39.dmp >> gi-tid.txt ;;                               # |
        esac ;;                                                                                             # |   
    4*) case "${line}" in                                                                                   # |               
          4)  grep -m1 "^\<$line\>" ${taxonomieiddb}-4.dmp >> gi-tid.txt ;;                                 # |
          40*) grep -m1 "^\<$line\>" ${taxonomieiddb}-40.dmp >> gi-tid.txt ;;                               # |
          41*) grep -m1 "^\<$line\>" ${taxonomieiddb}-41.dmp >> gi-tid.txt ;;                               # |
          42*) grep -m1 "^\<$line\>" ${taxonomieiddb}-42.dmp >> gi-tid.txt ;;                               # |
          43*) grep -m1 "^\<$line\>" ${taxonomieiddb}-43.dmp >> gi-tid.txt ;;                               # |
          44*) grep -m1 "^\<$line\>" ${taxonomieiddb}-44.dmp >> gi-tid.txt ;;                               # |
          45*) grep -m1 "^\<$line\>" ${taxonomieiddb}-45.dmp >> gi-tid.txt ;;                               # |
          46*) grep -m1 "^\<$line\>" ${taxonomieiddb}-46.dmp >> gi-tid.txt ;;                               # |
          47*) grep -m1 "^\<$line\>" ${taxonomieiddb}-47.dmp >> gi-tid.txt ;;                               # |
          48*) grep -m1 "^\<$line\>" ${taxonomieiddb}-48.dmp >> gi-tid.txt ;;                               # |
          49*) grep -m1 "^\<$line\>" ${taxonomieiddb}-49.dmp >> gi-tid.txt ;;                               # |
        esac ;;                                                                                             # |   
    5*) case "${line}" in                                                                                   # |               
          5)  grep -m1 "^\<$line\>" ${taxonomieiddb}-5.dmp >> gi-tid.txt ;;                                 # |
          50*) grep -m1 "^\<$line\>" ${taxonomieiddb}-50.dmp >> gi-tid.txt ;;                               # |
          51*) grep -m1 "^\<$line\>" ${taxonomieiddb}-51.dmp >> gi-tid.txt ;;                               # |
          52*) grep -m1 "^\<$line\>" ${taxonomieiddb}-52.dmp >> gi-tid.txt ;;                               # |
          53*) grep -m1 "^\<$line\>" ${taxonomieiddb}-53.dmp >> gi-tid.txt ;;                               # |
          54*) grep -m1 "^\<$line\>" ${taxonomieiddb}-54.dmp >> gi-tid.txt ;;                               # |
          55*) grep -m1 "^\<$line\>" ${taxonomieiddb}-55.dmp >> gi-tid.txt ;;                               # |
          56*) grep -m1 "^\<$line\>" ${taxonomieiddb}-56.dmp >> gi-tid.txt ;;                               # |
          57*) grep -m1 "^\<$line\>" ${taxonomieiddb}-57.dmp >> gi-tid.txt ;;                               # |
          58*) grep -m1 "^\<$line\>" ${taxonomieiddb}-58.dmp >> gi-tid.txt ;;                               # |
          59*) grep -m1 "^\<$line\>" ${taxonomieiddb}-59.dmp >> gi-tid.txt ;;                               # |
        esac ;;                                                                                             # |   
    6*) case "${line}" in                                                                                   # |               
          6)  grep -m1 "^\<$line\>" ${taxonomieiddb}-6.dmp >> gi-tid.txt ;;                                 # |
          60*) grep -m1 "^\<$line\>" ${taxonomieiddb}-60.dmp >> gi-tid.txt ;;                               # |
          61*) grep -m1 "^\<$line\>" ${taxonomieiddb}-61.dmp >> gi-tid.txt ;;                               # |
          62*) grep -m1 "^\<$line\>" ${taxonomieiddb}-62.dmp >> gi-tid.txt ;;                               # |
          63*) grep -m1 "^\<$line\>" ${taxonomieiddb}-63.dmp >> gi-tid.txt ;;                               # |
          64*) grep -m1 "^\<$line\>" ${taxonomieiddb}-64.dmp >> gi-tid.txt ;;                               # |
          65*) grep -m1 "^\<$line\>" ${taxonomieiddb}-65.dmp >> gi-tid.txt ;;                               # |
          66*) grep -m1 "^\<$line\>" ${taxonomieiddb}-66.dmp >> gi-tid.txt ;;                               # |
          67*) grep -m1 "^\<$line\>" ${taxonomieiddb}-67.dmp >> gi-tid.txt ;;                               # |
          68*) grep -m1 "^\<$line\>" ${taxonomieiddb}-68.dmp >> gi-tid.txt ;;                               # |
          69*) grep -m1 "^\<$line\>" ${taxonomieiddb}-69.dmp >> gi-tid.txt ;;                               # |
        esac ;;                                                                                             # |   
    7*) case "${line}" in                                                                                   # |             
          7)  grep -m1 "^\<$line\>" ${taxonomieiddb}-7.dmp >> gi-tid.txt ;;                                 # |
          70*) grep -m1 "^\<$line\>" ${taxonomieiddb}-70.dmp >> gi-tid.txt ;;                               # |
          71*) grep -m1 "^\<$line\>" ${taxonomieiddb}-71.dmp >> gi-tid.txt ;;                               # |
          72*) grep -m1 "^\<$line\>" ${taxonomieiddb}-72.dmp >> gi-tid.txt ;;                               # |
          73*) grep -m1 "^\<$line\>" ${taxonomieiddb}-73.dmp >> gi-tid.txt ;;                               # |
          74*) grep -m1 "^\<$line\>" ${taxonomieiddb}-74.dmp >> gi-tid.txt ;;                               # |
          75*) grep -m1 "^\<$line\>" ${taxonomieiddb}-75.dmp >> gi-tid.txt ;;                               # |
          76*) grep -m1 "^\<$line\>" ${taxonomieiddb}-76.dmp >> gi-tid.txt ;;                               # |
          77*) grep -m1 "^\<$line\>" ${taxonomieiddb}-77.dmp >> gi-tid.txt ;;                               # |
          78*) grep -m1 "^\<$line\>" ${taxonomieiddb}-78.dmp >> gi-tid.txt ;;                               # |
          79*) grep -m1 "^\<$line\>" ${taxonomieiddb}-79.dmp >> gi-tid.txt ;;                               # |
        esac ;;                                                                                             # |   
    8*) case "${line}" in                                                                                   # |           
          8)  grep -m1 "^\<$line\>" ${taxonomieiddb}-8.dmp >> gi-tid.txt ;;                                 # |
          80*) grep -m1 "^\<$line\>" ${taxonomieiddb}-80.dmp >> gi-tid.txt ;;                               # |
          81*) grep -m1 "^\<$line\>" ${taxonomieiddb}-81.dmp >> gi-tid.txt ;;                               # |
          82*) grep -m1 "^\<$line\>" ${taxonomieiddb}-82.dmp >> gi-tid.txt ;;                               # |
          83*) grep -m1 "^\<$line\>" ${taxonomieiddb}-83.dmp >> gi-tid.txt ;;                               # |
          84*) grep -m1 "^\<$line\>" ${taxonomieiddb}-84.dmp >> gi-tid.txt ;;                               # |
          85*) grep -m1 "^\<$line\>" ${taxonomieiddb}-85.dmp >> gi-tid.txt ;;                               # |
          86*) grep -m1 "^\<$line\>" ${taxonomieiddb}-86.dmp >> gi-tid.txt ;;                               # |
          87*) grep -m1 "^\<$line\>" ${taxonomieiddb}-87.dmp >> gi-tid.txt ;;                               # |
          88*) grep -m1 "^\<$line\>" ${taxonomieiddb}-88.dmp >> gi-tid.txt ;;                               # |
          89*) grep -m1 "^\<$line\>" ${taxonomieiddb}-89.dmp >> gi-tid.txt ;;                               # |
        esac ;;                                                                                             # |       
    9*) case "${line}" in                                                                                   # |                   
          9)  grep -m1 "^\<$line\>" ${taxonomieiddb}-9.dmp >> gi-tid.txt ;;                                 # |
          90*) grep -m1 "^\<$line\>" ${taxonomieiddb}-90.dmp >> gi-tid.txt ;;                               # |
          91*) grep -m1 "^\<$line\>" ${taxonomieiddb}-91.dmp >> gi-tid.txt ;;                               # |
          92*) grep -m1 "^\<$line\>" ${taxonomieiddb}-92.dmp >> gi-tid.txt ;;                               # |
          93*) grep -m1 "^\<$line\>" ${taxonomieiddb}-93.dmp >> gi-tid.txt ;;                               # |
          94*) grep -m1 "^\<$line\>" ${taxonomieiddb}-94.dmp >> gi-tid.txt ;;                               # |
          95*) grep -m1 "^\<$line\>" ${taxonomieiddb}-95.dmp >> gi-tid.txt ;;                               # |
          96*) grep -m1 "^\<$line\>" ${taxonomieiddb}-96.dmp >> gi-tid.txt ;;                               # |
          97*) grep -m1 "^\<$line\>" ${taxonomieiddb}-97.dmp >> gi-tid.txt ;;                               # |
          98*) grep -m1 "^\<$line\>" ${taxonomieiddb}-98.dmp >> gi-tid.txt ;;                               # |
          99*) grep -m1 "^\<$line\>" ${taxonomieiddb}-99.dmp >> gi-tid.txt ;;                               # |
        esac ;;                                                                                             # |       
  esac                                                                                                      # V
}                                                                                                   

function GiDownload()                                                                                       # download sequences for given gi list
{
  more $i | $blastdir/blastdbcmd -db $ntdb -entry_batch - -out ${i}.fna 1>/dev/null 2>&1                         # use blastdbcmd do download sequences from nucleotide database and write into file.fna, $i is part-* file
}


function DoubbleKill()                                                                                      # removes duplicated gis after GiDownload
{
csplit -f ${out}/split ${out}/TID-${parentTax}.fna '/^>/' {*} & wait 1>/dev/null 2>&1      # split the sequence file for Tax-id by header so each file contains one sequence

i=0                                                                                                         # set counter to 0
while (( i < 10 ))                                                                                          # as long as i is smaller than 10, do....
    do 
        echo ${out}/split*${i} | xargs md5sum >> ${out}/md5.txt 2>/dev/null                                                         # calculate md5sum for all split-files that end with i
        let i=i+1                                                                                           # add 1 to i
    done 1>/dev/null 2>&1 ; wait 1>/dev/null 2>&1                                                                           # wait for all background processes to finish
cut -d " " -f3 ${out}/md5.txt | sort > ${out}/all.txt                                                       # cut the 3rd column of md5.txt (two spaces between columns) to get all filenames
sort ${out}/md5.txt | uniq -w32 | cut -d " " -f3 | sort  > ${out}/uniq.txt                                  # sort md5.txt and make uniq by comparing the first 32 characters (md5sum) and cut filenames again and sort

comm -13 ${out}/uniq.txt ${out}/all.txt | sed -e 's/^[ \t]*//' > ${out}/to_remove.txt                       # compare uniq files and all files and get all files that are only in all (-13), delete tabs at beginning of line and write them to to_remove
comm -12 ${out}/uniq.txt ${out}/all.txt | sed -e 's/^[ \t]*//' > ${out}/to_keep.txt                         # compare uniq files and all files and get all files that are in both lists (-12), delete tabs at beginning of line and write them to to_keep

cd ${out}                                                                                                   # change working directory to ${out}
xargs rm < ${out}/to_remove.txt 1>/dev/null 2>&1                                                                 # remove all files which are listed in to_remove, xargs accepts a list of arguments and runs command (here rm) for this list  
                                                                                                            # faster than rm alone and can also use large lists of arguments, only rm might give error for to long list
i=0                                                                                                         # set counter back to 0
while (( i < 10 ))                                                                                          # as long as i is less than 10, do ...
    do 
        echo ${out}/split*${i} | xargs cat #2>/dev/null                                                                              # cat all split files that end with $i
        echo                                                                                                # echo to start new line after each cat
        let i=i+1                                                                                           # add 1 to i
    done >> ${out}/TID-${parentTax}_${name}.fna                                                  # write everythink into one file in Refseqdir

sed '/>/!s/[KMRYSWBVHDNkmryswbvhdn]//g' ${out}/TID-${parentTax}_${name}.fna > ${refseqdir}/TID-${parentTax}_${name}.fna     #remove all ambiguities in reference sequences except in lines containing an ">" (header)
    
all=`grep "^>" ${refseqdir}/TID-${parentTax}_${name}.fna | sort | wc -l`                                    # grep all header, sort them and count the lines (equal to count all header)
uniq=`grep "^>" ${refseqdir}/TID-${parentTax}_${name}.fna | sort | uniq | wc -l`                            # grep all header, sort them uniq and count the lines (equal to count all uniq header)
if (( $all > $uniq ))                                                                                       # if number of all header is greater than number of uniq header, then ...
    then                                                                                                    
        echo "DAMN IT"                                                                                      # user info 
fi

xargs rm < ${out}/to_keep.txt 1>/dev/null 2>&1                                                                   # remove all split files that are left using xargs (see above)
rm ${out}/md5.txt;  rm ${out}/all.txt; rm ${out}/uniq.txt  1>/dev/null 2>&1                                                 # remove files that will not be needed anymore
rm ${out}/to_remove.txt; rm ${out}/to_keep.txt 1>/dev/null 2>&1                                                  # remove files that will not be needed anymore

cd ${arbeitsvz}                                                                                             # change to working directory
rm -r ${refseqdir}/reftmp_$USER  1>/dev/null 2>&1                                                                                # and remove reftmp folder
}

# Schnellere Funktion zum Auffinden der Family- und Superkingdom Taxids
function FamSkTaxDetermination()                                                                            # function to get family- and superkingdom tax-ids
{
taxid=$tid                                                                                                  # assign tid to taxid (basically the same)
famtax=NA
sktax=NA
isFamily=FALSE                                                                                              # set isFamily to FALSE to start with until loop
isSK=FLASE                                                                                                  # set isSK to FALSE to start with until loop

until [[ ${isFamily} == TRUE ]] && [[ ${isSK} == TRUE ]]                                                    # until $isSFamily and $isSK are TRUE, do ...
    do
        taxTemp=`grep "^\<${taxid}\>" ${taxdir}/nodes.dmp`                                                  # retrieve the line with information on the current TAX-ID from the nodes.dmp database
        if (( `echo ${taxTemp} | grep -wc "\<family\>"` > 0 ))                                              # check how often you can grep "family" for current Tax-id and if it more 0 times, then ...
            then                                                                                            
                famtax=`echo ${taxTemp} | cut -d "|" -f1 | tr -d [:blank:]`                                 # get the famtax in the first column
                taxid=`echo ${taxTemp} | cut -d "|" -f2 | tr -d [:blank:]`                                  # set taxid to taxid in second column
                isFamily=TRUE                                                                               # and set isFamily to TRUE (so one part of until loop is already TRUE
        elif (( `echo ${taxTemp} | grep -wc "\<superkingdom\>"` > 0 ))                                      # if you can grep "superkingdom" more than 0 times, then ...
            then
                sktax=`echo ${taxTemp} | cut -d "|" -f1 | tr -d [:blank:]`                                  # get sktax from first column of taxTemp
                taxid=`echo ${taxTemp} | cut -d "|" -f2 | tr -d [:blank:]`                                  # set taxid to second column of taxTemp
                isSK=TRUE                                                                                   # and set isSK to TRUE (so that part for until loop is TRUE)
        elif (( ${taxid} == 1 ))                                                                              # else if tax-id is 1 (reached root), then ...
            then
                famtax=NA                                                                                   # set famtax to NA (no family found in tax-tree)
                isFamily=TRUE                                                                               # and set isFamily=TRUE (otherwise until loop cannot be left)
                isSK=TRUE
        elif [[ -z $taxTemp ]]                                                                              # if taxTemp is empty is empty, then ...
            then
                break                                                                                       # leave loop
        else                                                                                                # if you can do neither of the above mentioned ...
            taxid=`echo ${taxTemp} | cut -d "|" -f2 | tr -d [:blank:]`                                      # then retrieve parent-TAX-ID of taxid (parent in column 2 -> cut -f2)
        fi
    done                    
#gi=`echo $line | cut -d " " -f1`
if [[ $i == tid_part-* ]]                                                                                   # if function was started with part-files, then ...
    then
        echo -ne "${line}\t$tid\t${famtax}\t${sktax}\n" >> ${i}-tid-fam-sk-tax.txt                          # add taxids in subfiles, 
    else                                                                                                    # else ...
        echo -ne "${line}\t$tid\t${famtax}\t${sktax}\n" >> tid-fam-sk-tax.txt                               # paste line (gi and tid to start with), and found famtax and sktax and write into one file
fi
        }

