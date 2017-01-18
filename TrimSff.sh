#    TrimSff.sh - part of RIEMS - Reliable Information Extraction from Metagenomic Sequence datasets
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

######################################################################### Quality and Adapter Trimming #######################################################################

cd ${arbeitsvz}                                                                                    

echo -ne "\nRead quality trimming... "                                                                  # user info

i=0
while (( $i < 1 ))
    do
        if [[ -n ${unspecified_input} ]]
            then
                ${gsflxdir}/runMapping -force -m -cpu ${threads} -ud -np -mi 100 -ml 100 -no -tr -o ${arbeitsvz}/D-Mapping ${referenz_trimming} ${unspecified_input[@]} 1>/dev/null 
        fi #Quality Mapping against dummy-reference without adapter trimming
        
        if [[ -n ${illumina_input} ]]    
            then
                ${gsflxdir}/runMapping -force -m -cpu ${threads} -ud -np -no -tr -mi 100 -ml 100 -vt ${il_adapter} -o ${arbeitsvz}/Dil-Mapping ${referenz_trimming} ${illumina_input[@]} 1>/dev/null  # Dummy-Mapping Illumina
        fi #Quality Mapping of Illumina data input with adapter trimming
        
        if [[ -n ${ionTorrent_input} ]]
            then
                ${gsflxdir}/runMapping -force -m -cpu ${threads} -ud -np -no -tr -mi 100 -ml 100 -vt ${it_adapter} -o ${arbeitsvz}/Dit-Mapping ${referenz_trimming} ${ionTorrent_input[@]} 1>/dev/null # Dummy-Mapping Iontorrent
        fi #Quality Mapping of IonTorrent data input with adapter trimming
        let i=i+1
    done ; wait

for i in D*-Mapping
    do
        cat ${i}/454TrimmedReads.fastq >> ${arbeitsvz}/TrimmedReads.fastq
        cat ${i}/454TrimStatus.txt >> ${arbeitsvz}/TrimStatus.txt
        cat ${i}/454ReadStatus.txt >> ${arbeitsvz}/ReadStatus.txt
done                                                                                                    # cat of fastq-files, trimstatus and readstatus from all mappings for further metagenomic analysis
  
cut -f3 ${arbeitsvz}/ReadStatus.txt | sort | uniq > ${arbeitsvz}/restAcc.txt                            # get list of all read-accession with high quality
grep -vw "Mapped" ${arbeitsvz}/restAcc.txt | grep -v "Accuracy" > ${arbeitsvz}/tmp1.txt                 # | remove header
mv ${arbeitsvz}/tmp1.txt ${arbeitsvz}/restAcc.txt                                                       # V and move file
cut -f2 ${arbeitsvz}/TrimStatus.txt | grep -v "Accno" > ${arbeitsvz}/allAcc.txt                         # get list of all input read-accessions (includes low quality reads)
grep -vw "^Accno" ${arbeitsvz}/allAcc.txt > ${arbeitsvz}/tmp1.txt                                       # and remove header
sort ${arbeitsvz}/tmp1.txt | uniq > ${arbeitsvz}/allAcc.txt
rm ${arbeitsvz}/tmp1.txt 


comm -3 ${arbeitsvz}/restAcc.txt ${arbeitsvz}/allAcc.txt | sed -e 's/^[ \t]*//' | grep -vw "Mapped"  | sed '1d' > ${arbeitsvz}/trimmedAcc.txt 
                                                                                                        # compare restAcc and allAcc to get list (-3) of all read-accession that were low quality
                        
i=`wc -l ${arbeitsvz}/trimmedAcc.txt | cut -d " " -f1`                                                  # get number of low qualtity reads
if (( $i > 0 ))                                                                                         # and prepare data-output of asignedAcc.txt
    then
        printf 'low quality\n%.0s' $(seq 1 $i) > ${arbeitsvz}/tmp1.txt                                  # write "low quality" for each low quality read
        printf '\t\t\t\t\t\n%.0s' $(seq 1 $i) > ${arbeitsvz}/tmp3.txt                                   # columns have to be added to extend list with further results of mapping and blast
        paste ${arbeitsvz}/tmp3.txt ${arbeitsvz}/trimmedAcc.txt ${arbeitsvz}/tmp1.txt > ${arbeitsvz}/tmp2.txt   #paste temporary files and trimmed Accnos
        mv ${arbeitsvz}/tmp2.txt ${arbeitsvz}/asignedAcc.txt                                            # move tmp2 to asignedAcc.txt
        rm ${arbeitsvz}/tmp1.txt ${arbeitsvz}/tmp3.txt                                                  # and delete temporary files
fi

echo "finished"                                                                                         # user info






