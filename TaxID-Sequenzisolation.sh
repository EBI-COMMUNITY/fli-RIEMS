#!/bin/bash

#    TaxID-Sequenzisolation.sh - part of RIEMS - Reliable Information Extraction from MEtagenomic Sequence datasets
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




####################################################### Voreinstellungen bei direkter Benutzung ##########################################################

. ${installdir}/Config.txt                                                                                                                  # <- Path to config file (unchanged)
. ${installdir}/TaxidDetermination.sh                                                                                                       # Path to functions

out=${installdir}/reftmp_$USER                                                                                                                    # define out-directory

################################################################## Sequenzisolierung #####################################################################

if [[ ! -d ${out} ]]
    then
        mkdir -p ${out}; chmod -R 777 ${out}                                                                                                           # creat out-directory
fi
# code in the while-loop will determine the taxonomy level for the respective TAX-ID until the species level is reached; therefore, repeat the loop until the species level was reached

echo -ne "\n[$(date)] Obtaining scientific name for ${parentTax} and downloading genome sequences."
name=`grep "^\<$parentTax\>" ${taxdir}/names.dmp | grep '\<scientific name\>' | cut -f3 | tr -s [:blank:] "-" | tr -d [=\'=],[=\(=],[=\)=],[=\[=],[=\]=] | tr -s [=\/=] "." | tr -d "." | sed 's/[[:punct:]]/-/g'`    # get name for parentTax-id
if [ -z "$name" ]                                                                                                                           # if no name exists, then ...
    then                            
        echo -ne "\nNo Name to TaxID $parentTax found. No sequences will be downloaded.\n"                                                          # user info and no further try to download sequences 
    else
        if [ -s ${refseqdir}/TID-${parentTax}_*.fna ]                                                                                               # if taxid is already available in RefSequenzen, then ...
            then
                echo -ne "\n\n$name was identified. \nSequences for $name are already available as reference sequences \n"                          # user info
            else                                                                                                                                    # else ...
                echo -ne "\n[$(date)] ${name} was identified, attempting to download sequences as references."
                grep "^\<${parentTax}\>" ${fasttaxid}/parentChildTaxid/parentChildTaxid-`echo $parentTax | cut -c1-2`.dmp | cut -f2 >> ${out}/tmptaxchild.txt # get all ChildTax-ids assigned to parentTax, only grep in file were tax-ids begin with first two characters of parenttax, get 2nd column (childTax)
                cat ${out}/tmptaxchild.txt >> ${out}/taxchild.txt                                                                                   # cat temporary taxchild to taxchild

                START_TAXTREEINFO=$(date +%s)
                echo -ne "\n[$(date)] Obtaining taxonomic tree of information"
                until ! [ -s ${out}/tmptaxchild.txt ]                                                                                               # until the tmptaxchild.txt is empty, do ...
                    do
                        while read line                                                                                                             # while reading tmptaxchild.txt, do ...
                            do
                                m=`echo $line | cut -c1-2`                                                                                          # get the first 2 characters of the taxid
                                grep "^\<${line}\>" ${fasttaxid}/parentChildTaxid/parentChildTaxid-`echo $line | cut -c1-2`.dmp | cut -f2 >> ${out}/tmp2taxchild.txt   # grep grep child-tax-ids in respective parentChildTaxid*.tmp, cut 2nd column (child-tax) and save in 2nd temporary file
                            done < ${out}/tmptaxchild.txt & wait                                                                   # in background to grep all at once and wait til they finished
                        mv ${out}/tmp2taxchild.txt ${out}/tmptaxchild.txt                                                                           # move 2nd to 1st temporary file (now tmptaxchild is not empty)
                        cat ${out}/tmptaxchild.txt >> ${out}/taxchild.txt                                                                           # cat tmptaxchild to taxchild.txt
                    done
                echo $parentTax | cat >> ${out}/taxchild.txt                                                                                        # finally add parentTax to taxchild
                cd ${out}
                END_TAXTREEINFO=$(date +%s)
                DIFF_TAXTREEINFO=$(( $END_TAXTREEINFO - $START_TAXTREEINFO ))
                echo -ne "\nTIMING\t${DIFF_TAXTREEINFO}\tObtaining taxonomic tree of information"

                # Potential issue with this if statement???
                if ! [ -s ${refseqdir}/TID-${parentTax}_*.fna ]                                                                                     # if the referenz does not already exist in Refseq, then ...
                    then
                        START_ORGDOWN=$(date +%s)
                        echo -ne "\n[`date`] $name was identified. \nSequences for $name will be downloaded"
                        z=`wc -l < ${out}/taxchild.txt`                                                                                             # count number of taxids
                        if (( $z > $threadsSplit ))                                                                                                 # if $z is bigger thant the set threads, then ...
                            then 
                                (( p=$z/$threadsSplit ))                                                                                            # calculate p by dividing number of taxids by number of cores available
                                split -l $p ${out}/taxchild.txt ${out}/part- ; wait 1>/dev/null 2>&1                                                                # split childtax.txt by p
                            else
                                split -l 1 ${out}/taxchild.txt ${out}/part- ; wait 1>/dev/null 2>&1                                                                              # just move file to part-aa
                        fi
                        for i in part-*                                                                                                             # for each part-* file, do ...
                            do 
                                while read line                                                                                                     # while reading part-* file, do ...
                                    do
                                        k=`echo $line | cut -c1-2`                                                                                  # save first to characters of taxid in k (for grep)
                                        grep "^\<${line}\>" ${fasttaxid}/taxid_gi_nucl/taxid_gi_nucl-`echo $line | cut -c1-2`.dmp | cut -d " " -f2 >> ${out}/gi_${i}.txt  # grep taxid and associated gi in taxid_gi_nucl (only grep in file where taxids start with first to characters of $line), cut 2nd column (gi)
                                    done < ${i}  &                                                                                      # in background to start 23 processes at once
                            done ; wait 1>/dev/null 2>&1                                                                                                # and wait til all of them finish
                        rm ${out}/part-*                                                                                                            # remove all part files
                        cat gi_* > ${out}/gis.txt                                                                                                   # cat all gi files
                        rm gi_*                                                                                                                     # remove subset gi files
                        z=`wc -l < ${out}/gis.txt`                                                                                                  # count number of all gis found
                        if (( $z > $threadsSplit ))                                                                                                 # if number of gis is bigger than number of set threads, then ...
                            then 
                                (( p=$z/$threadsSplit ))                                                                                            # calculate p by dividing number of taxids by number of cores available
                                split -l $p ${out}/gis.txt ${out}/gi_part- ; wait  1>/dev/null 2>&1                                                 # split childtax.txt by p
                            else
                                split -l 1 ${out}/gis.txt ${out}/gi_part- ; wait  1>/dev/null 2>&1                                                  # else split by line
                        fi
                        (for i in ${out}/gi_* ; do GiDownload $i & done 1>/dev/null 2>&1; wait)                                                     # download sequences for gis found (GiDownload function in TaxidDetermination.sh)
                        cat ${out}/gi_*.fna > ${out}/TID-${parentTax}.fna                                                                           # get all subset Sequence files
                        rm ${out}/gi_*                                                                                                   # and remove all subset sequence files
                        rm ${out}/taxchild.txt ; rm ${out}/tmptaxchild.txt                                                                          # remove taxchild and tmptax
                        if [[ -s ${out}/TID-${parentTax}.fna ]]
                            then
                                START_X2KILL=$(date +%s)
                                DoubbleKill                                                                                                                 # doubbleKill function to remove duplicated gis (TaxidDetermination.sh)
                                END_X2KILL=$(date +%s)
                                DIFF_X2KILL=$(( $END_X2KILL - $START_X2KILL ))
                                echo -ne "\nTIMING\t${DIFF_X2KILL}\tRemoving duplicated sequences from the download"
                            else
                                echo -ne "\nNo Sequences found for TaxID $parentTax $name.\n"
                        fi
                        END_ORGDOWN=$(date +%s)
                        DIFF_ORGDOWN=$(( $END_ORGDOWN - $START_ORGDOWN ))
                        echo -ne "\nTIMING\t${DIFF_ORGDOWN}\tDownloading sequences for organism (includes other steps)"
                    fi
                if [ -s ${refseqdir}/TID-${parentTax}_$name.fna ]                                                                                   # if the TID-$parentTax exists and is not empty, then ...
                    then 
                        echo                                                                                                                        # print empty string
                    else                                                                                                                            # else, ...
                        echo -e "\nNo sequences found\n"                                                                                              # user info
                fi 
        fi    
fi

