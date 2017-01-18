#    TaxID_exclude.sh - part of RIEMS - Reliable Information Extraction from Metagenomic Sequence datasets
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
#this subscript of the RIEMS workflow generates a gi-list to exclude during the analysis
#by default following Tax-ID are excluded: 48479, 28384, 12908

. ${installdir}/Config.txt                                                                                                                          # <- path to Contig-file

excludeTax=(48479 28384 12908)                                                                                                                      # List of Tax-IDs to exclude... the List might be appended
out=${installdir}/taxtmp                                                                                                                            # set directory for output
mkdir -p ${out}                                                                                                                                     # make the output directory
cd ${out}                                                                                                                                           # change working directory to output
for i in ${excludeTax[*]}                                                                                                                           # for each taxid in the exclude list, do...
    do  
        echo $i                                                                                                                                     # give user info about current taxid
        parentTax=$i                                                                                                                                # save current taxid as parentTax
        grep "^\<${parentTax}\>" ${fasttaxid}/parentChildTaxid/parentChildTaxid-`echo $parentTax | cut -c1-2`.dmp | cut -f2 >> ${out}/tmptaxchild-${i}.txt    # get all ChildTax-ids assigned to parentTax, only grep in file were tax-ids begin with first two characters of parenttax, get 2nd column (childTax)                        
        cat ${out}/tmptaxchild-${i}.txt >> ${out}/taxchild-${i}.txt                                                                                 # cat temporary taxchild to taxchild
        
        until ! [ -s ${out}/tmptaxchild-${i}.txt ]                                                                                                  # until the tmptaxchild.txt is empty, do ...
            do
                while read line                                                                                                                     # while reading tmptaxchild.txt, do ...
                    do
                        m=`echo $line | cut -c1-2`                                                                                                  # save first two characters of taxid in m
                        grep "^\<${line}\>" ${fasttaxid}/parentChildTaxid/parentChildTaxid-`echo $line | cut -c1-2`.dmp | cut -f2 >> ${out}/tmp2taxchild-${i}.txt   # grep grep child-tax-ids in respective parentChildTaxid*.tmp, cut 2nd column (child-tax) and save in 2nd temporary file
                done < ${out}/tmptaxchild-${i}.txt 1>/dev/null 2>&1 & wait                                                                          # in background to grep all at once and wait til they finished
                mv ${out}/tmp2taxchild-${i}.txt ${out}/tmptaxchild-${i}.txt                                                                         # move 2nd to 1st temporary file (now tmptaxchild is not empty)
                cat ${out}/tmptaxchild-${i}.txt >> ${out}/taxchild-${i}.txt                                                                         # cat tmptaxchild to taxchild.txt
        done
        echo $parentTax | cat >> ${out}/taxchild-${i}.txt                                                                                           # finally add parentTax to taxchild
        cat ${out}/taxchild-${i}.txt >> ${installdir}/exclude_tid.txt                                                                               # cat all taxchilds of each given taxid in one file
        z=`wc -l < ${out}/taxchild-${i}.txt`                                                                                                        # count number of taxids
        if (( $z > 23 ))                                                                                                                            # if number is greater than 23, then ..
            then        
                (( p=$z/23 ))                                                                                                                       # calculate p by dividing number of taxids by number of cores available
                split -l $p ${out}/taxchild-${i}.txt ${out}/part-${i}- ; wait                                                                       # split childtax.txt by p, else ...
            else        
                cp ${out}/taxchild-${i}.txt ${out}/part-${i}-aa                                                                                     # move taxid-file to part-file
        fi      
        for j in part-${i}-*                                                                                                                        # for each part-* file, do ...
            do      
                while read line                                                                                                                     # while reading part-* file, do ...
                    do
                        k=`echo $line | cut -c1-2`                                                                                                  # save first two characters of cutrrent line (taxid)
                        grep "^\<${line}\>" ${fasttaxid}/taxid_gi_nucl/taxid_gi_nucl-`echo $line | cut -c1-2`.dmp | cut -d " " -f2 >> ${out}/gi_${i}-${j}.txt  # grep taxid and associated gi in taxid_gi_nucl (only grep in file where taxids start with first to characters of $line), cut 2nd column (gi)
                    done < ${j} &  1>/dev/null 2>&1                                                                                                 # in background to start 23 processes at once
        done ; wait 1>/dev/null 2>&1                                                                                                                # wait to finish processes
        wait
        cat gi_*-part-* > gi_TID-${i}.txt                                                                                                           # cat all sub-files to one file
        #rm tmptaxchild.txt ; rm taxchild.txt ; rm *part*
    done 
    
cat gi_TID* > ${installdir}/gi_exclude.txt                                                                                                          # cat gis for each taxid into one file
cd ${arbeitsvz}                                                                                                                                     # change working directory
rm -r ${installdir}/taxtmp                                                                                                                          # remove folder