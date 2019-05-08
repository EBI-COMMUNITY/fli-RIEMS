#    Blast_functions.sh - part of RIEMS - Reliable Information Extraction from Metagenomic Sequence datasets
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

function sortBlastResults()
{
awk '$9 > 70' Hits.txt > sorted-Hits.txt
mv sorted-Hits.txt Hits.txt

awk '$8 != "N/A"' Hits.txt > non-NA.txt     # Remove hits that did not present with a taxonomic ID
mv non-NA.txt Hits.txt
}

function BlastHitTaxid() 
{
    cut -f8 Hits.txt > tid.txt                                                                              # cut 8th column of Hits.txt to get all taxids
    sort tid.txt | uniq > uniq-tids.txt                                                                     # sort them and make them uniq
    tid_num=`wc -l < $blastordner/uniq-tids.txt`                                                            # count numbers of taxids detected
    ((split_threads=$threads*4))                                                                            # set threads for splitting to 4 times threads
    if [[ $tid_num -gt $split_threads ]]                                                                    # if number of taxids is greater than split_threads, then ...
        then
            ((p=$tid_num/$split_threads))                                                                   # calculate the number of lines to split taxids
            split -l $p uniq-tids.txt tid_part-                                                             # and split them by $p (calculated)
        else                                                                                                # else, ...
            split -l 1 uniq-tids.txt tid_part-                                                              # write one taxid per file
    fi
    
    if [[ $method == Blastn_vs_Organism ]]                                                                  # if the method was Blast vs Organisms, then ...
        then
            cut -f2-4 ${arbeitsvz}/asignedAcc.txt | sort | uniq > skfamtax.txt                              # cut column 2-4 of assigned Reads to get all previous identified taxids (faster)
            for i in tid_part-*                                                                             # for each taxid-file, do ...
                do
                    while read line                                                                         # while reading each file per line, do ....
                        do 
                            tid=`echo $line`                                                                # assign line to tid
                            #gi=`echo $line | cut -d " " -f1`
                            get_species                                                                     # and get species-tid for tid (TaxidDetermination.sh)
                            sktax=`grep -w "${tid}$" $blastordner/skfamtax.txt | cut -f 1`                  # grep species-tid in previous found taxids and get sktax (1st column)
                            famtax=`grep -w "${tid}$" $blastordner/skfamtax.txt | cut -f 2`                 # and famtax (2nd column)
                            echo -ne "${line}\t${tid}\t${famtax}\t${sktax}\n" >> ${i}-tid-fam-sk-tax.txt    # write results tabseperated to file
                        done < $i &                                                                         # Fuktionsaufruf FamSkTaxDetermination + thread verteilung
                done ; wait
            cat tid_part-*-tid-fam-sk-tax.txt >> tid-fam-sk-tax.txt                                         # cat all sk_fam_tax files
        else                                                                                                # for any other method, ...
            for i in tid_part-*                                                                             # for each tid file, do ...
                do
                    while read line                                                                         # while reading each file per line, do ....
                    do 
                        tid=`echo $line`                                                                    # assign line to tid
                        get_species                                                                         # get species-tid for tid (TaxidDetermination.sh)
                        FamSkTaxDetermination $tid                                                          # and determine sk and famtax (TaxidDetermination.sh)
                    done < $i &                                                                             # Fuktionsaufruf FamSkTaxDetermination + thread verteilung
                done ; wait
                cat tid_part-*-tid-fam-sk-tax.txt >> tid-fam-sk-tax.txt
    fi
    
    (while read line 
        do                                                                                                  # line: Gi's
            echo -e "`grep "^\<${line}\>" $blastordner/tid-fam-sk-tax.txt 2>/dev/null `" >> tid-fam-sk-tax-ununiq.txt
        done < tid.txt & wait)
    sort -n tid-fam-sk-tax-ununiq.txt | cut -f2- > tmp.txt                                                  # numerische Sortierung nach der 1. Spalte (gis)
    paste AllHits.txt tmp.txt > AllHits-TaxID-${rank}TaxID.txt                                              # spaltenweise Verknüpfung der beiden Dateien
    
    sort AllHits-TaxID-${rank}TaxID.txt >> Blast-Hits.txt                                                   # Sortierung anhand der Accno -> Hits, Partial wieder untereinander   
}