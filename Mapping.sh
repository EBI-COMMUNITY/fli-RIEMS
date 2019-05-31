#    Mapping.sh - part of RIEMS - Reliable Information Extraction from Metagenomic Sequence datasets
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

#################################################################### Voreinstellungen ####################################################################

. ${installdir}/Config.txt                                                                                                    # <- Path of config file (unchanged)

referenz=${refseqdir}/TID-${parentTax}_*.fna                                                                              # Default setting for referenz

if ls ${referenz} 1> /dev/null 2>&1; then
    ziel=${arbeitsvz}/TID-${parentTax}_$(grep "^\<$tid\>" ${taxdir}/names.dmp | grep '\<scientific name\>' | cut -f3 | tr -s [:blank:] "-" | tr -d [=\'=],[=\(=],[=\)=],[=\[=],[=\]=] | tr -s [=\/=] "." | tr -d "." | sed 's/[[:punct:]]/-/g')
                                                                                                                            # define folder_name
    zielordner=${ziel}/MappingFiles                                                                                         # Default setting for folder

    ######################################################################## Mapping #########################################################################

    START_MAPVSORG=$(date +%s)
    mkdir -p ${ziel}/MappingFiles                                                                                           # create folder
    organism=`basename $referenz`                                                                                           # get name of organism
    echo -ne "\n[$(date)] Running Mapping vs ${organism}... "                                                                         # user info
    num=`wc -l < ${arbeitsvz}/restAcc.txt`                                                                                  # count remaining reads
    if (( num > 10000 ))                                                                                                    # if more than 10000 reads remaining do...
        then
            mem=`free -bt | tail -n 1 | cut -f 5 -d " "`                                                                    # get free memory
            size=`du -b $referenz | cut -f 1`                                                                               # get size of reference
            ((j=$mem/(2*$size)))                                                                                            # calculate how much mappings can be startet in parallel (twice the size because reads have to loaded and software needs memory itself)
            if (( j > $threads ))                                                                                           # if more processes could be started than cores available, then ...
                then
                    j=$threads                                                                                              # set processes to start to cores available
            fi
            if (( size > 100000000 ))
                then
                ((j=$j/4))
                if (( $j == 0 ))
                    then
                        j=1
                fi
            fi
            cd ${zielordner}                                                                                                # change directory
            ((p=$num/$j))                                                                                                   # calculate p by number of remaining reads and core to be used
            split -a 5 -d -l $p ${arbeitsvz}/restAcc.txt part-                                                                      # split file of remaining accno. by p lines

            START_ACTUALMAPPING=$(date +%s)
            ###### ATTEMPT GNU-PARALLEL
            #ls part-* | parallel -j0 "${gsflxdir}/runMapping -o map-{} -no -force -noace -nobam -noinfo -rst 0 -m -n -np -mi $identity -ml 95% -ud -cpu 10 -tr -notrim -fi {} ${referenz} ${arbeitsvz}/TrimmedReads.fastq 1>/dev/null"

            for i in part-*                                                                                                 # for each file do...
                do
                    ${gsflxdir}/runMapping -o map-${i} -no -force -noace -nobam -noinfo -rst 0 -m -n -np -mi $identity -ml 95% -ud -cpu 24 -tr -notrim -fi ${i} ${referenz} ${arbeitsvz}/TrimmedReads.fastq 1>/dev/null &    # map reads to Reference
                done        #runMapping against reference sequences for TaxID ($referenz) with given minimum identity (-mi) and a minimum overlap (-ml) of 95%, only with reads that are included in accession-file (-fi)
            wait   1>/dev/null 2>&1                                                                                     # and wait until all mappings have finished
            ##########################
            END_ACTUALMAPPING=$(date +%s)
            DIFF_ACTUALMAPPING=$(( $END_ACTUALMAPPING - $START_ACTUALMAPPING ))
            echo -ne "\nTIMING\t${DIFF_ACTUALMAPPING}\tCommand for mapping remaining reads against ${organism}"

            START_MAPPROC=$(date +%s)
            echo -ne "\n[$(date)] Post-processing mapping results..."
            for i in map-part-*                                                                                             # for each mapping, do ....
                do
                    cat ${i}/454ReadStatus.txt >> 454ReadStatus.txt &                                                       # get the ReadStatus and cat them together
                done; wait
            END_MAPPROC=$(date +%s)
            DIFF_MAPPROC=$(( $END_MAPPROC - $START_MAPPROC ))
            echo -ne "\nTIMING\t${DIFF_MAPPROC}\tCreating one file of combined mapping results"

            rm -r ${zielordner}/map-* ; rm ${zielordner}/part-*
            cd ${arbeitsvz}                                                                                                 # change directory
        else
            START_LESSMAP=$(date +%s)
            echo -ne "\n --> Less than 10000 reads remaining to map. "
            ${gsflxdir}/runMapping -force -noace -no -nobam -noinfo -m -n -np -mi $identity -ml 95% -rst 0 -ud -cpu $threads -notrim -tr -o $zielordner -fi ${arbeitsvz}/restAcc.txt $referenz ${arbeitsvz}/TrimmedReads.fastq 1>/dev/null  # else map all unassigned reads against reference
            END_LESSMAP=$(date +%s)
            DIFF_LESSMAP=$(( $END_LESSMAP - $START_LESSMAP ))
            echo -ne "\nTIMING\t${DIFF_LESSMAP}\tCommand for mapping remaining < 10000 reads against ${organism}"
    fi                      #runMapping against reference sequences for TaxID ($referenz) with given minimum identity (-mi) and a minimum overlap (-ml) of 95%, only with reads that are included in accession-file (-fi)

    echo $referenz >> ${arbeitsvz}/references.txt                                                                           # write the referenz to a file for further use in Mapping2
    ln -s $referenz ${arbeitsvz}/DB/`basename $referenz`                                                                    # create a link of to reference to the DB-directory
    grep "^>" $referenz | cut -f2 -d "|" >> ${arbeitsvz}/Orgs_gis.txt                                                       # get all gis/accessions from reference for further use in Blastn_vs_Organisms
    END_MAPVSORG=$(date +%s)
    DIFF_MAPVSORG=$(( $END_MAPVSORG - $START_MAPVSORG ))
    echo -ne "\nTIMING\t${DIFF_MAPVSORG}\tMapping reads against an identified organism (includes other steps)"

    START_MAPRET=$(date +%s)
    get_mapping_results                                                                                                     # get mapping results (see functions.sh)
    END_MAPRET=$(date +%s)
    DIFF_MAPRET=$(( $END_MAPRET - $START_MAPRET ))
    echo -ne "\nTIMING\t${DIFF_MAPRET}\tRetrieving mapping results"
else
    echo -ne "\nNo sequence file found for $referenz. Cannot run mapping..."
fi

