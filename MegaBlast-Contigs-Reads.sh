#!/bin/bash

#    MegaBlast-Contigs-Reads.sh - part of RIEMS - Reliable Information Extraction from Metagenomic Sequence datasets
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


. ${installdir}/Config.txt                                                                                              # <- Path of config file (unchanged)
. ${installdir}/functions_megablast_contig_reads.sh


multiblastdir=${arbeitsvz}/MultiBlast
cd ${arbeitsvz}/MultiBlast                                                                                              # cd: change directory to working directory
case "$mbcr" in                                                                                                         # get $mbcr and decide if contig or read
    contig)                                                                                                                    # if $mbcr is contig
        #echo -ne "\nRemaining reads will be assigned to assembled contigs... "
        #. ${installdir}/MegaBlast-Mapping.sh
        assigning=Assembly
        echo -ne "\n$(date) --> MEGABLAST of Contigs vs nucleotide database... "                                                     # user info
        cp ${arbeitsvz}/assembly/454AllContigs.fna ${multiblastdir}/454AllContigs.fna                                   # copy contig-file to working directory
        cut -d " " -f 1 ${arbeitsvz}/assembly/ContigInfo.txt | tr -d ">" > ${multiblastdir}/ContigInfo.txt              # cut the first column (contig names) of contigInfo and delete ">"

        START_4=$(date +%s)
        #${blastdir}/blastn -db $ntdb -num_threads $threads -evalue $evalue -max_hsps 1 -max_target_seqs 1 -negative_gilist ${installdir}/gi_exclude.txt -out ${multiblastdir}/Blast-Hits.txt -query ${multiblastdir}/454AllContigs.fna -outfmt '6 qseqid sseqid pident staxid qcovs'    # blast contigs and get query-id and hit-id tab-seperated
        ${blastdir}/blastn -db $ntdb -num_threads $threads -evalue $evalue -max_hsps 1 -max_target_seqs 1 -out ${multiblastdir}/Blast-Hits.txt -query ${multiblastdir}/454AllContigs.fna -outfmt '6 qseqid sseqid pident staxid qcovs'    # blast contigs and get query-id and hit-id tab-seperated
        END_4=$(date +%s)
        DIFF_4=$(( $END_4 - $START_4 ))
        echo -ne "\nTIMING\t${DIFF_4}\tMegablast of the assembled contigs during subset analysis"

        START_5=$(date +%s)
        echo -ne "\n[$(date)] Post processing contig megablast results carried out during subset analysis"
        awk '$5 > 70' ${multiblastdir}/Blast-Hits.txt > ${multiblastdir}/sorted-Hits.txt
        cut -f1-4 ${multiblastdir}/sorted-Hits.txt > ${multiblastdir}/Blast-Hits.txt

        #awk '{sub("N/A", "NA", $4);print}' ${multiblastdir}/Blast-Hits.txt > ${multiblastdir}/non-NA.txt     # Substitute N/A with NA
        awk '$4 != "N/A"' ${multiblastdir}/Blast-Hits.txt > ${multiblastdir}/non-NA.txt
        mv ${multiblastdir}/non-NA.txt ${multiblastdir}/Blast-Hits.txt

        if ! [ -s ${multiblastdir}/Blast-Hits.txt ]                                                                     # if Blast-Hits.txt is empty, then ...
            then
                rm ${multiblastdir}/454AllContigs.fna                                                                   # delete 454AllContigs.fna
                echo No hits found!                                                                                     # user info
                let rndseqextract=rndseqextract+1                                                                       # add 1 to rndseqextract-counter
                break                                                                                                   # and leave case-loop
        fi

        echo -ne "\n[$(date)] Combining contig information with contig Megablast results."
        mv ${multiblastdir}/Blast-Hits.txt ${multiblastdir}/tmp.txt                                                                                       # move Blast-Hits.txt to temporary file
        while read line                                                                                                 # while reading ContigInfo.txt, do...
            do
                grep ${line} ${multiblastdir}/tmp.txt >> ${multiblastdir}/Blast-Hits.txt                                                                  # grep the contig name in temp-file and write it to Blast-Hits.txt do change order of hits so contig with most reads is first
            done < ${multiblastdir}/ContigInfo.txt
        END_5=$(date +%s)
        DIFF_5=$(( $END_5 - $START_5 ))
        echo -ne "\nTIMING\t${DIFF_5}\tPost processing subset analysis contig megablast results"

        START_6=$(date +%s)
        echo -ne "\n[$(date)] Integrating taxonomic information with reads that formed contigs during subset analysis."
        get_tid                                                                                                         # get the tax-id for each hit by gi
        contigs_to_tid
        ###### ATTEMPT GNU PARALLEL HERE
        taxClassificationInfo
        #split -l 1 uniq-tids.txt tiduniq-part-
        #parallel -a uniq-tids.txt -k -l 1 -j0 taxClassificationInfo {}           # Obtain information on taxonomic IDs for assignment
        #rm tiduniq-part*
        #sort ${arbeitsvz}/identifiedTaxClassifications.txt | uniq > ${arbeitsvz}/tmpTaxClass.txt
        #mv ${arbeitsvz}/tmpTaxClass.txt ${arbeitsvz}/identifiedTaxClassifications.txt
        ################################
        assign_reads
        END_6=$(date +%s)
        DIFF_6=$(( $END_6 - $START_6 ))
        echo -ne "\nTIMING\t${DIFF_6}\tIntegrating taxonomic information with reads during subset analysis"

        START_7=$(date +%s)
        echo -ne "\n[$(date)] Running mappings of all remaining reads against contig megablast identified organisms."
        if [ -s ${multiblastdir}/uniq-tids.txt ]                                                                        # if uniq-tids.txt exists and is not empty, then ...
            then
                method=Mapping
                asign_to_tid                                                                                            # assign reads to tax-id
        fi                                                                                                              # end of case "contig"
        END_7=$(date +%s)
        DIFF_7=$(( $END_7 - $START_7 ))
        echo -ne "\nTIMING\t${DIFF_7}\tRunning mappings of remaining reads against contig megablast identified organisms"
        rm -r ${arbeitsvz}/MultiBlast/* ;;

    read)
        assigning=Megablast
        echo -ne "\n$(date) --> MEGABLAST of random Reads vs nucleotide database... "                                               # user info
        if [[ ! -d ${multiblastdir}/RndReads ]]                                                                         # if a RndReads directory does not exist already, then ...
            then
                mkdir ${multiblastdir}/RndReads                # creat the directory
                multiblastrnddir=${multiblastdir}/RndReads
        fi

        START_9=$(date +%s)
        #${blastdir}/blastn -db $ntdb -num_threads $threads -max_hsps 1 -evalue $evalue -max_target_seqs 1 -negative_gilist ${installdir}/gi_exclude.txt -out ${arbeitsvz}/MultiBlast/RndReads/Blast-Hits.txt -query ${arbeitsvz}/pick.fna -outfmt '6 qseqid sseqid pident staxid qcovs'   # blast contigs and get query-id and hit-id tab-seperated
        ${blastdir}/blastn -db $ntdb -num_threads $threads -max_hsps 1 -evalue $evalue -max_target_seqs 1 -out ${arbeitsvz}/MultiBlast/RndReads/Blast-Hits.txt -query ${arbeitsvz}/pick.fna -outfmt '6 qseqid sseqid pident staxid qcovs'   # blast contigs and get query-id and hit-id tab-seperated
        END_9=$(date +%s)
        DIFF_9=$(( $END_9 - $START_9 ))
        echo -ne "\nTIMING\t${DIFF_9}\tMegablast of random reads against the nt database"

        START_10=$(date +%s)
        echo -ne "\n[$(date)] Post processing random reads megablast results."
        awk '$5 > 70' ${arbeitsvz}/MultiBlast/RndReads/Blast-Hits.txt > ${arbeitsvz}/MultiBlast/RndReads/sorted-Hits.txt
        cut -f1-4 ${arbeitsvz}/MultiBlast/RndReads/sorted-Hits.txt > ${arbeitsvz}/MultiBlast/RndReads/Blast-Hits.txt

        awk '$4 != "N/A"' ${arbeitsvz}/MultiBlast/RndReads/Blast-Hits.txt > ${arbeitsvz}/MultiBlast/RndReads/non-NA.txt     # Remove N/A hits
        mv ${arbeitsvz}/MultiBlast/RndReads/non-NA.txt ${arbeitsvz}/MultiBlast/RndReads/Blast-Hits.txt
        END_10=$(date +%s)
        DIFF_10=$(( $END_10 - $START_10 ))
        echo -ne "\nTIMING\t${DIFF_10}\tPost-processing random reads megablast results"

        # megablast of random reads vs nucleotid database
        cd ${arbeitsvz}/MultiBlast/RndReads
        if [ -s ${multiblastrnddir}/Blast-Hits.txt ]                                                                    # if Blast-Hits.txt exists and is not empty, then ...
            then
                START_14=$(date +%s)
                echo -ne "\n[$(date)] Integrating taxonomic information to blast results and reads."
                get_tid_read                                                                                            # get the tax-id for each hit by gi
                END_14=$(date +%s)
                DIFF_14=$(( $END_14 - $START_14 ))
                echo -ne "\nTIMING\t${DIFF_14}\tObtaining taxonomic ID for each hit by GI accession - involves many different steps"
            else                                                                                                        # else...
                echo -e "No hits found.\n"                                                                              # user info
                let rndseqextract=rndseqextract+1                                                                       # add 1 to rndseqextract-counter
        fi
        rm -r ${arbeitsvz}/MultiBlast/* ;;                                                                             # end of case "read"
esac