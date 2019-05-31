#!/bin/bash
#    Assemnly-ContigBlast.sh - part of RIEMS - Reliable Information Extraction from Metagenomic Sequence datasets
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


. ${installdir}/Config.txt                                                                                                            # <- Path of config file (unchanged)

#chmod a+rxw ${arbeitsvz}/Assembly
outputordner=${arbeitsvz}/Assembly                                                                                                    # Defining working directory

########################################################### functions ####################################################################################

########################################################### Assembly + Read-Contig Zuordnung #############################################################
# Ermittlung der Dateiendung; falls vom User eine sff-datei eingegeben wurde, muss diese noch in eine fna-Datei umgewandelt werden

mkdir -p $outputordner                                                                                                          # creat directory
mkdir -p $outputordner/Assembly
chmod a+rxw $outputordner/Assembly

allaccnos=`wc -l < ${arbeitsvz}/restAcc.txt `                                                                                   # count all unassigned reads
        
echo -ne "\n\n$(date) --> Running ASSEMBLY of all $allaccnos unassigned Reads... \n"                                                        # user info

process=newbler

#cp ${arbeitsvz}/restAcc.txt ${arbeitsvz}/restAcc_`date +"%T"`
START_14=$(date +%s)
echo -ne "Carrying out assembly of unassigned reads: \n"
#echo -ne "${gsflxdir}/runAssembly -rip -o ${outputordner}/Assembly -fi ${arbeitsvz}/restAcc.txt -notrim -tr -force -noinfo -cpu ${threads} -acedir -ml 20 ${arbeitsvz}/TrimmedReads.fastq"
#echo -ne "${gsflxdir}/fnafile -o ${arbeitsvz}/restReads.fasta -i ${arbeitsvz}/restAcc.txt ${arbeitsvz}/TrimmedReads.fasta"
${gsflxdir}/fnafile -o ${arbeitsvz}/restReads.fasta -i ${arbeitsvz}/restAcc.txt ${arbeitsvz}/TrimmedReads.fasta                 # creat fasta-file of all unassigned reads
${gsflxdir}/runAssembly -rip -o ${outputordner}/Assembly -fi ${arbeitsvz}/restAcc.txt -notrim -tr -force -noinfo -cpu ${threads} -acedir -ml 20 ${arbeitsvz}/TrimmedReads.fastq 1>/dev/null & # runAssembly with all unassigned reads
wait
#cp ${outputordner}/Assembly/454AllContigs.fna ${arbeitsvz}/454AllContigs.fna.`date| perl -lne '~s/ //g;print;'`
   # -force: force overwriting ; -acedir: directory with info about read assignment to contigs; -ml min length of read ( in %) to be assembled
   # -ml: wahrscheinlich minimale L�nge eines Reads in % um Bestandteil eines Contigs zu werden
cd $outputordner
END_14=$(date +%s)
DIFF_14=$(( $END_14 - $START_14 ))
echo -ne "\nTIMING\t${DIFF_14}\tAssembly of remaining reads during assembly contig blast"

echo -ne "\n$(date) --> Data processing... "                                                                                                # user info
if [[ -s Assembly/454ReadStatus.txt ]] && [[ -s Assembly/454AllContigs.fna ]]
    then
        START_15=$(date +%s)
        egrep "Assembled|Repeat" Assembly/454ReadStatus.txt > Assembly/assembledAcc.txt                                         # grep readstatus for assembled reads
        egrep -v "Assembled|Repeat|Accno" Assembly/454ReadStatus.txt | cut -f 1 > Assembly/unassembledAcc.txt                   # grep readstatus for all unassembled reads (also excludes header line)
        assembledaccnos=`wc -l < Assembly/assembledAcc.txt`                                                                     # count all assembled reads
        gesamt=`grep -c '>' Assembly/454AllContigs.fna`                                                                         # count number of contigs
        echo "completed & $assembledaccnos Reads assembled into $gesamt Contig(s)."                                             # user info
        END_15=$(date +%s)
        DIFF_15=$(( $END_15 - $START_15 ))
        echo -ne "\nTIMING\t${DIFF_15}\tPost processing of reads assembled into contigs"

################################################################ Megablast Rest vs Contigs ################################################################
        echo -ne "\n\n-----{Megablast Unassembled Reads against Contigs - Start [`date`]}-----"
        START_16=$(date +%s)
        echo -ne "\n[$(date)] Preparing Megablast for all unassembled reads... "                                                             # user info
        grep -v 'Accno' Assembly/unassembledAcc.txt | cut -f1 | sort | uniq > Assembly/All-Accnos.txt                           # get all Accnos sort them and make them uniq
        allaccnos=`wc -l < Assembly/All-Accnos.txt`                                                                             # count number of Reads
        #rm Assembly/All-Accnos.txt
        
        ${gsflxdir}/fnafile -o Assembly/sff/Unassigned-Accnos.fasta -i Assembly/unassembledAcc.txt ${arbeitsvz}/TrimmedReads.fasta   # write all unassembled Reads to fasta-file
        echo -ne "\nCreating Blast database... "                                                                                # user info
        infilesize=`du -sbL Assembly/454AllContigs.fna | cut -f1`                                                               # get size of file
        dbsize=`echo "$infilesize*0.26/${threads}" | bc`                                                                        # divide size to all CPUs
        ${blastdir}/makeblastdb -dbtype nucl -in Assembly/454AllContigs.fna -max_file_sz ${dbsize}                              # and make a blastdb of the file
        END_16=$(date +%s)
        DIFF_16=$(( $END_16 - $START_16 ))
        echo -ne "\nTIMING\t${DIFF_16}\tPreparation for megablast of unassembled reads during assembly contig blast"

        echo -ne "\n[$(date)] Running Megablast... "                                                                                      # user info
        START_17=$(date +%s)
        ${blastdir}/blastn -db Assembly/454AllContigs.fna -num_threads ${threads} -max_hsps 1 -max_target_seqs 1 -out Assembly/Megablast_Unassigned.txt -query Assembly/sff/Unassigned-Accnos.fasta -outfmt '6 pident qseqid sseqid qstart qend qlen evalue '  # wg. Problemen mit BLAST bei culling_limit entfernt:
                # megablast of Reads vs Contigs
        END_17=$(date +%s)
        DIFF_17=$(( $END_17 - $START_17 ))
        echo -ne "\nTIMING\t${DIFF_17}\tMegablast of unassembled reads during assembly contigs blast"

        if [[ -s Assembly/Megablast_Unassigned.txt ]]                                                                           # if Assembly/Megablast_Unassigned.txt exists (Blast was successful and had hits), then ...
            then
                START_18=$(date +%s)
                echo -ne "\n[`date`] Post-processing unassembled reads megablast results."
                #cp Assembly/Megablast_Unassigned.txt ${arbeitsvz}/Megablast_Unassigned.txt.`date| perl -lne '~s/ //g;print;'`
                #echo -ne "\rRunning Megablast... data processing... "                                                          # user info
                gesamt=`wc -l < Assembly/Megablast_Unassigned.txt`                                                              # count number of Hits
                while read line                                                                                                 # while reading Megablast_Unassigned.txt, do ...
                    do                                                                                                          
                        #echo -ne "\rRunning Megablast... data processing... "                                                  # user info
                        grep "$line" Assembly/Megablast_Unassigned.txt > Assembly/tmp.txt                                       # write line to tmp-file
                        x1=$(echo $line | cut -d " " -f4)                                                                       # get start of hit
                        x2=$(echo $line | cut -d " " -f5)                                                                       # end of hit
                        len=$(echo $line | cut -d " " -f6)                                                                      # read length
                        ((diff=(${x2}-${x1})*100/${len}))                                                                       # calculate how much of the read matched in blast
                        if (( $diff >= 60 ))                                                                                    # if more than 60% were aligned, then ...
                            then                                                                                                
                                cat Assembly/tmp.txt >> Assembly/Megablast_Assigned.txt                                         # write tmp (containing Hit-Read) to Megablast_Assigned.txt
                        fi
                    done < Assembly/Megablast_Unassigned.txt & wait
                END_18=$(date +%s)
                DIFF_18=$(( $END_18 - $START_18 ))
                echo -ne "\nTIMING\t${DIFF_18}\tProcessing unassembled reads megablast results"
        fi
        echo -ne "\n\n-----{Megablast Unassembled Reads against Contigs - End [`date`]}-----"

###################################################################### Blastn Contigs #####################################################################
        echo -ne "\n\n-----{Megablast Contigs against Nucleotide Database - Start [`date`]}-----"
        # Aufruf des Blast Programms
        echo -ne "\n\nRunning Megablast to identify the Contigs\n"                                                                  # user info
        
        START_19=$(date +%s)
        blastvek=(0 0 tax 0 0)                                                                                                  # set blastvek for Blast.sh
        query=${outputordner}/Assembly/454AllContigs.fna                                                                        # set query for Blast.sh
        blastordner=${outputordner}/Blast                                                                                       # set blastfolder for Blast.sh
        ntblast=megablast                                                                                                       # set which blast for Blast.sh
        referenz=$ntdb                                                                                                          # set reference for Blast.sh
        mkdir -p $blastordner
        cd $blastordner
        . $installdir/Blast.sh                                                                                                  # Blast.sh
        END_19=$(date +%s)
        DIFF_19=$(( $END_19 - $START_19 ))
        echo -ne "\nTIMING\t${DIFF_19}\tMegablast of assembled contigs during assembly contigs blast"

        # Zu den Family Tax-id zugeordneten Contigs aus dem Blast, werden die Accnos herausgesucht
        START_20=$(date +%s)
        echo -ne "\n[$(date)] Final data processing... "                                                                                    # user info
        cd $outputordner                                                                                                        # change back to working directory (changed in Blast.sh)
        cp ${arbeitsvz}/Assembly/Assembly/454AllContigs.fna .                                                                   # copy Contig-file to folder
        sed '1,2d' Blast/Blast-Hits.txt > 454AllContigs-BlastnHits.txt                                                          # move blast results and delete first two rows (header and one empty row)
        cut -f1-3 Assembly/454ReadStatus.txt > Assembly/454ReadInfo.txt                                                         # cut column 1-3 from ReadStatus to avoid two results for one read because assembled in both
        
        grep ">" 454AllContigs.fna > allContigs.txt                                                                             # grep for header in AllContigs to get contig names
        cut -f1 454AllContigs-BlastnHits.txt > assignedContigs.txt                                                              # cut first column of Hits to cet all column that were assigned
        grep -v -f assignedContigs.txt allContigs.txt > unassignedContigs.txt                                                   # grep all except the contigs that were assigned (header) (-f to take a file)
        ${gsflxdir}/fnafile -i unassignedContigs.txt -o unassignedContigs.fna 454AllContigs.fna                                 # write all unassigned Contigs to new file

        # ==================================================== #
        #cp -R ${outputordner}/Blast/ ${arbeitsvz}/AssemblyBlast_`date| perl -lne '~s/ //g;print;'`
        # ==================================================== #
        rm ${outputordner}/Blast/*                                                                                              # remove all input in Blastfolder for next blast
        END_20=$(date +%s)
        DIFF_20=$(( $END_20 - $START_20 ))
        echo -ne "\nTIMING\t${DIFF_20}\tProcessing assembled contigs blastn results"
        echo -ne "\n\n-----{Megablast Contigs against Nucleotide Database - End [`date`]}-----"

        echo -ne "\n\n-----{Blastn of Unassigned Contigs against Nucleotide Database - Start [`date`]}-----"
        echo -ne "\n[$(date)] Running Blastn of unassigned contigs against ntdb..."
        START_21=$(date +%s)
        blastvek=(0 0 tax 0 0)                                                                                                  # set blastvek for Blast.sh
        query=${outputordner}/unassignedContigs.fna                                                                             # set query for Blast.sh
        blastordner=${outputordner}/Blast                                                                                       # set blastfolder for Blast.sh
        ntblast=blastn                                                                                                          # set which blast for Blast.sh
        referenz=$ntdb                                                                                                          # set reference for Blast.sh 
        . $installdir/Blast.sh                                                                                                  # Blast.sh
        END_21=$(date +%s)
        DIFF_21=$(( $END_21 - $START_21 ))
        echo -ne "\nTIMING\t${DIFF_21}\tBlastn of unassigned contigs during assembly contigs blast"

        START_ASSEMBLYPROC=$(date +%s)
        echo -ne "\n[`date`] Post processing unassigned contigs blastn results."
        cd $outputordner
        sed '1,2d' Blast/Blast-Hits.txt | cat >> 454AllContigs-BlastnHits.txt                                                   # move Blast-Hits and delete the first to lines (header and empty)
        
        i=`wc -l < 454AllContigs-BlastnHits.txt`                                                                                # count blast hits (one line one hit)
        if (( $i > 0 ))                                                                                                         # if at least one blast hit, then ...
            then
                START_22=$(date +%s)
                printf "Assembly\n%.0s" $(seq 1 $i) > tmp1.txt                                                                  # print the method as often as blast hits in temporary file
                    
                paste tmp1.txt 454AllContigs-BlastnHits.txt > contig_method.txt                                                 # paste temp file and blast hits to assign method to results
                cut -f 8 454AllContigs-BlastnHits.txt > taxids.txt                                                              # cut 8th column of blast hits to get taxid
                sort taxids.txt | uniq > taxids-uniq.txt
                while read line                                                                                                 # while reading taxids.txt, do ...
                    do 
                        if [[ `grep -w "^${line}" ${taxdir}/names.dmp | grep '\<scientific name\>'` != "" ]]                    # if "scientific name" for taxid can be grepped and is not an empty string, then ...
                            then
                                grep -w "^${line}" ${taxdir}/names.dmp | grep '\<scientific name\>' | cut -f3 | tr -d [=\'=],[=\(=],[=\)=],[=\[=],[=\]=] | tr -s [=\/=] "." | tr -s "-" " " | tr -d "." | sed 's/[[:punct:]]/-/g' >> names.txt   # grep $line in names.dmp, grep "scientific name" (if possible) and cut 3rd column for name
                        elif [[ `grep -w "^${line}" ${taxdir}/names.dmp | grep -m 1 '\<synonym\>'` != "" ]]                                                                                                # else, ...
                            then
                                grep -w "^${line}" ${taxdir}/names.dmp | grep -m 1 '\<synonym\>' | cut -f3 | tr -d [=\'=],[=\(=],[=\)=],[=\[=],[=\]=] | tr -s [=\/=] "." | tr -s "-" " " | tr -d "." | sed 's/[[:punct:]]/-/g' >> names.txt      # grep the first synonym you find as a name
                        else                                                                                                                                                                                                                                                            
                                echo >> names.txt                                                                               # write empty line if no name is available (important to keep order)                                                        
                        fi                                                                                                                                      
                    done < taxids-uniq.txt                                                                                                                                      

                paste taxids-uniq.txt names.txt > names-taxid.txt                                                               # combine tax-ids and names                                                                                                                                       and names

                while read line                                                                                                 # while reading each line of names-taxid.txt, do...                                   
                    do                                                                                                                                      
                        id=`echo $line | cut -d " " -f1`                                                                        # id is first column                                                                
                        name=`echo $line | cut -d " " -f2-`                                                                     # name is second column                                                                
                        awk -v tax="$id" '$9 == tax' contig_method.txt > tmp.txt                                                # get all lines where the 9th column equals the current taxid                                                                                        
                        num=`wc -l < tmp.txt`                                                                                   # count lines                                                    
                        if (( $num > 0 ))                                                                                       # if number of lines is greater than 0, then ...                                               
                            then                                                                                                                                        
                                printf "$name\n%.0s" $(seq 1 $num) > tmp_name.txt                                               # print the species name as often as lines were counted                                                                                        
                                paste tmp.txt tmp_name.txt > ${id}_names.txt                                                    # and paste the lines and names                                                                                   
                        fi                                                                                                                                      
                    done < names-taxid.txt                                                                                                                                      
                cat *_names.txt > contig_names_method.txt & wait                                                                # cat all names files                                                                        
                rm *_names.txt                                                                                                  # and remove them subsequently (not needed anymore)
                END_22=$(date +%s)
                DIFF_22=$(( $END_22 - $START_22 ))
                echo -ne "\nTIMING\t${DIFF_22}\tObtaining species information during assembly contigs blast"

                awk 'BEGIN {FS="\t"} {x="for_asigning.txt"} {print $2 FS $1 FS $11 FS $10 FS $9 FS $12 FS $8 > x}' < contig_names_method.txt    # change order of columns of the contig_names_method file                                                                                                                                    
                cut -f 1 454AllContigs-BlastnHits.txt > contigs.txt                                                             # cut first column to get contig names                                                                        
                #wait

                START_23=$(date +%s)
                if [[ -s Assembly/Megablast_Assigned.txt ]]                                                                     # if the file Megablast_Assigned.txt exists and is not empty, then ...
                    then                                                                                                                                        
                        while read line                                                                                         # while reading the contig-file, do...                                               
                            do                                                                                                                                      
                                    grep "$line" Assembly/454ReadInfo.txt | cut -f1,2 > ${line}.txt                             # grep $line (contig name) in 454ReadInfo and cut columns 1 and 2, to get all Reads assembled in contig                                                                                                        
                                    grep "$line" Assembly/Megablast_Assigned.txt | cut -f2-3 > read-${line}.txt                 # grep $line (contig name) in Megablast_Assigned.txt and cut columns 2-3 to get all Reads blasted to contig                                                                                                                       
                            done < contigs.txt                                                                                                                                                                      
                    else                                                                                                        # if the Megablast_Assigned.txt does not exists or is empty, then ...                                
                        while read line                                                                                                                                         
                            do                                                                                                                                      
                                grep "$line" Assembly/454ReadInfo.txt | cut -f1,2 > ${line}.txt                                 # grep $line (contig name) in 454ReadInfo and cut columns 1 and 2, to get all Reads assembled in contig                                                                                                     
                            done < contigs.txt                                                                                                                                              
                fi                                                                                                                                      

                while read line;                                                                                                # while reading the asigning file                                        
                    do                                                                                                                                      
                        t=${line}                                                                                               # save line in $t                                       
                        m=`echo $t | cut -d " " -f1` ;                                                                          # get the first column of $t (contig_name)                                                           
                        u=`wc -l < ${m}.txt` ;                                                                                  # count the number of reads belonging to contig $t                                                    
                        printf "$line\n%.0s" $(seq 1 $u) > info                                                                 # print $line as often as reads belonging to contig                                                                    
                        paste info ${m}.txt | cut -f2- > ${m}.info                                                              # paste info and reads and subsequently delete first column (contig name)                                                                        
                    done < for_asigning.txt                                                                                     

                if [[ -s Assembly/Megablast_Assigned.txt ]]                                                                     # if any Reads were assigned by Megablast to the Contigs, then ...
                    then                                                                                                      
                        while read line                                                                                         # while reading each line for "for_asigning.txt"
                            do                                                                                             
                                t=${line}                                                                                       # save line in $t
                                m=`echo $t | cut -d " " -f1` ;                                                                  # cut first column (contig name)
                                u=`wc -l < read-${m}.txt`                                                                       # count how many reads were blasted to contig
                                if (( $u > 0 ))                                                                                 # if more than 0 reads were matched to contig, then ...
                                    then                                                                        
                                        printf "$line\n%.0s" $(seq 1 $u) > readinfo                                             # print line as often as reads were assigned
                                        paste readinfo read-${m}.txt | cut -f2-8 > read-${m}.info                               # paste info and reads and cut columns 2-8 (1st is contig name, 9- are read blast infos)
                                fi                                                                                            
                            done < for_asigning.txt                                                                          
                fi                                                                                                        

                END_23=$(date +%s)
                DIFF_23=$(( $END_23 - $START_23 ))
                echo -ne "\nTIMING\t${DIFF_23}\tPost processing with species information during assembly contigs blast"

                cat *.info >> Assembly_Reads.txt                                                                                # cat all info files (info for all assigned Reads)
                cat Assembly_Reads.txt >> ${arbeitsvz}/asignedAcc.txt                                                           # add assigned Reads to asignedAcc.txt
                cut -f7 ${arbeitsvz}/asignedAcc.txt | sort | uniq > ${arbeitsvz}/excludeAcc.txt                                 # get a list of reads to exclude for further steps
                comm -3 ${arbeitsvz}/excludeAcc.txt ${arbeitsvz}/allAcc.txt | sed -e 's/^[ \t]*//' > ${arbeitsvz}/restAcc.txt   # get rest Reads
                restreads=`wc -l ${arbeitsvz}/restAcc.txt | cut -d " " -f1`                                                     # count rest Reads
        fi
        END_ASSEMBLYPROC=$(date +%s)
        DIFF_ASSEMBLYPROC=$(( $END_ASSEMBLYPROC - $START_ASSEMBLYPROC ))
        echo -ne "\nTIMING\t${DIFF_ASSEMBLYPROC}\tPost processing unassigned contigs blastn results"
        #rm *info; rm *contig*
#        if [[ -s Assembly_Reads.txt ]]
#            then
#                assfailcounter=0                                                                                                # set counter back to 0
#           else
#                let assfailcounter=assfailcounter+1
#        fi
        
        echo -ne "finished! --> $(date)\n"                                                                                                  # user info
        echo -ne "\n\n-----{Blastn of Unassigned Contigs against Nucleotide Database - Start [`date`]}-----"
fi