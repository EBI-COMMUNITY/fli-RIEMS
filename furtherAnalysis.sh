#    furtherAnalysis.sh - part of RIEMS - Reliable Information Extraction from Metagenomic Sequence datasets
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

. ${installdir}/Config.txt 
. ${installdir}/functions.sh  
. ${installdir}/Blast_functions.sh

cd ${arbeitsvz}
dateFurther1=`date`
scriptBlastp=${installdir}/Blastp.R
scriptFurther=${installdir}/resultprotocol_further.R
scriptTFurther=${installdir}/resultprotocol_further.Rnw

################################### Assembly of RestReads + Blastp ######################################
echo -ne "$(date) --> $restreads will be assembled\n"
SECONDS=0
${gsflxdir}/runAssembly -o assemblyFA -fi ${arbeitsvz}/restAcc.txt -notrim -tr -force -noinfo -cpu ${threads} -acedir -ml 20 ${arbeitsvz}/TrimmedReads.fastq 1>/dev/null
ContigRead=Contig
method=Blastp_vs_protdb
if [[ -s ${arbeitsvz}/assemblyFA/454AllContigs.fna ]]
    then
        cd ${arbeitsvz}/assemblyFA/
        ${embossdir}/getorf -sequence 454AllContigs.fna -find 0 -minsize 90 -outseq 454AllContigs.orf -sformat1 pearson
        numorf=`grep -c "^>" 454AllContigs.orf`
        if [[ -z numorf ]]
            then
                numorf=0
        fi
        echo -ne "\nBLAST (blastp) of $numorf potential ORFs protein database...\n"
        #${blastdir}/blastp -db ${protdb} -num_threads ${threads} -negative_gilist ${installdir}/gi_exclude.txt -evalue ${evalue} -max_target_seqs 1 -out Hits.txt -query 454AllContigs.orf -max_hsps 1 -outfmt '6 qseqid sseqid qstart qend qlen evalue pident staxid length qcovs' 
        ${blastdir}/blastp -db ${protdb} -num_threads ${threads} -evalue ${evalue} -max_target_seqs 1 -out Hits.txt -query 454AllContigs.orf -max_hsps 1 -outfmt '6 qseqid sseqid qstart qend qlen evalue pident staxid length qcovs' 
        if [[ -s Hits.txt ]]
            then
                echo -e "Accno \tRefsequence-ID \tHit-Range \t\tLength \tE-value \tperc.identity \tTax-ID \t$rank Tax-ID \tSuperkingdom Tax-ID \n" > Blast-Hits.txt   # create a file with header for blast-results
                export arbeitsvz
                export ContigRead
                ${R}/Rscript --vanilla ${scriptBlastp} 12 &>>${arbeitsvz}/console.log
                numHits=`wc -l < Hits.txt`
                #dos2unix Hits.txt 1>/dev/null 2&>1
                cut -f1-5 454ReadStatus.txt > ReadStatus.txt
                while read line
                    do
                        contig=`echo $line | cut -f1 -d " "`
                        grep -w "$contig" ReadStatus.txt | cut -f 1 > ${contig}_reads.txt
                        r=`wc -l < ${contig}_reads.txt`
                        printf "$line\n%.0s" $(seq 1 $r) | cut -f2- > info.txt
                        paste ${contig}_reads.txt info.txt > ${contig}_assigned.txt
                    done < Hits.txt
            
                cat *_assigned.txt > Hits.txt
                cat Hits.txt | sort -n -t$'\t' -k 8 | grep -vw "NoHit" | cut -f1-7 > AllHits.txt                                    # get Hits and sort them by taxid (8th column, tab separated) and remove of NoHits (probably not needed) and 
                numHits=`wc -l < Hits.txt`
                echo -ne "$numHits protein sequence(s) could be identified by Blastp\n Data processing ...\n"
        
                BlastHitTaxid
                getOutput
        fi
fi
duration=$SECONDS
echo -ne "\nTotal time taken for assembly and Blastp of unassigned contigs was $(($duration / 60)) minutes and $(($duration % 60)) seconds."
echo -ne "\nTIMING $duration Unassigned contigs assembly and blastp" 

################################### Blastp of RestReads ##############################################
ContigRead=Read        
echo -ne "\n$(date) --> protein sequences of $restreads sequences will be analysed\n"
SECONDS=0

method=Blastp_vs_protdb                                                                                         # set method (needed for output)

    ${gsflxdir}/fnafile -o ${arbeitsvz}/Ausgangsdateien/Blastp.fna -i ${arbeitsvz}/restAcc.txt ${arbeitsvz}/TrimmedReads.fasta   # write all unassigned reads into fasta-file
    ${embossdir}/getorf -sequence ${arbeitsvz}/Ausgangsdateien/Blastp.fna -find 0 -minsize 90 -outseq ${arbeitsvz}/Ausgangsdateien/Blastp.orf -sformat1 pearson 
    blastvek=(0 0 tax 0 0)                                                                                      # set blastvector for blast.sh              
    query=${arbeitsvz}/Ausgangsdateien/Blastp.orf                                                               # set query for blast   
    referenz=${protdb}                                                                                          # set reference for blast 
    blastordner=${arbeitsvz}/MultiBlast/Blastp                                                                  # define blast-folder  
    ntblast=blastp                                                                                              # set further parameters for blastn (e.g. set blast to blastn)
    . ${installdir}/Blast.sh                                                                                    # call blast.sh
    getOutput                                                                                                   # get Output from blast and assign it to reads
echo -ne "$restreads\t$method\n" >> ${arbeitsvz}/report.txt

cd ${arbeitsvz}
duration=$SECONDS
echo -ne "\nTotal time taken for Blastp of rest unassigned reads was $(($duration / 60)) minutes and $(($duration % 60)) seconds."
echo -ne "\nTIMING $duration Blastp unassigned reads" 
                                                         
################################### Blastx of RestReads ##############################################

ContigRead=Read        
echo -ne "\n$(date) --> protein sequences of $restreads sequences will be analysed\n"
SECONDS=0

method=Blastx_vs_protdb                                                                                         # set method (needed for output)

    ${gsflxdir}/fnafile -o ${arbeitsvz}/Ausgangsdateien/Blastx.fna -i ${arbeitsvz}/restAcc.txt ${arbeitsvz}/TrimmedReads.fasta   # write all unassigned reads into fasta-file
    blastvek=(0 0 tax 0 0)                                                                                      # set blastvector for blast.sh              
    query=${arbeitsvz}/Ausgangsdateien/Blastx.fna                                                               # set query for blast   
    referenz=${protdb}                                                                                          # set reference for blast 
    blastordner=${arbeitsvz}/MultiBlast/Blastx                                                                  # define blast-folder  
    ntblast=blastx                                                                                              # set further parameters for blastn (e.g. set blast to blastn)
    . ${installdir}/Blast.sh                                                                                    # call blast.sh
    getOutput                                                                                                   # get Output from blast and assign it to reads
echo -ne "$restreads\t$method\n" >> ${arbeitsvz}/report.txt

cd ${arbeitsvz}
duration=$SECONDS
echo -ne "\nTotal time taken for Blastx of remaining reads was $(($duration / 60)) minutes and $(($duration % 60)) seconds."
echo -ne "\nTIMING $duration Blastx unassigned reads" 

 ############################### tBlastx of RestReads ######################################################

ContigRead=Read        
echo -ne "\n$(date) --> protein sequences of $restreads sequences will be analysed\n"
SECONDS=0

method=tBlastx_vs_ntdb                                                                                          # set method (needed for output)

    ${gsflxdir}/fnafile -o ${arbeitsvz}/Ausgangsdateien/tBlastx.fna -i ${arbeitsvz}/restAcc.txt ${arbeitsvz}/TrimmedReads.fasta   # write all unassigned reads into fasta-file
    blastvek=(0 0 tax 0 0)                                                                                      # set blastvector for blast.sh              
    query=${arbeitsvz}/Ausgangsdateien/tBlastx.fna                                                              # set query for blast   
    referenz=${ntdb}                                                                                            # set reference for blast 
    blastordner=${arbeitsvz}/MultiBlast/tBlastx                                                                 # define blast-folder  
    ntblast=tblastx                                                                                             # set further parameters for blastn (e.g. set blast to blastn)
    . ${installdir}/Blast.sh                                                                                    # call blast.sh
    getOutput                                                                                                   # get Output from blast and assign it to reads
echo -ne "$restreads\t$method\n" >> ${arbeitsvz}/report.txt

cd ${arbeitsvz}
duration=$SECONDS
echo -ne "\nTotal time taken for tBlastx of unassigned reads was $(($duration / 60)) minutes and $(($duration % 60)) seconds."
echo -ne "\nTIMING $duration tBlastx unassigned reads"

################################# Get Family Names for further analysis #####################################

restreads=`wc -l < ${arbeitsvz}/restAcc.txt`                                                                    # and rename tmp.txt file
cut -f 3 ${arbeitsvz}/asignedAcc_pep.txt > ${arbeitsvz}/famtax_pep.txt                                          # cut the 3rd column of asignedAcc.txt to get all famtax that were detected (necessary for resultprotocol) 
sort -n ${arbeitsvz}/famtax_pep.txt | uniq | sed '/^$/d' > ${arbeitsvz}/uniq_famtax_pep.txt                     # sort, uniq and delete all empty lines in famtax.txt
while read line                                                                                                 # while reading uniq_famtax.txt
    do
        famtax=$line                                                                                            # assign line to famtax
            if [[ $famtax == NA ]]                                                                              # if famtax is NA (no family found), then ...
                then 
                    echo "NA" >> ${arbeitsvz}/names_pep.txt                                                     # print "NA" to names.txt, else ...
                else
                    grep "^\<$famtax\>" ${taxdir}/names.dmp | grep '\<scientific name\>' | cut -f3 >> ${arbeitsvz}/names_pep.txt    # grep famtax in names.dmp and grep "scientific name" (3rd column) and write it to names.txt
            fi
    done < ${arbeitsvz}/uniq_famtax_pep.txt
paste ${arbeitsvz}/uniq_famtax_pep.txt ${arbeitsvz}/names_pep.txt > ${arbeitsvz}/famtax_names_pep.txt           # paste famtax and names to get a table for the resultprotocol
sort ${arbeitsvz}/asignedAcc_pep.txt | uniq > ${arbeitsvz}/tmp.txt                                              # sort and uniq all asigned Acc to avoid potential duplicates
mv ${arbeitsvz}/tmp.txt ${arbeitsvz}/asignedAcc_pep.txt

################################### Prepare Result report ##############################################

dateFurther2=`date`                                                                                             # get time when analysis finished

export dateFurther1
export dateFurther2

echo -ne "\n --- Further METAGENOMIC ANALYSIS FINISHED --- `date`\n"                                            # user info
echo -ne "\nAn additional result report will be created...\n" 

echo "Construct: *.txt" >> zeiten.log
date >> zeiten.log
${R}/Rscript --vanilla ${scriptFurther} 12 &>>${arbeitsvz}/console.log
wait
date >> zeiten.log

md5completepep=`md5sum ${arbeitsvz}/resultprotocol-complete-pep.txt | cut -d " " -f1`
md5compactpep=`md5sum ${arbeitsvz}/resultprotocol-compact-pep.txt | cut -d " " -f1`

if [[ $resultprotocol == "y" ]]
    then
        
        restreads=`cut -f8 ${arbeitsvz}/asignedAcc.txt | grep -c "unclassified"`
        lowqual=`cut -f8 ${arbeitsvz}/asignedAcc.txt | grep -c "low quality"`
        totalreads=`wc -l < ${arbeitsvz}/asignedAcc.txt`
        
        export md5completepep
        export md5compactpep
        #export restreads
        #export lowqual
        #export totalreads
        #export projektname
        
        echo "Construct: *.tex" >> zeiten.log
        date >> zeiten.log
        ${R}/Rscript -e "library(knitr); knit('${scriptTFurther}')" 12 &>>${arbeitsvz}/console.log
        date >> zeiten.log
        
        tex2pdf=${arbeitsvz}/resultprotocol_further.tex
        echo "Construct: *.pdf" >> zeiten.log
        date >> zeiten.log
        pdflatex ${tex2pdf} 12 &>>${arbeitsvz}/console.log
        pdflatex ${tex2pdf} 12 &>>${arbeitsvz}/console.log
        date >> zeiten.log
        
        mv resultprotocol_further.pdf ${projektname}-resultprotocol_further.pdf
        mv resultprotocol_further.tex ${projektname}-resultprotocol_further.tex
        mv resultprotocol_further.aux ${projektname}-resultprotocol_further.aux
        mv resultprotocol_further.log ${projektname}-resultprotocol_further.log
        mv resultprotocol-compact-pep.txt ${projektname}-resultprotocol-compact_pep.txt
        mv resultprotocol-complete-pep.txt ${projektname}-resultprotocol-complete_pep.txt
        mv resultprotocol-compact.txt ${projektname}-resultprotocol-compact.txt
        mv resultprotocol-complete.txt ${projektname}-resultprotocol-complete.txt
        echo -ne "Done --> $(date)"
            
fi
