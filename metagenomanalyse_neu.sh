#    Metagenomanalyse.sh - part of RIEMS - Reliable Information Extraction from Metagenomic Sequence datasets
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
    
date1=`date`
################################################################# Parameterdefinitionen ##################################################################

echo -e "\n
***************** RIEMS 4.0 - Reliable Information Extraction of Metagenomic Sequence datasets ******************    
*                                                                                                               *
*                 Copyright (C) 2009-2016  Ariane Belka, Maria Jenckel, Matthias Scheuch, Dirk Hoeper           *
*                 This program comes with ABSOLUTELY NO WARRANTY                                                *
*                 This is free software, and you are welcome to redistribute it                                 *
*                 under certain conditions.                                                                     *
*                                                                                                               *
*****************************************************************************************************************\n    
RIEMS 4.0 Workflow started at `date`\n\n"                                                                           # user info

# see Riems.sh
        
####################################################################### Funktionen ######################################################################

. ${installdir}/TrimSff.sh                                                                                      # starts the trimming for more see TrimSff.sh                                              
    
trim=`wc -l < ${arbeitsvz}/asignedAcc.txt `                                                                     # count number of reads that were discarded after trimming
untrim=`wc -l < ${arbeitsvz}/allAcc.txt `                                                                       # count number of total reads
${embossdir}/seqret -sequence ${arbeitsvz}/TrimmedReads.fastq -outseq ${arbeitsvz}/TrimmedReads.fasta           # convert the fastq-file into a fasta-file    
echo -ne "\nA total of $untrim Reads was added to the metagenomic workflow \n     
$trim Reads were marked as low quality Reads \n\n"                                                              # user info

######################################################################### ViRIEMS ########################################################################

if [[ $viriems == "y" ]]
    then
        . ${installdir}/ViRiems.sh
fi

######################################################################### Assembly #######################################################################

restreads=`wc -l < ${arbeitsvz}/restAcc.txt`                                                                    # count number of remaining reads

if [ -s ${arbeitsvz}/uniq-tids.txt ]                                                                            # if taxid were given in advance
    then    
        echo -ne "Following tax-ids have been given for prior assignment: ${tax[*]} "                           # user info
        echo -ne "\nStarting with assigning the reads to given tax-ids "                                        # user info
        identity=90                                                                                             # set identity to 90 for pre-screening
        method=Pre-Screening                                                                                    # set method to pre-screening
        asign_to_tid                                                                                            # reads will be assigned to given taxid in advance to RIEMS workflow (see functions.sh)
        rm ${arbeitsvz}/uniq-tids.txt                                                                           # remove uniq-tids.txt after assigning reads
        cd ${arbeitsvz}                                                                                         # change to working directory
        identity=97                                                                                             # set identity to 97 for first Mapping
        method=Mapping                                                                                          # set method to Mapping
        subset_analyses                                                                                         # start with the analysis of subsets (functions.sh)
        echo $restreads >> ${arbeitsvz}/report.txt
        rndseqextract=0         # rndseqextract holds Info about status of assembly... set to 0  possible status for rndseqextract.txt:
                                # 0 beginning of mapping and assembly & after success # 1 no contigs in Assembly or assembly failed or no new species in Contigs  # 2 failed twice (see status 1)
    else
        method=Mapping                                                                                          # if no taxids were given, set method to mapping
        identity=97                                                                                             # set identity to 97
        subset_analyses                                                                                         # and start the analysis of subsets
        echo -ne "$restreads\t$method\n" >> ${arbeitsvz}/report.txt
 fi        
################################################################# Megablast der Rnd. Reads ###############################################################

restreads=`wc -l ${arbeitsvz}/restAcc.txt | cut -d " " -f1`                                                     # count all remaing reads

rndseqextract=0                                                                                                 # set the rndseqextract-counter to 0
cd ${arbeitsvz}                                                                                                 # change working directory
shuf -n 1000 ${arbeitsvz}/restAcc.txt > ${arbeitsvz}/tmp1.txt                                                   # get 5000 random accnos
${gsflxdir}/fnafile -i ${arbeitsvz}/tmp1.txt -o pick.fna ${arbeitsvz}/TrimmedReads.fasta         # and use them to create fna sequences for these reads
mbcr=read                                                                                                       # set mbcr to read
. ${installdir}/MegaBlast-Contigs-Reads.sh                                                                      # call MegaBlast-Contigs-Reads.sh
restreads=`wc -l < ${arbeitsvz}/restAcc.txt `                                                                   # Zählen der übrig gebliebenen Reads für die nächste Runde
echo -ne "$restreads\t$method\n" >> ${arbeitsvz}/report.txt

while (( ${restreads} > $cutoffRndReads && rndseqextract < 2 ))                                                 # as long as more than 250000 (cutoff, see Config.txt) reads are unassigned and the rndseqextract-counter is less than 2, do ...
    do
        cd ${arbeitsvz}                                                                                         # change to working directory
        shuf -n 1000 ${arbeitsvz}/restAcc.txt > ${arbeitsvz}/tmp1.txt                                           # get 5000 random accnos
        ${gsflxdir}/fnafile -i ${arbeitsvz}/tmp1.txt -o pick.fna ${arbeitsvz}/TrimmedReads.fasta                # and use them to create fna sequences for these reads
        mbcr=read                                                                                               # set mbcr to read
        . ${installdir}/MegaBlast-Contigs-Reads.sh                                                              # call MegaBlast-Contigs-Reads.sh
        restreads=`wc -l < ${arbeitsvz}/restAcc.txt `                                                           # count all remaining reads
        echo -ne "$restreads\t$method\n" >> ${arbeitsvz}/report.txt
    done

################################################################# Mapping 2 #############################################################################

if [[ -s ${arbeitsvz}/references.txt ]]                                                                         # if a references.txt file exists in the working directory (if previous analysis steps were successful), then ...
    then
        echo -ne "\n\nWe will map again with a lower cutoff (90%) against all previous detected species\n"      # user info
        method=Mapping2                                                                                         # set method to Mapping2
        identity=90                                                                                             # set identity to 90 
        while read line                                                                                         # while reading the references.txt file, do ...
            do
                i=`echo $line | sed 's/.*\///'`                                                                 # i is the reference without the whole path
                name=`echo "${line##*/}" | sed 's/.fna//'`                                                      # get name from line without the ending .fna
                organism=`basename $line`                                                                       # get organism as basename from line
                name2=`echo $name | sed 's/.*_//'`                                                              # get name2 by deleting TID-XXX_
                tid=`echo $name | sed 's/_.*//' | sed 's/.*-//' `                                               # extract TID from line
                ziel=${arbeitsvz}/${name}                                                                       # set working directory (combination of working directory and name)
                zielordner=${ziel}/Mapping2                                                                     # set new directory as Mapping2
                cd ${ziel}                                                                                      # change directory
                num=`wc -l < ${arbeitsvz}/restAcc.txt`                                                          # count number of remaining reads
                echo -ne "\n\nMapping again against $name2\n"                                                   # user info
                if (( $num > 20000 ))                                                                            # if more than 10000 reads are left, then ...
                    then
                        mkdir Mapping2                                                                          # create directory Mapping2
                        cd Mapping2                                                                             # and change directory to Mapping2
                        mem=`free -bt | tail -n 1 | cut -f 5 -d " "`                                                                    # get free memory
                        size=`du -b $line | cut -f 1`                                                                               # get size of reference
                        ((j=$mem/(2*$size)))                                                                                            # calculate how much mappings can be startet in parallel (twice the size because reads have to loaded and software needs memory itself)
                        if (( $j > $threads ))                                                                                           # if more processes could be started than cores available, then ...
                            then
                                j=$threads                                                                                              # set processes to start to cores available
                        fi
                        if (( size > 100000000 ))
                            then
                                ((j=$j/4))
                                if (( j == 0 ))
                                    then
                                        j=1
                                fi
                        fi
                        ((p=$num/$j))                                                                                                   # calculate p by number of remaining reads and core to be used
                        split -a 5 -d -l $p ${arbeitsvz}/restAcc.txt part-                                                                      # split file of remaining accno. by p lines
                        for j in ${ziel}/Mapping2/part-*                                                                         # for each accession subset, do...
                            do
                                ${gsflxdir}/runMapping -o ${ziel}/Mapping2/map-`basename ${j}` -no -noace -nobam -noinfo -notrim -force -m -n -ml 95% -np -mi $identity -ud -cpu 24 -tr -fi ${j} ${line} ${arbeitsvz}/TrimmedReads.fastq  1>/dev/null & 
                            done                                                                                # runMapping including only the subset against the organism ($line)
                            wait 1>/dev/null 2>&1                                                               # and wait until all Mappings have finished
                        for k in ${ziel}/Mapping2/map-part-*                                                                     # for each mapping directory, do...
                            do  
                                cat ${k}/454ReadStatus.txt >> ${ziel}/Mapping2/454ReadStatus.txt                      # get the ReadStatus of each subdirectory and cat it to on file in Mapping2
                            done
                    else                                                                                        # if you have less than 10000 reads left, just run the mapping with all reads
                        ${gsflxdir}/runMapping -force -noace -nobam -noinfo -no -m -n -np -mi $identity -ml 95% -ud -cpu $threads -notrim -tr -o ${ziel}/Mapping2 -fi ${arbeitsvz}/restAcc.txt ${line} ${arbeitsvz}/TrimmedReads.fastq 1>/dev/null
                        if [[ -s ${ziel}/Mapping2/454ReadStatus.txt ]]                                                  # if the mapping was successfull and a 454ReadStatus.txt exists, then ...
                            then
                            cp ${ziel}/Mapping2/454ReadStatus.txt .                                                     # copy it to Mapping2
                        fi
                fi
                if [[ -s ${ziel}/Mapping2/454ReadStatus.txt ]]                                                                   # if you have a 454ReadStatus.txt in Mapping2 (Mapping was successful), then ...
                    then
                    get_mapping_results                                                                         # get the maaping results and assign them (see functions.sh)
                fi
            done < ${arbeitsvz}/references.txt
fi           

echo -ne "$restreads\t$method\n" >> ${arbeitsvz}/report.txt                                                     # write restreads and mehtod to the report.txt
################################################################# Rest Assembly ##########################################################################

. ${installdir}/Assembly-ContigBlast.sh                                                                         # start Assembly-ContigBlast.sh

assfailcounter=0                                                                                                # set assfailcounter to 0 (possible values 0, 1 ,2 for more info see rndseqextract counter)
until (( ${assfailcounter} > 4 ))                                                      # until you cannot assign any more reads or the assfailcounter is 2, do ...
    do
        rm -r ${arbeitsvz}/Assembly/                                                                            # remove Assembly directory
        let assfailcounter=assfailcounter+1                                                                     # count assfailcounter plus 1        
        . ${installdir}/Assembly-ContigBlast.sh                                                                 # start Assembly-ContigBlast.sh
    done
echo -ne "$restreads\t$method\n" >> ${arbeitsvz}/report.txt

############################################################### Identify all Orgs and download sequences ###################################################
echo -ne "Organisms identified within the final Assembly will be identified and reference sequences will be downloaded for Blastn vs all organisms.\n"                                          # user info
cut -f4 ${arbeitsvz}/asignedAcc.txt | sort | uniq > ${arbeitsvz}/uniq_tid.txt                                   # get all uniq taxids from asignedAcc.txt (all previous identified species)
while read line                                                                                                 # while reading the uniq_tid.txt, do...
    do
        tid=$line                                                                                               # line equals taxid
        parentTax=$line                                                                                         # as well as parentTax
        if [[ `grep -w "$parentTax" ${arbeitsvz}/AnBlast.txt` == "" ]]                                          # if the found taxid is not already present in AnBlast.txt (all previous taxids are saved in the file), then ...
            then
                . ${installdir}/TaxID-Sequenzisolation.sh                                                       # start SequenceIsolation.sh to get reference sequences
                echo ${refseqdir}/TID-${parentTax}_$name.fna >> ${arbeitsvz}/references.txt                     # add the reference to references.txt
                grep "^>" ${refseqdir}/TID-${parentTax}_$name.fna | cut -f2 -d "|" >> ${arbeitsvz}/Orgs_gis.txt
                ln -s ${refseqdir}/TID-${parentTax}_$name.fna ${arbeitsvz}/DB/`basename ${refseqdir}/TID-${parentTax}_$name.fna`       #and set a link to the DB-directory     
        fi
    done < ${arbeitsvz}/uniq_tid.txt

################################################################## Rest Megablast #########################################################################

method=Megablast_vs_ntdb                                                                                        # set method (needed for output)

    ${gsflxdir}/fnafile -o ${arbeitsvz}/Ausgangsdateien/Megablast.fna -i ${arbeitsvz}/restAcc.txt ${arbeitsvz}/TrimmedReads.fasta    # write all unassigned reads into fasta-file
    blastvek=(0 0 tax 0 0)                                                                                      # set blastvector for blast.sh              
    query=${arbeitsvz}/Ausgangsdateien/Megablast.fna                                                            # set query for blast   
    referenz=${ntdb}                                                                                            # set reference for blast 
    blastordner=${arbeitsvz}/MultiBlast/Megablast                                                               # define blast-folder  
    ntblast=megablast                                                                                           # set further parameters for blastn (e.g. set blast to megablast)
    . ${installdir}/Blast.sh                                                                                    # call blast.sh
    getOutput                                                                                                   # get Output from blast and assign it to reads
echo -ne "$restreads\t$method\n" >> ${arbeitsvz}/report.txt
################################################################ Blastn vs Orgs #######################################################################
method=Blastn_vs_Organism                                                                                       # set method (needed for output)

    ${gsflxdir}/fnafile -o ${arbeitsvz}/Ausgangsdateien/Blastn-vs-orgs.fna -i ${arbeitsvz}/restAcc.txt ${arbeitsvz}/TrimmedReads.fasta    # write all unassigned reads into fasta-file
    ${blastdir}/blastdb_aliastool -db ${ntdb} -gilist ${arbeitsvz}/Orgs_gis.txt -dbtype nucl -out ${arbeitsvz}/DB/AllOrgs.DB -title "AllOrgs"
    #cat ${arbeitsvz}/DB/TID-*.fna > ${arbeitsvz}/DB/AllOrgs.fna                                                # cat all previous identified organism-fasta fils into AllOrgs.fna
    blastvek=(0 db tax 0 0)                                                                                     # set blastvector for blast.sh
    query=${arbeitsvz}/Ausgangsdateien/Blastn-vs-orgs.fna                                                       # set query for blast
    referenz=${arbeitsvz}/DB/AllOrgs.DB                                                                         # set reference for blast
    blastordner=${arbeitsvz}/MultiBlast/Blastn-vs-Orgs                                                          # define blast-folder 
    ntblast=blastn                                                                                              # set further parameters for blastn (e.g. set blast to megablast)
    . ${installdir}/Blast.sh                                                                                    # call blast.sh
    getOutput                                                                                                   # get Output from blast and assign it to reads
echo -ne "$restreads\t$method\n" >> ${arbeitsvz}/report.txt    
#################################################################### Rest Blastn #########################################################################

method=Blastn_vs_ntdb                                                                                           # set method (needed for output)

    ${gsflxdir}/fnafile -o ${arbeitsvz}/Ausgangsdateien/Blastn.fna -i ${arbeitsvz}/restAcc.txt ${arbeitsvz}/TrimmedReads.fasta   # write all unassigned reads into fasta-file
    blastvek=(0 0 tax 0 0)                                                                                      # set blastvector for blast.sh              
    query=${arbeitsvz}/Ausgangsdateien/Blastn.fna                                                               # set query for blast   
    referenz=${ntdb}                                                                                            # set reference for blast 
    blastordner=${arbeitsvz}/MultiBlast/Blastn                                                                  # define blast-folder  
    ntblast=blastn                                                                                              # set further parameters for blastn (e.g. set blast to blastn)
    . ${installdir}/Blast.sh                                                                                    # call blast.sh
    getOutput                                                                                                   # get Output from blast and assign it to reads
echo -ne "$restreads\t$method\n" >> ${arbeitsvz}/report.txt

#################################################################### finish Analysis #########################################################################
    
restreads=`wc -l < ${arbeitsvz}/restAcc.txt`                                                                    # count lines of restAcc.txt to get number of unclassified reads

echo -ne "\n $restreads reads couldn't be assigned\n"

if (( $restreads > 0 ))                                                                                         # if number of restReads is greater 0, then...
    then 
        printf "\t\t\t\t\t\n%.0s" $(seq 1 $restreads) > zwischendurch1.txt                                      # print tabs as often as unclassified reads were counted 
        printf "unclassified\n%.0s" $(seq 1 $restreads) > zwischendurch2.txt                                    # print "unclassified" as often as unclassified reads were counted
        paste zwischendurch1.txt ${arbeitsvz}/restAcc.txt zwischendurch2.txt > zwischendurch3.txt               # and combine both files with unclassified read accessions
        mv zwischendurch3.txt asigned.txt                                                                       # and move back to assigned
        rm zwischendurch*.txt                                                                                   # remove temporary file
        cat asigned.txt >> ${arbeitsvz}/asignedAcc.txt                                                          # and add unclassified reads to results in asignedAcc
fi

cut -f 3 ${arbeitsvz}/asignedAcc.txt > ${arbeitsvz}/famtax.txt                                                  # cut the 3rd column of asignedAcc.txt to get all famtax that were detected (necessary for resultprotocol) 
sort -n ${arbeitsvz}/famtax.txt | uniq | sed '/^$/d' > ${arbeitsvz}/uniq_famtax.txt                             # sort, uniq and delete all empty lines in famtax.txt
while read line                                                                                                 # while reading uniq_famtax.txt
    do
        famtax=$line                                                                                            # assign line to famtax
            if [[ $famtax == NA ]]                                                                              # if famtax is NA (no family found), then ...
                then 
                    echo "NA" >> ${arbeitsvz}/names.txt                                                         # print "NA" to names.txt, else ...
                else
                    grep "^\<$famtax\>" ${taxdir}/names.dmp | grep '\<scientific name\>' | cut -f3 >> ${arbeitsvz}/names.txt    # grep famtax in names.dmp and grep "scientific name" (3rd column) and write it to names.txt
            fi
    done < ${arbeitsvz}/uniq_famtax.txt
paste ${arbeitsvz}/uniq_famtax.txt ${arbeitsvz}/names.txt > ${arbeitsvz}/famtax_names.txt                       # paste famtax and names to get a table for the resultprotocol
sort ${arbeitsvz}/asignedAcc.txt | uniq > ${arbeitsvz}/tmp.txt                                                  # sort and uniq all asigned Acc to avoid potential duplicates
mv ${arbeitsvz}/tmp.txt ${arbeitsvz}/asignedAcc.txt                                                             # and rename tmp.txt file
date2=`date`                                                                                                    # get time when analysis finished
echo -ne "\n --- BASIC METAGENOMIC ANALYSIS FINISHED --- `date`\n"                                              # user info

echo -ne "\nA result report will be created...\n"                                                               # user info

emboss=`${embossdir}/embossversion 2>/dev/null`                                                                              # get emboss version
newbler=`${gsflxdir}/newbler -version | head -n1`                                                                          # get newbler version (1st row)
newblertail=`${gsflxdir}/newbler -version | head -n2 | tail -n1`                                                            # get newbler version (2nd row)
catnewbler="$newbler $newblertail"                                                                              # cat both lines of newbler version to one line
blastv=`${blastdir}/blastn -version | tail -n1`                                                                   # get blastversion
blastdbv=`${blastdir}/blastdbcmd -db /home/blast_db/ncbi/nt -info | grep "Date" | cut -f1 | sed 's/Date: //g'`              # get version of blast database
texversion=`${latexdir}/tex --version | head -n 1`

export arbeitsvz                                                                                                # export variables to use in R to create resultprotocol
#export input                                                                                                   # |
export emboss                                                                                                   # |
export catnewbler                                                                                               # | 
export date1                                                                                                    # |
export date2                                                                                                    # |
export blastv                                                                                                   # |
export blastdbv                                                                                                 # V
export texversion
cd $arbeitsvz                                                                                                   # change bask to working directory

############################################# Create resultprotocol ################################################################################
echo "Construct: *.txt" >> zeiten.log                                                                           # save info
date >> zeiten.log                                                                                              # save date
${R}/Rscript --vanilla ${scriptR} 12 &>>${arbeitsvz}/console.log                                                # run R script to get result tables (see resultprotocol.R)
wait                                                                                                            # wait until R anaylsis has finished
date >> zeiten.log                                                                                              # save date

. ${installdir}/knitr_result.sh                                                                                 # run knitr_result.sh to get result report PDF

echo ${arbeitsvz}/TID*/Mapping2/map-part* | xargs rm -r                                                                       # remove files that will no longer be used
echo ${arbeitsvz}/TID*/MappingFiles/map-part* | xargs rm -r                                                                   # |                          
rm -r ${arbeitsvz}/DB/                                                                                          # |   
rm -r ${arbeitsvz}/D*Mapping/                                                                                   # |           
#rm -r ${arbeitsvz}/MultiBlast/Blastn-vs-Orgs/DB                                                                 # V                           
echo -ne "Done\n"                                                                                               # user info
#furtherAnalysis=y                                                                                           
if [[ $furtherAnalysis == y ]]                                                                                  # if a further analysis should be performed, then ...
    then
        echo -ne "--- START OF FURTHER ANALYSIS ---\n"                                                          # user info
        echo -ne "Subsequently, the result report will be updated\n"                                            # user info

        . ${installdir}/furtherAnalysis.sh                                                                  # run further anaylsis 
fi
