#    functions_megablast-contigs-reads.sh - part of RIEMS - Reliable Information Extraction from Metagenomic Sequence datasets
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


function get_tid() {
    gesamt=`wc -l < Blast-Hits.txt`                                                             # get number of blast-hits
    cut -f4 Blast-Hits.txt > tid.txt                                                            # get 4th column of Blast-Hits to get taxid
    cat tid.txt | sort -u > uniq-tids.txt                                                       # sort uniq by 2nd column with hold tax-ids
    uniq uniq-tids.txt > tmp                                                                    # can not be sorted because this will destroy order of taxids respective to reads in contigs
    mv tmp uniq-tids.txt                                                                        # rename file back
    sort -n uniq-tids.txt | uniq > tmp_tax                                                      # write temporary file with uniq taxids to print on screen
    echo -e "\rData processing... 100% completed"                                               # user info
    echo -ne "Following Tax-IDs found:\n\n"                                                     # user info                                                                                                                                                                                                                                                                                                        # user info
    cat tmp_tax                                                                                 # print uniq tax-id on screen
    let rndseqextract=rndseqextract+1                                                           # add 1 do rndseqextract
}

function contigs_to_tid() {                                                                     # functions to assign contigs to tids
    cut -f1 Blast-Hits.txt > contig_order.txt                                                   # get first column of Blast-Hits to get contig name and order
    cut -f2 tid.txt > tmp.txt                                                                   # cut second column of tid to get tid
    paste contig_order.txt tmp.txt > contig_tid.txt                                             # combine contig_name with tids
    rm tmp.txt                                                                                  # remove tmp-file
    }

function reads_to_tid() {                                                                       # functions to assign contigs to tids
    cut -f1,3 Blast-Hits.txt > reads_order.txt                                                  # cut first and 3rd column to get read names and order
    cut -f2 tid.txt > tmp.txt                                                                   # cut second column of tid to get tid
    sort tmp.txt | uniq > tid_uniq.txt                                                          # get uniq taxids
    paste reads_order.txt tmp.txt > reads_tid.txt                                               # combine reads with tids
    rm tmp.txt                                                                                  # remove tmp-file
    }
    
function assign_reads() {                                                                       # assign assembled reads to contigs and therefore to organism
    if [[ -s contig_tid.txt ]]                                                                  # if assmbled contigs resulted in Blast-Hits and contigs were assigned to tids, then...
        then
            linenum=`wc -l < contig_tid.txt`
            split -l 1 contig_tid.txt contig_part-
            for i in contig_part-*
                do
                    while read line                                                                     # while reading the contig_tid.txt, do ...
                        do
                            num=`echo $line | cut -f 1 -d " "`                                          # get the contig name (1st column)
                            tid=`echo $line | cut -f 2 -d " "`                                          # get taxid (2nd column)
                            ident=`grep ${num} Blast-Hits.txt | cut -f 3`                               # get Blast Identity (3rd column of Blast-Results)
                            get_species                                                                 # get species tid (see TaxidDetermination.sh)
                            FamSkTaxDetermination                                                       # get superkingdom- and family-tid (see TaxidDetermination.sh)
                            if [[ `grep -w "^${tid}" ${taxdir}/names.dmp | grep '\<scientific name\>'` != "" ]]          # if "scientific name" for taxid can be grepped and is not an empty string, then ...
                                then
                                    name=`grep -w "^${tid}" ${taxdir}/names.dmp | grep '\<scientific name\>' | cut -f3 | tr -d [=\'=],[=\(=],[=\)=],[=\[=],[=\]=] | tr -s [=\/=] "." | tr -s "-" " " | tr -d "." | sed 's/[[:punct:]]/-/g'`    # grep $line in names.dmp, grep "scientific name" (if possible) and cut 3rd column for name
                            elif [[ `grep -w "^${tid}" ${taxdir}/names.dmp | grep -m 1 '\<synonym\>'` != "" ]]                                                                                                # else, ...
                                then
                                    name=`grep -w "^${tid}" ${taxdir}/names.dmp | grep -m 1 '\<synonym\>' | cut -f3 | tr -d [=\'=],[=\(=],[=\)=],[=\[=],[=\]=] | tr -s [=\/=] "." | tr -s "-" " " | tr -d "." | sed 's/[[:punct:]]/-/g'`       # grep the first synonym you find as a name
                                else
                                    name=`echo`                                                         # if no name is found write and empty line (important because otherwise the order is messed up)
                            fi
                            cat ${arbeitsvz}/assembly/${num}* > ${num}.txt                              # cat all files named after the contig name from the assembly folder
                            length=`wc -l < ${num}.txt`                                                 # get number of reads per contig
                            if (( $length > 0 ))                                                        # if the number of reads is greater 0, then...
                                then
                                    method=Assembly
                                    printf "${method}\t${sktax}\t${famtax}\t${tid}\t${name}\t${ident}\n%.0s" $(seq 1 $length) > ${num}_tmp1.txt    # print length-times basic info for asignedAcc.txt
                                    paste ${num}_tmp1.txt ${num}.txt > info_${num}.txt                         # and combine basic info with read-info
                            fi
                        done < $i & 
                done
            wait
            cat info* > assembledReads.txt                                                      # combine all assigned Reads into one file
            cat assembledReads.txt >> ${arbeitsvz}/asignedAcc.txt                               # add the assigned Reads to asignedAcc.txt file
            get_rest_reads                                                                      # get all unassigned reads (see functions.sh)
            restreads=`wc -l ${arbeitsvz}/restAcc.txt | cut -d " " -f1`                         # count number of rest reads
            assigned=`wc -l < assembledReads.txt`                                               # count number of reads that could be assigned by assembly
            echo -ne "\n\n$assigned could be assigned within the assembly\n                     
            and $restreads are left\n"                                                          # user info
            #rm info* 1>/dev/null 2>&1
    fi
}

function assign_rnd_reads() {    
    method=Megablast_vs_ntdb                                                                    # set method to Megablast_vs_ntdb for assignment
    if [[ -s tid_uniq.txt ]]                                                                    # if assmbled contigs resulted in Blast-Hits and contigs were assigned to tids, then...
        then
            num=`wc -l < ${multiblastrnddir}/tid_uniq.txt`                                      # get number of uniq taxids
            if (( $num < $threadsMax ))                                                         # if the number is smaller than the maximum number of threads used, then ....
                then 
                    split -l 1 tid_uniq.txt part-                                               # split the file by line
                else                                                                            # else, ...
                    ((p=$num/$threadsMax))                                                      # calculate the number of lines per file ($p)
                    split -l $p tid_uniq.txt part-                                              # and split the file by $p
            fi
            for i in part-*                                                                     # for each part-file (containing taxid(s)), do ...
                do
                    while read line                                                             # while reading the contig_tid.txt, do ...
                        do
                            tid=$line                                                           # get taxid (2nd column)
                            get_species                                                         # get species tid
                            FamSkTaxDetermination                                               # get superkingdom- and family-tid
                            if [[ `grep -w "^${tid}" ${taxdir}/names.dmp | grep '\<scientific name\>'` != "" ]]          # if "scientific name" for taxid can be grepped and is not an empty string, then ...
                                then
                                    name=`grep -w "^${tid}" ${taxdir}/names.dmp | grep '\<scientific name\>' | cut -f3 | tr -d [=\'=],[=\(=],[=\)=],[=\[=],[=\]=] | tr -s [=\/=] "." | tr -s "-" " " | tr -d "." | sed 's/[[:punct:]]/-/g'`    # grep $line in names.dmp, grep "scientific name" (if possible) and cut 3rd column for name
                            elif [[ `grep -w "^${tid}" ${taxdir}/names.dmp | grep -m 1 '\<synonym\>'` != "" ]]                                                                                                # else, ...
                                then
                                    name=`grep -w "^${tid}" ${taxdir}/names.dmp | grep -m 1 '\<synonym\>' | cut -f3 | tr -d [=\'=],[=\(=],[=\)=],[=\[=],[=\]=] | tr -s [=\/=] "." | tr -s "-" " " | tr -d "." | sed 's/[[:punct:]]/-/g'`       # grep the first synonym you find as a name
                            else
                                name=`echo`                                                     # if no name is found write and empty line (important because otherwise the order is messed up)
                            fi
                            awk -v tax="$line" '$3 == tax' reads_tid.txt > ${line}_tmp.txt      # get all lines were the 3rd column equals the current taxid and write them to file
                            length=`wc -l < ${line}_tmp.txt`                                    # count lines of file belonging to taxid
                            if (( $length > 0 )) && [[ -n $name ]]                              # if the length is greater 0 and name is not am empty string, then ...
                                then
                                    printf "${method}\t${sktax}\t${famtax}\t${tid}\t${name}\n%.0s" $(seq 1 $length) > ${line}_tmp_name.txt  # print info to file for assignment
                                    paste ${line}_tmp_name.txt ${line}_tmp.txt > ${line}_names.txt                                          # paste info and reads
                            fi
                        done < $i & 
                done ; wait                                                                     # wait for loop to finish
            cat *_names.txt > reads_tid_names.txt & wait                                        # cat all files for assigned reads
            rm part-*; rm *_tmp_name.txt; rm *_tmp.txt                                          # remove tmp-files
            awk 'BEGIN {FS="\t"} {x="blastedReads.txt"} {print $1 FS $2 FS $3 FS $4 FS $5 FS $7 FS $6 > x}' < reads_tid_names.txt # change order of read_names_method.txt and write it to new file 
            cat blastedReads.txt >> ${arbeitsvz}/asignedAcc.txt                                 # add the assigned Reads to asignedAcc.txt file
            rm *_names.txt                                                                      # remove names-files
            get_rest_reads                                                                      # get all unassigned reads
            restreads=`wc -l ${arbeitsvz}/restAcc.txt | cut -d " " -f1`                         # count number of rest reads
            assigned=`wc -l < blastedReads.txt`                                                 # count number of assigned read during round of blastedReads
            echo -ne "\n\n$assigned of 1000 random reads could be assigned by Megablast \n
            and $restreads are left\n"                                                          # user info
            rm info*                                                                            # remove files with no use anymore
    fi
    }

function get_krit() {                                                                           # function to get krit-value
    if (( $restreads >= 100000 ))                                                               # if more than 100000 read left, then ...
        then                                                                                    
            krit=2                                                                              # krit = 2    
    elif (( $restreads >= 80000 && $restreads < 100000 ))                                       # if more than 80000 and less than 100000, then ...                                            
        then                                                                                    
            krit=3                                                                              # krit = 3    
    elif (( $restreads >= 60000 && $restreads < 80000 ))                                        # if more than 60000 and less than 80000, then ...                                            
        then                                                                                    
            krit=4                                                                              # krit = 4    
    elif (( $restreads >= 40000 && $restreads < 60000 ))                                        # if more than 40000 and less than 60000, then ...                                            
        then                                                                                    
            krit=5                                                                              # krit = 5    
        else                                                                                    # else ...
            krit=6                                                                              # krit = 6    
    fi                                                                                  
}

function get_tid_read() {
    cut -f4 Blast-Hits.txt 2>/dev/null  > tid.txt                                               # get second column of Blast-Hits and then get second column again by "|" to get gis
    reads_to_tid                                                                                # assign reads to taxid (see above)
    assign_rnd_reads                                                                            # assign random reads (see above)
    cut -f2 tid.txt | sort | uniq > uniq-tids.txt                                               # sort uniq tax-ids
    gesamt=`wc -l < ${multiblastrnddir}/uniq-tids.txt`                                                              # get number of uniq tax-ids
    paste Blast-Hits.txt tid.txt > gis.txt                                                      # combine blast-hits and tax-ids
    while read line                                                                             # while reading uniq-tids.txt, do ...
        do                                                                                      
            count=`grep -c "$line" gis.txt`                                                     # count the appearance of tax-id in gis-file
            echo "$line $count" >> Tax-Count.txt                                                # write tax-id and count in Tax-Cound.txt
        done < uniq-tids.txt 
    echo -e "\rData processing... 100% completed"                                               # user info about progress
    echo -ne "Following Tax-IDs found:  \n"                                                     # user info about found tax-ids
    cat uniq-tids.txt                                                                           # print tax-id on screen
    restreads=`wc -l < ${arbeitsvz}/restAcc.txt`                                                # count unassigned reads

    get_krit                                                                                    # get krit-value (see above)
# get only tax-ids that are great or equal to a given krit-value
    if [ -s Tax-Count.txt ]                                                                     # if Tax-Count.txt exists and is not empty, then ...
        then
            while read line                                                                     # while reading Tax-Count.txt, do ...
                do 
                    if (( `echo "$line" | cut -d " " -f2` >= $krit ))                           # if 2nd column of Tax-Count.txt line is great or equal to $krit, then ...
                        then 
                            echo "`echo "$line" | cut -d " " -f1`" >> Tax-Count-bigger10.txt    # add tax-id to Tax-Count-bigger10.txt
                    fi  
                done < Tax-Count.txt
            if [ -s Tax-Count-bigger10.txt ]                                                    # if Tax-Count-bigger10.txt exists and is not empty, then ...
                then 
                    echo -ne "following Tax-IDs will be further used:  \n"                      # user info
                    cat Tax-Count-bigger10.txt                                                  # about tax-ids
                    mv Tax-Count-bigger10.txt uniq-tids.txt                                     # and move Tax-Count-bigger10.txt to uniq-tids.txt
                    identity=97                                                                 # set identity to assign reads
                    method=Mapping                                                              # set method to assign reads
                    asign_to_tid                                                                # to assign reads to tax-id (see above)
            fi
            
    fi 
    echo
    rm MultiBlast/RndReads/gis.txt; rm MultiBlast/RndReads/tid.txt; rm MultiBlast/RndReads/uniq-tids.txt # remove files
    let rndseqextract=rndseqextract+1                                                           # add 1 to rndseqextract
}
