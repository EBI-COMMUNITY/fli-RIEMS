#    functions.sh - part of RIEMS - Reliable Information Extraction from Metagenomic Sequence datasets
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

. ${installdir}/Config.txt                                                                                                  # <- Path of config file (unchanged)           
    
#aus Megablast-Contigs_Reads.sh 

function get_rest_reads() {

if [[ $method == Blastp_vs_protdb ]]
    then
        cut -f7 ${arbeitsvz}/asignedAcc_pep.txt | sort | uniq > ${arbeitsvz}/excludeAcc.txt                                 # get all assigned Accnos (column 7), sort them uniq and write to exclude list
        diff --speed-large-files ${arbeitsvz}/excludeAcc.txt ${arbeitsvz}/restAcc.txt | grep "^>" | sed 's/^> //g' > ${arbeitsvz}/tmp           # compare exclude list and allAcc-list, take all that are not included in exclude list (-3), remove tab-stops at beginning of list and write to restAcc.txt
        mv ${arbeitsvz}/tmp ${arbeitsvz}/restAcc.txt
    else
        cut -f7 ${arbeitsvz}/asignedAcc.txt | sort | uniq > ${arbeitsvz}/excludeAcc.txt                                     # get all assigned Accnos (column 7), sort them uniq and write to exclude list
        diff ${arbeitsvz}/excludeAcc.txt ${arbeitsvz}/restAcc.txt | grep "^>" | sed 's/^> //g' > ${arbeitsvz}/tmp           # compare exclude list and allAcc-list, take all that are not included in exclude list (-3), remove tab-stops at beginning of list and write to restAcc.txt
        mv ${arbeitsvz}/tmp ${arbeitsvz}/restAcc.txt
fi
}
    
function asign_to_tid() {   
    while read line                                                                                                         # while reading uniq-tids.txt, do...
        do              
            tid=`echo "$line" | cut -d " " -f1`                                                                             # cut 1st colume (tax-id) and assign it to $tid
            get_species                                                                                                     # get taxid of species
            if [[ $tid == 77133 ]] || [[ $tid == 32630 ]]                                                                   # is species is uncultured bacterium (77133) or synthetic contruct (32630), then ...
                then                
                    echo $tid >> ${arbeitsvz}/skiped.txt                                                                    # only write taxid to skip-file
                else                                                                                                        # else, ...
                    grep -q "TID-${tid}" ${arbeitsvz}/AnBlast.txt                                                           # check if reference was already detected (-q quiet)
                    if [[ `echo $?` != 0 ]]                                                                                 # if grep was unsuccessful, then ... ($? holds information about last command and if it was successful. 0=successful, 1=failed)
                        then                
                            . ${installdir}/TaxID-Sequenzisolation.sh                                                       # start TaxID-Sequenzisolation.sh  
                            echo "TID-${tid}" >> ${arbeitsvz}/AnBlast.txt                                                   # and write tax-id to Anblast.txt
                            . ${installdir}/Mapping.sh                                                                      # and start Mapping.sh
                            rndseqextract=0                                                                                 # set rndseqextract back to 0
                    fi  
            fi  
        done < uniq-tids.txt
}   

function map_vs_contigs()
{
num=`wc -l < ${mappingdir}/restAcc.txt`                                                                                     # count remaining reads
if (( num > 50000 ))                                                                                                        # if more than 50000 reads remaining do...
    then                
        #cd ${zielordner}                                                                                                   # change directory 
        ((p=$num/$threadsSplit))                                                                                            # calculate p by number of remaining reads and core to be used
        split -l $p ${mappingdir}/restAcc.txt ${mappingdir}/part-                                                           # split file of remaining accno. by p lines
        for i in ${mappingdir}/part-*                                                                                       # for each file do...
            do
                ${gsflxdir}/runMapping -o ${mappingdir}/map_vs_contigs_`basename ${i}` -force -noace -nobam -noinfo -rst 0 -m -n -np -mi 97 -ml 97% -ud -cpu 24 -tr -notrim -fi $i ${mappingdir}/454AllContigs.fna ${arbeitsvz}/TrimmedReads.fastq 1>/dev/null & 
            done ; wait
        cat ${mappingdir}/map_vs_contigs_*/454ReadStatus.txt > ${mappingdir}/MapResults.txt 
    else
        ${gsflxdir}/runMapping -o ${mappingdir}/map_vs_contigs -force -noace -nobam -noinfo -rst 0 -m -n -np -mi 97 -ml 97% -ud -cpu 24 -tr -notrim -fi ${mappingdir}/restAcc.txt ${mappingdir}/454AllContigs.fna ${arbeitsvz}/TrimmedReads.fastq 1>/dev/null
        cp ${mappingdir}/map_vs_contigs/454ReadStatus.txt ${mappingdir}/MapResults.txt 
fi    
rm -r ${mappingdir}/map_vs_contigs*  
if [ -s ${mappingdir}/MapResults.txt ]                                                                                      # if 454AllContigs.fna existes and is not empty, then ...
    then                                                                                            
        while read line                                                                                                     # while going through all the contigs in contigs.txt, do ...
            do
                grep "${line}" ${mappingdir}/MapResults.txt | grep -w "Full" | cut -f3-6 > ${mappingdir}/${line}_map.txt    # grep all the reads that were assembled into the contig
            done < ${mappingdir}/contigs.txt
fi
}


function subset_analyses()
{
while (( ${restreads} > $cutoffAssembly && rndseqextract < 2 )) || [[ ! -s ${arbeitsvz}/references.txt ]]                   # as long as more than 10000 reads are remaining and assembly has not failed twice in a row continue
    do                                                                                                                      
        identity=97         
        echo -ne "\nA subset of 50000 reads will be created... "            
        shuf -n 50000 ${arbeitsvz}/restAcc.txt > ${arbeitsvz}/tmp1.txt                                                      # get a list of 5000 random accnos
        echo "assembly of subset... "
        ${gsflxdir}/runAssembly -o ${arbeitsvz}/assembly -fi ${arbeitsvz}/tmp1.txt -force -cpu ${threads} -mi 97 -notrim -nobig -noinfo -a 100  ${arbeitsvz}/TrimmedReads.fastq 1>/dev/null ; wait   # and assemble them
        # Assembly of subset; -a 100 (min length of contig) -force overwrite if existing; -nobig, -noinfo suppress big data files, -mi minimum identity of reads to be assembled
        echo "finished"
# get contigs by numreads or depth ...
        mappingdir=${arbeitsvz}/assembly
        cd ${mappingdir}                                                                                                    # change directory
    if [ -s ${mappingdir}/454AllContigs.fna ]                                                                               # if assembly was successful and contigs are available, then ...
        then            
            cut -f1-3 ${mappingdir}/454ReadStatus.txt > ${mappingdir}/454ReadInfo.txt           
            grep "contig" ${mappingdir}/454ContigGraph.txt | cut -f2- > ${mappingdir}/contigGraph.txt                       # get all contigs from contigGraph.txt and delete first column
            csplit -s ${mappingdir}/454AllContigs.fna /\>/ {*} -f "contig-" -n 5 -z -q                                      # split the 454AllContigs-file at the header so each contig is in a seperate file
            mv ${mappingdir}/454AllContigs.fna ${mappingdir}/454AllContigs-original.fna                                     # move 454AllContigs.fna to 454AllContigs-original.fna
            for i in ${mappingdir}/contig-*                                                                                 # for each file containing a contig, do ...
                do                                                                                                          
                numread=`grep 'numreads*' $i | cut -d "=" -f3`                                                              # get the number of reads used to build the contig
                contig=`grep "^>" ${i} | cut -d " " -f1 | tr -d ">"`                                                        # get the contig name (e.g contig00001)
                depth=`grep "${contig}" contigGraph.txt | cut -f3`                                                          # get depth for this contig from contigGraph
                round=`printf "%.0f" $depth`                                                                                # round depth value
                if (( ${numread} >= 10 || ${round} >= 4 ))                                                                  # if more than 10 reads were used or depth is higher than 4, then ...
                    then                                                                                                        
                        cat ${i} >> ${mappingdir}/454AllContigs.fna                                                         # add the contig to 454AllContig.fna 
                fi          
                done            
            grep ">" ${mappingdir}/454AllContigs.fna > ${mappingdir}/tmp.txt                                                # get all header and write them in temporary file
            echo -ne `wc -l < tmp.txt` contigs were generated by assembly         
            sort -t "=" -k 3 -nr ${mappingdir}/tmp.txt > ${mappingdir}/ContigInfo.txt                                       # sort tmp-file by number of reads in contig (-t "=" set delimiter to "=", -k 3 column 3, -nr sort numeric in decreasing oder)
            cut -f1 -d  " " ${mappingdir}/tmp.txt | tr -d ">" > ${mappingdir}/contigs.txt           
                        
            if [ -s ${mappingdir}/454AllContigs.fna ]                                                                       # if 454AllContigs.fna existes and is not empty, then ...
                then                                                                                                        
                    while read line                                                                                         # while going through all the contigs in contigs.txt, do ...
                        do
                            grep "${line}" ${mappingdir}/454ReadInfo.txt | grep -w "Assembled" | cut -f1-2 > ${mappingdir}/${line}.txt                  # grep all the reads that were assembled into the contig
                        done < contigs.txt
                    cat ${mappingdir}/contig0*.txt | sort > ${mappingdir}/assembledReads.txt
                    cut -f 1 ${mappingdir}/assembledReads.txt | sort | uniq > ${mappingdir}/assembledAcc.txt
                    diff ${mappingdir}/assembledAcc.txt ${arbeitsvz}/restAcc.txt | grep "^>" | sed 's/^> //g' > ${mappingdir}/restAcc.txt               # compare exclude list and allAcc-list, take all that are not included in exclude list (-3), remove tab-stops at beginning of list and write to restAcc.txt   # compare exclude list and allAcc-list, take all that are not included in exclude list (-3), remove tab-stops at beginning of list and write to restAcc.txt
                    map_vs_contigs          
                    echo -ne "\nSpecies will be identified for assembled contigs... "                                       # user info
                    mbcr=contig                                                                                             # set mbrc to contig                                                                                                                                                                    
                    . ${installdir}/MegaBlast-Contigs-Reads.sh                                                              # start MegaBlast-Contig-Read.sh tool    
                else                                                                                                
                    let rndseqextract=rndseqextract+1                                                                       # else 1 is added to the rndseqextract-counter,
                    #rm ${arbeitsvz}/Blastn/* 2>/dev/null                                                                   # and all filed in the Blastn-folder will be deleted
            fi                                                                                                      
        else                                                                                                    
            let rndseqextract=rndseqextract+1                                                                               # if no contigs are available add 1 to the rndseqextract-counter
            #rm ${arbeitsvz}/Blastn/* 2>/dev/null                                                                           # and all filed in the Blastn-folder will be deleted
        fi          
                    
    rm ${arbeitsvz}/assembly/contig*                                                                                        # delete all files in the assembly folder
    cd ${arbeitsvz}/Ausgangsdateien                                                                                         # change directory
    rm -r ${arbeitsvz}/assembly                                                                                             # remove the assembly directory
done            
}   

function get_mapping_results() 
{
if [[ -s ${zielordner}/454ReadStatus.txt ]]                                                                                 # if ReadStatus exisits (Mapping was successful), then ...
    then
        egrep "Full|Partial|Repeat|Chimeric" ${zielordner}/454ReadStatus.txt | cut -f3-7 > asigned.txt                      # grep all reads that mapped (somehow) and get column 3-7 (1-2 have no useful info)
        i=`wc -l < asigned.txt`                                                                                             # count number of mapped reads
        
        FamSkTaxDetermination                                                                                               # determine fam- and sk-tax (TaxidDetermination.sh)
        
        name=`echo $organism | sed 's/.fna//' | sed 's/.*_//' | sed 's/-/ /g' | tr -d "." | sed 's/[[:punct:]]/-/g'`        # get name by delete ending, hyphens and punctuation from reference file
        echo -ne "\n$i additional reads could be assigned to $name\n"                                                       # user info
        if (( $i > 0 ))                                                                                                     # if at least one read was asigned, then ...
            then
                printf "$method\t$sktax\t$famtax\t$tid\t$name\t$identity\n%.0s" $(seq 1 $i) > zwischendurch1.txt            # print read assignment info as often are reads were assigned to reference
                paste zwischendurch1.txt asigned.txt > zwischendurch3.txt                                                   # paste info an reads
                mv zwischendurch3.txt asigned.txt                                                                           # rename file
                rm zwischendurch1.txt                                                                                       # remove tmp files
                cat asigned.txt >> ${arbeitsvz}/asignedAcc.txt                                                              # add reads to asignedAcc.txt
                get_rest_reads                                                                                              # get rest read (functions.sh)
        fi

        restreads=`wc -l ${arbeitsvz}/restAcc.txt | cut -d " " -f1`                                                         # count number of rest Reads
        echo -ne "and $restreads are still waiting"                                                                         # user info
    else
        name=`echo $organism | sed 's/.fna//' | sed 's/.*_//' | sed 's/-/ /g' | tr -d "." | sed 's/[[:punct:]]/-/g'`        # get name by delete ending, hyphens and punctuation from reference file
        echo -ne "\nno Reads could be assigned to $name"                                                                    # user info
fi
}

function getOutput()                                                                                                        # get output after blast
{                                                                                           
if [[ -s Blast-Hits.txt ]]                                                                                                  # if Blast-Hits were found, then ...
    then    
        sed '1,2d' Blast-Hits.txt | grep -v "NA-" > tmp.txt                                                                 # delete first to lines of Blast-Hits.txt (Header) and all "NA-" tax-ids and write to temporary file
        i=`wc -l < tmp.txt`                                                                                                 # count blast-hits
        if (( $i > 0 ))                                                                                                     # if at least one blast-hit was found, then ...
            then                                                                                                                
                printf "$method\n%.0s" $(seq 1 $i) > tmp1.txt                                                               # prepare temporatry file with method (repeat as often as hits were found)        
                paste tmp1.txt tmp.txt > read_method.txt                                                                    # and combine temporary blast-file and temporary method-file       
        fi      
        grep "NA-" Blast-Hits.txt > NA-hits.txt                                                                             # write all hits with "NA-" tax-id to separate file
        i=`wc -l < NA-hits.txt`                                                                                             # and count them
        head -n $i tmp1.txt > forNAHits.txt                                                                                 # take $i line of temporary method-file
        paste forNAHits.txt NA-hits.txt > NA-read_method.txt                                                                # and combine it with NA-hits
                
        cut -f 8 Blast-Hits.txt | sed '1,2d' | grep -v "NA-" > taxids.txt                                                   # get all tax-ids from Blast-Hits (column 8), exclude first two lines and NA hits              
        sort taxids.txt | uniq | sed '/^$/d' > taxids-uniq.txt                                                              # sort them uniq and remove empty rows
                
        while read line                                                                                                     # while reading taxids-uniq.txt, do ...
            do  
                if [[ `grep -w "^${line}" ${taxdir}/names.dmp | grep '\<scientific name\>'` != "" ]]                        # if "scientific name" for taxid can be grepped and is not an empty string, then ...
                    then
                        grep -w "^${line}" ${taxdir}/names.dmp | grep '\<scientific name\>' | cut -f3 | tr -d [=\'=],[=\(=],[=\)=],[=\[=],[=\]=] | tr -s [=\/=] "." | tr -s "-" " " | tr -d "." | sed 's/[[:punct:]]/-/g' >> names2.txt   # grep $line in names.dmp, grep "scientific name" (if possible) and cut 3rd column for name
                elif [[ `grep -w "^${line}" ${taxdir}/names.dmp | grep -m 1 '\<synonym\>'` != "" ]]                                                                                                # else, ...
                    then
                        grep -w "^${line}" ${taxdir}/names.dmp | grep -m 1 '\<synonym\>' | cut -f3 | tr -d [=\'=],[=\(=],[=\)=],[=\[=],[=\]=] | tr -s [=\/=] "." | tr -s "-" " " | tr -d "." | sed 's/[[:punct:]]/-/g' >> names2.txt      # grep the first synonym you find as a name
                    else
                        echo >> names2.txt
                fi
            done < taxids-uniq.txt
        
        paste taxids-uniq.txt names2.txt > names-taxid.txt                                                                  # combine tax-ids and names
                                
        while read line                                                                                                     # while reading the names-taxid-file, do...
            do
                id=`echo $line | cut -d " " -f1`                                                                            # get taxid (first column)
                name=`echo $line | cut -d " " -f2-`                                                                         # get name (beginning second column)
                awk -v tax="$id" '$9 == tax' read_method.txt > tmp.txt                                                      # get all lines that have the taxid in the 9th column
                num=`wc -l < tmp.txt`                                                                                       # and count lines
                if (( $num > 0 ))                                                                                           # if at least one hit was found, then ...
                    then
                        printf "$name\n%.0s" $(seq 1 $num) > tmp_name.txt                                                   # print name as often as hits were found
                        paste tmp.txt tmp_name.txt > ${id}_names.txt                                                        # paste hits and names per taxid
                fi
            done < names-taxid.txt
        cat *_names.txt > read_names_method.txt & wait                                                                      # cat all hit_names subsets
        rm *_names.txt                                                                                                      # and delete subsets
        awk -v var="$method" 'BEGIN {FS="\t"} {x=var".txt"} {print $1 FS $11 FS $10 FS $9 FS $12 FS $8 FS $2 > x}' < read_names_method.txt # change order of read_names_method.txt and write it to new file (name by previouly defined method)
        if [[ $method == Blastp_vs_protdb ]]                                                                                # if the method was blastp, then ... (needed because results are saves seperately)
            then
                cat ${method}.txt >> ${arbeitsvz}/asignedAcc_pep.txt                                                        # add results to asignedAcc_pep.txt
        elif [[ $method == Blastx_vs_protdb ]]                                                                              # if the method was blastp, then ... (needed because results are saves seperately)
            then
                cat ${method}.txt >> ${arbeitsvz}/asignedAcc_pep.txt
        elif [[ $method == tBlastx_vs_ntdb ]]                                                                               # if the method was blastp, then ... (needed because results are saves seperately)
            then
                cat ${method}.txt >> ${arbeitsvz}/asignedAcc_pep.txt
        elif [[ $method == ViRIEMS ]]
            then
                mv ${method}.txt ${arbeitsvz}/ViRIEMS_asigned.txt
            else                                                                                                            # else, ...
                cat ${method}.txt >> ${arbeitsvz}/asignedAcc.txt                                                            # add assigned Reads to asignedAcc.txt file
        fi
        get_rest_reads                                                                                                      # get rest reads (functions.sh)
fi
restreads=`wc -l ${arbeitsvz}/restAcc.txt | cut -d " " -f1`                                                                 # count number of rest reads
echo -ne "$restreads are still waiting to be assigned\n"                                                                    # user info
}
