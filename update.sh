#!/bin/bash

. ~/RIEMS/Contig.txt

### build optimized nodes.dmp for fast retrieval of all TAX-ID succeeding underneath the species level determined before
### for this purpose, the nodes.dmp is split like the gi-taxid-nucl database
cut -d "|" -f1 /home/blast_db/taxonomy/nodes.dmp | tr -d [:blank:] > childTaxid.tmp             # extract the first column holding the child-TAX-ID from nodes.dmp and remove any whitespace
cut -d "|" -f2 /home/blast_db/taxonomy/nodes.dmp | tr -d [:blank:] > parentTaxid.tmp            # extract the second column holding the parent-TAX-ID from nodes.dmp and remove any whitespace

paste childTaxid.tmp parentTaxid.tmp | sort > childParentTaxid.tmp                              # fuse child-TAX-ID and parent-TAX-ID columns into a single file and sort ascending by child-TAX-ID

mkdir sorted
a=childParentTaxid.tmp
i=10                                                                                # Initialisierung i=10
until (( $i == 100 )) 
    do 
        echo $i >> zahlen.txt 
        ((i=$i+1)) 
    done  
while read line 
    do 
        OrgDmpDataSplitting $line & 
    done < zahlen.txt 
    wait       # Funktionsaufruf für Splitten der .dmp Datei in 89 Einzelteile
    


cp childParentTaxid.tmp temp1                                                                   # make a temporary copy of the child-to-parent file for the extraction of all lines with child-TAX-IDs beginning with defined two digits
i=10                                                                                            # initialize two-digit counter
while (( i < 100 ))                                                                             # while loop to repeat the extraction with increasing first two digits of the TAX-ID
do  
    grep "^${i}" temp1 > childParentTaxid_${i}.tmp                                              # extract all lines beginning with digits specified by i and write to new file 
    grep -v "^${i}" temp1 > temp2                                                               # remove the lines beginning with digits specified by i and write to new temporary file
    mv temp2 temp1                                                                              # rename temporary file for next use
    let i=i+1                                                                                   # increase counter
done                                                                                            # end while-loop
    
paste parentTaxid.tmp childTaxid.tmp | sort > parentChildTaxid.tmp                              # fuse parent-TAX-ID and child-TAX-ID columns into a single file and sort ascending by parent-TAX-ID
cp parentChildTaxid.tmp temp1                                                                   # make a temporary copy of the parent-to-child file for the extraction of all lines with parent-TAX-IDs beginning with defined two digits
i=10                                                                                            # initialize two-digit counter
while (( i < 100 ))                                                                             # while loop to repeat the extraction with increasing first two digits of the TAX-ID
do  
    grep "^${i}" temp1 > parentChildTaxid_${i}.tmp                                              # extract all lines beginning with digits specified by i and write to new file
    grep -v "^${i}" temp1 > temp2                                                               # remove the lines beginning with digits specified by i and write to new temporary file
    mv temp2 temp1                                                                              # rename temporary file for next use
    let i=i+1                                                                                   # increase counter
done                                                                                            # end while-loop

#########################################################################################################
#this part will split the gi_taxid_nucl.dmp file into taxid_gi_nucl.tmp files for a faster search


awk '{x="taxid_gi_nucl.dmp"} {print $2,$1 > x}' < gi_taxid_nucl.dmp                             # faster awk alternative to cut twice and paste subsequently

cp taxid_gi_nucl.tmp temp1                                                                      # make a temporary copy of the TAX-ID-gi file for the extraction of all lines with TAX-IDs beginning with defined two digits
i=10                                                                                            # initialize two-digit counter
while (( i < 100 ))                                                                             # while loop to repeat the extraction with increasing first two digits of the TAX-ID
do
    grep "^${i}" temp1 > taxid_gi_nucl_${i}.tmp                                                 # extract all lines beginning with digits specified by i and write to new file 
    grep -v "^${i}" temp1 > temp2                                                               # remove the lines beginning with digits specified by i and write to new temporary file
    mv temp2 /temp1                                                                             # rename temporary file for next use
    let i=i+1                                                                                   # increase counter
done                                                                                            # end while-loop
mv temp1 taxid_gi_nucl_0.tmp                                                                    # move all rest data (all tax-ids less than 10) to taxid_gi_nucl_0.tmp

for i in taxid_gi_nucl_*                                                                        #start loop for all taxid_gi_nucl files
    do
        sort ${i} > sort_${i}                                                                   #sort each taxid_gi_nucl file and save in sort
        mv sort_${i} ${i}                                                                       #move sort back to original file
    done





################################################################
#this subscript of the RIEMS workflow generates a gi-list to exclude during the analysis
#by default following Tax-ID are excluded: 48479, 28384, 12908

excludeTax=(48479 28384 12908)                                                                                                              # List of Tax-IDs to exclude... the List might be appended
out=${installdir}/taxtmp
mkdir -p ${out}
cd ${out}
for i in ${excludeTax[*]}
    do  
        echo $i
        parentTax=$i
        grep "^\<${parentTax}\>" ${fasttaxid}/sorted/parentChildTaxid_`echo $parentTax | cut -c1-2`.tmp | cut -f2 >> ${out}/tmptaxchild-${i}.txt # get all ChildTax-ids assigned to parentTax, only grep in file were tax-ids begin with first two characters of parenttax, get 2nd column (childTax)                        
        cat ${out}/tmptaxchild-${i}.txt >> ${out}/taxchild-${i}.txt                                                                                   # cat temporary taxchild to taxchild
        
        until ! [ -s ${out}/tmptaxchild-${i}.txt ]                                                                                               # until the tmptaxchild.txt is empty, do ...
            do
                while read line                                                                                                             # while reading tmptaxchild.txt, do ...
                    do
                        m=`echo $line | cut -c1-2`
                        if (( $m < 10 ))
                            then
                                grep "^\<${line}\>" ${fasttaxid}/sorted/parentChildTaxid_0.tmp | cut -f2 >> ${out}/tmp2taxchild-${i}.txt 
                            else
                                grep "^\<${line}\>" ${fasttaxid}/sorted/parentChildTaxid_`echo $line | cut -c1-2`.tmp | cut -f2 >> ${out}/tmp2taxchild-${i}.txt   # grep grep child-tax-ids in respective parentChildTaxid*.tmp, cut 2nd column (child-tax) and save in 2nd temporary file
                        fi
                    done < ${out}/tmptaxchild-${i}.txt 1>/dev/null 2>&1 & wait                                                                                        # in background to grep all at once and wait til they finished
                mv ${out}/tmp2taxchild-${i}.txt ${out}/tmptaxchild-${i}.txt                                                                           # move 2nd to 1st temporary file (now tmptaxchild is not empty)
                cat ${out}/tmptaxchild-${i}.txt >> ${out}/taxchild-${i}.txt                                                                           # cat tmptaxchild to taxchild.txt
            done
        echo $parentTax | cat >> ${out}/taxchild-${i}.txt                                                                                        # finally add parentTax to taxchild
        cat ${out}/taxchild-${i}.txt >> ${installdir}/exclude_tid.txt
        z=`wc -l < ${out}/taxchild-${i}.txt`                                                                                                     # count number of taxids
        if (( $z > 23 ))        
            then        
                (( p=$z/23 ))                                                                                                               # calculate p by dividing number of taxids by number of cores available
                split -l $p ${out}/taxchild-${i}.txt ${out}/part-${i}- ; wait                                                                         # split childtax.txt by p
            else        
                cp ${out}/taxchild-${i}.txt ${out}/part-${i}-aa       
        fi      
        for j in part-${i}-*                                                                                                                     # for each part-* file, do ...
            do      
                while read line                                                                                                             # while reading part-* file, do ...
                    do
                        k=`echo $line | cut -c1-2`
                        if (( $k < 10 ))
                            then 
                                grep "^\<${line}\>" ${fasttaxid}/tax_id_nucl/taxid_gi_nucl_0.tmp | cut -d " " -f2 >> ${out}/gi_${i}-${j}.txt
                            else
                                grep "^\<${line}\>" ${fasttaxid}/tax_id_nucl/taxid_gi_nucl_`echo $line | cut -c1-2`.tmp | cut -d " " -f2 >> ${out}/gi_${i}-${j}.txt  # grep taxid and associated gi in taxid_gi_nucl (only grep in file where taxids start with first to characters of $line), cut 2nd column (gi)
                        fi
                    done < ${j} &  1>/dev/null 2>&1                                                                                         # in background to start 23 processes at once
            done ; wait 1>/dev/null 2>&1
        wait
        cat gi_*-part-* > gi_TID-${i}.txt
        #rm tmptaxchild.txt ; rm taxchild.txt ; rm *part*
    done 
    
cat gi_TID* > ${installdir}/gi_exclude.txt
cd ${arbeitsvz}
rm -r ${installdir}/taxtmp

################################################################
#this subscript of the RIEMS workflow updates reference sequences in refseqdir

for k in ~/RIEMS/RefSequenzen/TID-*.fna
    do
        cp $k .
        k=`echo $k | sed 's/.*\///g'`
        grep "^>" $k | cut -f2 -d "|" | sort > existing_gis.txt
        tid=`echo $k | sed 's/_.*//g' | sed 's/.*-//g'`
        grep "^\<${tid}\>" ${fasttaxid}/parentChildTaxid/parentChildTaxid-`echo $tid | cut -c1-2`.dmp | cut -f2 >> tmptaxchild-${tid}.txt # get all ChildTax-ids assigned to parentTax, only grep in file were tax-ids begin with first two characters of parenttax, get 2nd column (childTax)                        
        cat tmptaxchild-${tid}.txt >>taxchild-${tid}.txt                                                                                   # cat temporary taxchild to taxchild
        
        until ! [ -s tmptaxchild-${tid}.txt ]                                                                                               # until the tmptaxchild.txt is empty, do ...
            do
                while read line                                                                                                             # while reading tmptaxchild.txt, do ...
                    do
                        m=`echo $line | cut -c1-2`
                        grep "^\<${line}\>" ${fasttaxid}/parentChildTaxid/parentChildTaxid-`echo $line | cut -c1-2`.dmp | cut -f2 >> tmp2taxchild-${tid}.txt   # grep grep child-tax-ids in respective parentChildTaxid*.tmp, cut 2nd column (child-tax) and save in 2nd temporary file
                    done < tmptaxchild-${tid}.txt 1>/dev/null 2>&1 & wait                                                                                        # in background to grep all at once and wait til they finished
                    mv tmp2taxchild-${tid}.txt tmptaxchild-${tid}.txt                                                                           # move 2nd to 1st temporary file (now tmptaxchild is not empty)
                    cat tmptaxchild-${tid}.txt >> taxchild-${tid}.txt                                                                           # cat tmptaxchild to taxchild.txt
            done
        echo $tid | cat >> taxchild-${tid}.txt
        rm tmptaxchild-*
        z=`wc -l < taxchild-${tid}.txt`                                                                                                     # count number of taxids
        if (( $z > $threads ))        
            then        
                (( p=$z/$threads ))                                                                                                               # calculate p by dividing number of taxids by number of cores available
                split -l $p $taxchild-${tid}.txt part-${tid}- ; wait                                                                         # split childtax.txt by p
            else        
                split -l 1 taxchild-${tid}.txt part-${tid}- ; wait
        fi  
        
        for j in part-${tid}-*                                                                                                                     # for each part-* file, do ...
            do      
                while read line                                                                                                             # while reading part-* file, do ...
                    do
                        k=`echo $line | cut -c1-2`
                        grep "^\<${line}\>" ${fasttaxid}/taxid_gi_nucl/taxid_gi_nucl-`echo $line | cut -c1-2`.dmp | cut -d " " -f2 >> gi_${tid}-${j}.txt  # grep taxid and associated gi in taxid_gi_nucl (only grep in file where taxids start with first to characters of $line), cut 2nd column (gi)
                    done < ${j} &  1>/dev/null 2>&1                                                                                         # in background to start 23 processes at once
            done ; wait 1>/dev/null 2>&1
        wait
        cat gi_*-part-* | sort > gi_TID-${tid}_all.txt
        rm gi_*-part-*; rm part-${tid}-*
        comm -3 existing_gis.txt gi_TID-${tid}_all.txt | sed -e 's/^[ \t]*//' > gi_TID-${tid}_new.txt
        z=`wc -l < gi_TID-${tid}_new.txt`
        if (( $z > 23 )) 
            then 
                (( p=$z/23 ))                                                                                                       # calculate p by dividing number of taxids by number of cores available
                split -l $p gi_TID-${tid}_new.txt gi_part- ; wait  1>/dev/null 2>&1                                                                 # split childtax.txt by p
            else
                split -l 1 gi_TID-${tid}_new.txt gi_part- ; wait  1>/dev/null 2>&1
        fi
        (for i in gi_part* ; do GiDownload $i & done 1>/dev/null 2>&1; wait) 
        
       
        cat gi_part*.fna $k > tmp.fna
        mv tmp.fna $k.new
        rm gi_* 2>/dev/null
        rm taxchild*.txt ; rm existing_gis.txt
         DoubbleKill_update
        echo -e "$k\t`date`" >> date.txt
        rm $k.new ; rm $k.kill
        #rm ~/RIEMS/reftmp/*
    done


function DoubbleKill_update()                                                                                      # removes duplicated gis after GiDownload
{
csplit -f split $k.new '/^>/' {*} 1>/dev/null 2>&1 & wait 1>/dev/null 2>&1           # split the sequence file for Tax-id by header so each file contains one sequence

i=0                                                                                                         # set counter to 0
(while (( i < 10 )); do md5sum split*${i} >> md5.txt & let i=i+1 ; done 1>/dev/null 2>&1 ; wait)   # calculate md5sum for each split file and save info in md5.txt # wait for all background processes to finish
cut -d " " -f3 md5.txt | sort > all.txt                                                       # cut the 3rd column of md5.txt (two spaces between columns) to get all filenames
sort md5.txt | uniq -w32 | cut -d " " -f3 | sort  > uniq.txt                                  # sort md5.txt and make uniq by comparing the first 32 characters (md5sum) and cut filenames again and sort

comm -13 uniq.txt all.txt | sed -e 's/^[ \t]*//' > to_remove.txt                       # compare uniq files and all files and get all files that are only in all (-13), delete tabs at beginning of line and write them to to_remove
comm -12 uniq.txt all.txt | sed -e 's/^[ \t]*//' > to_keep.txt                         # compare uniq files and all files and get all files that are in both lists (-12), delete tabs at beginning of line and write them to to_keep
                                                                                                  # change working directory to ${out}
xargs rm < to_remove.txt 2>/dev/null                                                                 # remove all files which are listed in to_remove, xargs accepts a list of arguments and runs command (here rm) for this list  
                                                                                                            # faster than rm alone and can also use large lists of arguments, only rm might give error for to long list
j=0                                                                                                         # set counter back to 0
while (( j < 10 ))                                                                                          # as long as i is less than 10, do ...
    do 
        cat split*${j}                                                                               # cat all split files that end with $i
        echo                                                                                                # echo to start new line after each cat
        let j=j+1                                                                                           # add 1 to i
    done >>  $k.kill 2>/dev/null                                                 # write everythink into one file in Refseqdir

sed '/>/!s/[KMRYSWBVHDNkmryswbvhdn]//g'  $k.kill > $k.rename     #remove all ambiguities in reference sequences except in lines containing an ">" (header)
    
all=`grep "^>" $i | sort | wc -l`                                    # grep all header, sort them and count the lines (equal to count all header)
uniq=`grep "^>" $i | sort | uniq | wc -l`                            # grep all header, sort them uniq and count the lines (equal to count all uniq header)
if (( $all > $uniq ))                                                                                       # if number of all header is greater than number of uniq header, then ...
    then                                                                                                    
        echo "DAMN IT"                                                                                      # user info 
fi

xargs rm <  to_keep.txt 2>/dev/null                                                                   # remove all split files that are left using xargs (see above)
rm  md5.txt;  rm  all.txt; rm  uniq.txt                                                   # remove files that will not be needed anymore
rm  to_remove.txt; rm  to_keep.txt 2>/dev/null                                                  # remove files that will not be needed anymore

#cd ${arbeitsvz}                                                                                             # change to working directory
#rm -r ${installdir}/reftmp                                                                                  # and remove reftmp folder
}









