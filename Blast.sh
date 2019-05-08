#    Blast.sh - part of RIEMS - Reliable Information Extraction from MEtagenomic Sequence datasets
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

. ${installdir}/Config.txt                                                                                          # <- Path of config file (unchanged)
. ${installdir}/TaxidDetermination.sh                                                                               # path to Tax functions
. ${installdir}/Blast_functions.sh
# Explanation of blastvek = initialisation vector for blast
# ${blastvek[*]}: map db retro/tax expl par
# ${blastvek[*]}: 0 0 0 0 0 for deaktivation of all functions

rank=family                                                                                                         # Default setting for rank
taxonomieiddb=$taxdir/gi_taxid_nucl.dmp                                                                             # Default setting for original (non-splittet) gi_taxid_nucl-file              
#cdhit=/home/software/CD-HIT                                                    
#identity=0.98                                                                   
if [ -s ${fasttaxidntdb}.dmp ]                                                                                      # if fasttaxid files exist (splitted gi_taxid_nucl), then...
    then                                                                                                        
        taxonomieiddb=$fasttaxidntdb                                                                                # set taxonomydb to fasttaxid to make the identification of taxids faster
fi                                                                                                              
                                                                                                                                                                          
################################################################## Erstellen der Blast Datenbank aller detektierten Orgs #########################################################

mkdir -p $blastordner                                                                                               # creat the blastfolder
cd $blastordner                                                                                                     # and change to blastfolder
                                           
query2=`basename $query .fna`                                                                                       # ${zielordner}/Partial10-89Reads.fna -> Partial10-89Reads
ln -s $query $query2 
query=$query2

################################################################## Blast + Hits Bearbeitung ##############################################################

anf=`grep -c '^>' $query`                                                                                           # count number of sequences in query                                 
echo -e "Accno \tRefsequence-ID \tHit-Range \t\tLength \tE-value \tperc.identity \tTax-ID \t$rank Tax-ID \tSuperkingdom Tax-ID \n" > Blast-Hits.txt   # create a file with header for blast-results
if [ ${blastvek[0]} != map ] && [ $referenz != $protdb ]                                                                                       # if the 1st parameter in blastvek is unequal to map, then ...
    then 
        echo -ne "BLAST ($ntblast) of $anf sequences vs nucleotide database... "                                    # user info
elif [ $referenz == $protdb ]
    then
        echo -ne "$anf protein sequences found\nBLAST ($ntblast) of $anf sequences vs protein database... "
fi
if [ $ntblast == blastp ]                                                              # if the reference is the protein database and the 2nd parameter in blastvek is unequal to db, then ....
    then 
        $blastdir/blastp -db $protdb -num_threads $threads -evalue $evalue -max_target_seqs 1 -out $blastordner/Hits.txt -query $query -max_hsps 1 -outfmt '6 qseqid sseqid qstart qend qlen evalue pident staxid length qcovs'
        #$blastdir/blastp -db $protdb -num_threads $threads -negative_gilist ${installdir}/gi_exclude.txt -evalue $evalue -max_target_seqs 1 -out $blastordner/Hits.txt -query $query -max_hsps 1 -outfmt '6 qseqid sseqid qstart qend qlen evalue pident staxid length qcovs'
        # perform blastp; negative gilist includes gis for environment samples and synthetic constructs, max_target_seqs and max_hsps have to be set to 1 to get only one result per sequence; output="tab-seperated(6) query sequence id; subtject sequence id; query start; query end; query length; evalue; identity; taxid"
elif [[ $ntblast == blastx ]]                                                                                       # if ntblast is set to tblastx, then ...
    then
        $blastdir/blastx -db $protdb -num_threads $threads -evalue $evalue -max_target_seqs 1 -out $blastordner/Hits.txt -query $query -max_hsps 1 -outfmt '6 qseqid sseqid qstart qend qlen evalue pident staxid'
        #$blastdir/blastx -db $protdb -num_threads $threads -negative_gilist ${installdir}/gi_exclude.txt -evalue $evalue -max_target_seqs 1 -out $blastordner/Hits.txt -query $query -max_hsps 1 -outfmt '6 qseqid sseqid qstart qend qlen evalue pident staxid'
        # perform blastx; negative gilist includes gis for environment samples and synthetic constructs, max_target_seqs and max_hsps have to be set to 1 to get only one result per sequence; output="tab-seperated(6) query sequence id; subtject sequence id; query start; query end; query length; evalue; identity; taxid"
elif [[ $ntblast == tblastx ]]                                                                                      # if ntblast is set to tblastx, then ...
    then
        $blastdir/tblastx -db $referenz -num_threads $threads -evalue $evalue -max_target_seqs 1 -out $blastordner/Hits.txt -query $query -max_hsps 1 -outfmt '6 qseqid sseqid qstart qend qlen evalue pident staxid'
        #$blastdir/tblastx -db $referenz -num_threads $threads -negative_gilist ${installdir}/gi_exclude.txt -evalue $evalue -max_target_seqs 1 -out $blastordner/Hits.txt -query $query -max_hsps 1 -outfmt '6 qseqid sseqid qstart qend qlen evalue pident staxid'
        # perform tblastx; negative gilist includes gis for environment samples and synthetic constructs, max_target_seqs and max_hsps have to be set to 1 to get only one result per sequence; output="tab-seperated(6) query sequence id; subtject sequence id; query start; query end; query length; evalue; identity; taxid"
    else                                                                                                            # otherwise, ...
        $blastdir/blastn -task $ntblast -db $referenz -num_threads $threads -evalue $evalue -max_target_seqs 1 -out $blastordner/Hits.txt -query $query -max_hsps 1 -outfmt '6 qseqid sseqid qstart qend qlen evalue pident staxid qcovs'
        #$blastdir/blastn -task $ntblast -db $referenz -negative_gilist ${installdir}/gi_exclude.txt -num_threads $threads -evalue $evalue -max_target_seqs 1 -out $blastordner/Hits.txt -query $query -max_hsps 1 -outfmt '6 qseqid sseqid qstart qend qlen evalue pident staxid qcovs'
        # perform blastn (megablast depending on -task); negative gilist includes gis for environment samples and synthetic constructs, max_target_seqs and max_hsps have to be set to 1 to get only one result per sequence; output="tab-seperated(6) query sequence id; subtject sequence id; query start; query end; query length; evalue; identity; taxid"

        if [[ $method != Blastn_vs_ntdb ]]
            then
                sortBlastResults
            else
                awk '$8 != "N/A"' Hits.txt > non-NA-Blastn.txt     # Remove hits that did not present with a taxonomic ID
                mv non-NA-Blastn.txt Hits.txt
        fi
fi
############################# probably not of use anymore ##############################

echo -ne "Data processing... "                                                                                      # user info

if [[ $method == Blastp_vs_protdb ]]                                                                                # if the method is a blastp, then ...
    then
        export arbeitsvz                                                                                            # export workingdir to environment
        export blastordner
        export ContigRead                                                                                           # export contant of ContigRead variable to environment (necessary for R-script)
        ${R}/Rscript --vanilla ${scriptBlastp} 12 &>>${arbeitsvz}/console.log                                       # run R-script to sort Blastp results (because more than one hit per sequence is possible) (see Blastp.R)
fi

cat $blastordner/Hits.txt | sort -n -t$'\t' -k 8 | grep -vw "NoHit" | cut -f1-7 > $blastordner/AllHits.txt          # get Hits and sort them by taxid (8th column, tab separated) and remove of NoHits (probably not needed) and 
                                                                                                                    # only cut column 1-7 because taxids will be added again later in the script to be sure species tax is detected

############################################################## Taxonomy Determination ####################################################################


if [ -s $blastordner/AllHits.txt ]                                                                                  # if the AllHits.txt file exists (blast was successful), then ...
    then                                                                                                        
        BlastHitTaxid                                                                                               # process data (see Blast_functions.sh)
fi
#rm $query 2>/dev/null
rm $blastordner/tid*  ; rm $blastordner/uniq-tids.txt  ; rm $blastordner/AllHits*.txt                          # remove files not of use anymore
