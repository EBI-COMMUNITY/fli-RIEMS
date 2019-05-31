#!/usr/bin/env bash

#   ScaffoldBlast.sh
#   Author: Nadim Rahman


. ${installdir}/Config.txt
. ${installdir}/functions_megablast_contig_reads.sh
. ${installdir}/TaxidDetermination.sh


# Setup
mkdir ${arbeitsvz}/trimmed_reads
mkdir ${arbeitsvz}/spades_assembly
mkdir ${arbeitsvz}/bowtie_readtoscaffold


# Trimmomatic quality check
START_01=$(date +%s)
echo -ne "\n[`date`] Trimming RIEMS input with trimmomatic..."

if [[ -n ${unspecified_input} ]]
    then
        echo -ne "${unspecified_input[@]}"
fi #Quality Mapping against dummy-reference without adapter trimming

if [[ -n ${illumina_input} ]]
    then
        parentName=`basename ${illumina_input[0]} | cut -d '_' -f 1`    # Get the name of parent name of the fastq file (e.g. ERR2023578 from ERR2023578_1.fastq)
        trimmomatic PE ${illumina_input[0]} ${illumina_input[1]} -baseout ${arbeitsvz}/trimmed_reads/${parentName}.fastq ILLUMINACLIP:${trimmomaticdir}adapters/TruSeq2-PE.fa:2:30:10 >/dev/null 2>/dev/null
fi #Quality Mapping of Illumina data input with adapter trimming

if [[ -n ${ionTorrent_input} ]]
    then
        echo -ne "${ionTorrent_input[@]}"
fi #Quality Mapping of IonTorrent data input with adapter trimming
echo -ne "\n[`date`] Trimming RIEMS input with trimmomatic... COMPLETE"
END_01=$(date +%s)
DIFF_01=$(( $END_01 - $START_01 ))
echo -ne "\nTIMING\t${DIFF_01}\tTrimmomatic trimming RIEMS input"


# SPAdes assembly
echo -ne "\n[`date`] Carrying out genome assembly using SPAdes with trimmed reads..."
START_02=$(date +%s)
${spadesdir}/spades.py -1 ${arbeitsvz}/trimmed_reads/${parentName}_1P.fastq -2 ${arbeitsvz}/trimmed_reads/${parentName}_2P.fastq -o ${arbeitsvz}/spades_assembly >/dev/null 2>/dev/null
cp ${arbeitsvz}/spades_assembly/scaffolds.fasta ${arbeitsvz}
cd $arbeitsvz
grep ">" scaffolds.fasta | tr -d ">" > scaffoldsInfo.txt        # Obtain a list of all scaffolds
echo -ne "\n[`date`] Carrying out genome assembly using SPAdes with trimmed reads... COMPLETE. Obtained a list of scaffolds."
END_02=$(date +%s)
DIFF_02=$(( $END_02 - $START_02 ))
echo -ne "\nTIMING\t${DIFF_02}\tSPAdes genome assembly with trimmed reads"


# Megablast scaffolds
START_03=$(date +%s)
echo -ne "\n[`date`] Megablast scaffolds generated from assembly."
${blastdir}/blastn -db ${ntdb} -task megablast -num_threads 40 -evalue 0.001 -max_hsps 1 -max_target_seqs 1 -out Blast-Hits.txt -query scaffolds.fasta -outfmt '6 qseqid sseqid pident staxid qcovs'    # blast contigs and get query-id and hit-id tab-seperated
END_03=$(date +%s)
DIFF_03=$(( $END_03 - $START_03 ))
echo -ne "\nTIMING\t${DIFF_03}\tMegablast command for assembly"

START_04=$(date +%s)
awk '$5 > 70' Blast-Hits.txt > sorted-Hits.txt
cut -f1-4 sorted-Hits.txt > Blast-Hits.txt
awk '$4 != "N/A"' Blast-Hits.txt > non-NA-Hits.txt
mv non-NA-Hits.txt Blast-Hits.txt

mv Blast-Hits.txt tmp.txt
while read line
    do
        grep ${line} tmp.txt >> Blast-Hits.txt
    done < scaffoldsInfo.txt

cut -f1 Blast-Hits.txt | sort | uniq -u > assigned_scaffolds.txt        # Include a unique list of scaffold names which have a BLAST hit
END_04=$(date +%s)
DIFF_04=$(( $END_04 - $START_04 ))
echo -ne "\nTIMING\t${DIFF_04}\tPost processing BLAST hits"


# Mapping reads back to assembly
START_05=$(date +%s)
echo -ne "\n[`date`] Mapping reads back to the scaffold assembly."

# Remove spaces in the headers for each read
perl -p -i.bk -e '~s/ /\//g;' ${arbeitsvz}/trimmed_reads/${parentName}_1P.fastq
perl -p -i.bk -e '~s/ /\//g;' ${arbeitsvz}/trimmed_reads/${parentName}_2P.fastq

cp scaffolds.fasta ${arbeitsvz}/bowtie_readtoscaffold
cd ${arbeitsvz}/bowtie_readtoscaffold
bowtie2-build --threads ${threads} scaffolds.fasta scaffolds >/dev/null 2>/dev/null    # Required to build an index of scaffolds.fasta to run bowtie mapping
# Bowtie command:
#   --> -x is the base name for all the scaffolds indexes bowtie2-build generated
#   --> Output generated is SAM format
bowtie2 -x scaffolds -p ${threads} -1 ${arbeitsvz}/trimmed_reads/${parentName}_1P.fastq -2 ${arbeitsvz}/trimmed_reads/${parentName}_2P.fastq -S bowtie_output.txt
grep -v '^@SQ' bowtie_output.txt > read_scaffold_alignment.txt     # Remove the headers from output (lines begin with @)
awk '! ( $3 ~ /\*/ )' read_scaffold_alignment.txt | cut -f 1,2,3,5 > reads_to_scaffold_hits_only.txt       # Only include reads which mapped to a scaffold (* = no result was found in second column)
awk '( $3 ~ /\*/ )' read_scaffold_alignment.txt | cut -f 1,2,3,5 > non_scaffold_mapped_reads.txt       # Presents specific bowtie results for reads which did not map to a scaffold
END_05=$(date +%s)
DIFF_05=$(( $END_05 - $START_05 ))
echo -ne "\nTIMING\t${DIFF_05}\tMapping reads back to assembly"


# Assign a species to each read from the assembly
# First ensure that only reads from assigned scaffolds (scaffolds with BLAST hits) are included
START_06=$(date +%s)
echo -ne "\n[`date`] Preparing for taxonomic assignment of reads for each scaffold."
cp ${arbeitsvz}/bowtie_readtoscaffold/reads_to_scaffold_hits_only.txt $arbeitsvz
cp ${arbeitsvz}/bowtie_readtoscaffold/non_scaffold_mapped_reads.txt $arbeitsvz
cd $arbeitsvz
while read line
    do
        grep ${line} reads_to_scaffold_hits_only.txt >> reads_in_assigned_scaffolds.txt
    done < assigned_scaffolds.txt

get_tid     # Get the taxonomic IDs
scaffolds_to_tid        # Get a list of the taxonomic IDs for the scaffolds
END_06=$(date +%s)
DIFF_06=$(( $END_06 - $START_06 ))
echo -ne "\nTIMING\t${DIFF_06}\tPreparations for read taxonomic assignment"

echo -ne "\n[`date`] Obtaining taxonomic information for the taxonomic IDs identified in BLAST hits, prior to assignment."
START_PRETAX=$(date +%s)
split -l 1 uniq-tids.txt tid-part-
for i in tid-part-*
    do
        while read line
            do
                tid=${line}
                get_species
                FamSkTaxDetermination

                if [[ `grep -w "^${tid}" ${taxdir}/names.dmp | grep '\<scientific name\>'` != "" ]]          # if "scientific name" for taxid can be grepped and is not an empty string, then ...
                    then
                        name=`grep -w "^${tid}" ${taxdir}/names.dmp | grep '\<scientific name\>' | cut -f3 | tr -d [=\'=],[=\(=],[=\)=],[=\[=],[=\]=] | tr -s [=\/=] "." | tr -s "-" " " | tr -d "." | sed 's/[[:punct:]]/-/g'`    # grep $line in names.dmp, grep "scientific name" (if possible) and cut 3rd column for name
                elif [[ `grep -w "^${tid}" ${taxdir}/names.dmp | grep -m 1 '\<synonym\>'` != "" ]]                                                                                                # else, ...
                    then
                        name=`grep -w "^${tid}" ${taxdir}/names.dmp | grep -m 1 '\<synonym\>' | cut -f3 | tr -d [=\'=],[=\(=],[=\)=],[=\[=],[=\]=] | tr -s [=\/=] "." | tr -s "-" " " | tr -d "." | sed 's/[[:punct:]]/-/g'`       # grep the first synonym you find as a name
                    else
                        name=`echo`                                                         # if no name is found write and empty line (important because otherwise the order is messed up)
                fi
                printf "${sktax}\t${famtax}\t${name}\t${tid}\n%.0s" >> ${arbeitsvz}/identifiedTaxClassifications.txt         # Save the taxonomic information retrieved in a file
            done < $i &
    done; wait
rm tid-part*

# Remove duplicate taxonomic ID information
sort ${arbeitsvz}/identifiedTaxClassifications.txt | uniq > ${arbeitsvz}/tmpTaxClass.txt
mv ${arbeitsvz}/tmpTaxClass.txt ${arbeitsvz}/identifiedTaxClassifications.txt
END_PRETAX=$(date +%s)
DIFF_PRETAX=$(( $END_PRETAX - $START_PRETAX ))
echo -ne "\nTIMING\t${DIFF_PRETAX}\tObtaining taxonomic information prior to assignment"


echo -ne "\n[`date`] Running taxonomic assignment of reads for each assigned scaffold, this may take some time..."
if [[ -s scaffold_tid.txt ]]
    then
        linenum=`wc -l < scaffold_tid.txt`      # Number of scaffolds with taxonomic IDs
        split -l 1 Blast-Hits.txt scaffoldBlastHits_part-
        START_07=$(date +%s)
        for i in scaffoldBlastHits_part-*
            do
                while read line
                    do
                        #START_VARS=$(date +%s)
                        scaffold_name=`echo $line | cut -f 1 -d " "`        # Get the name of the scaffold
                        ident=`echo $line | cut -f 3 -d " "`     # Get the BLAST identity
                        tid=`echo $line | cut -f 4 -d " "`      # Get the taxonomic ID
                        taxInfoCheck=`grep -m 1 "${tid}$" ${arbeitsvz}/identifiedTaxClassifications.txt`
                        #END_VARS=$(date +%s)
                        #DIFF_VARS=$(( $END_VARS - $START_VARS ))
                        #echo -ne "\nTIMING\t${DIFF_VARS}\tDefining variables"

                        if [[ "${taxInfoCheck}" != "" ]]
                            then
                                #echo -ne "\n[$(date)] Taxonomic ID already has been identified, obtaining stored taxonomic information."
                                #START_TAXINFO=$(date +%s)
                                sktax=`echo ${taxInfoCheck} | awk '{ print $1 }'`
                                famtax=`echo ${taxInfoCheck} | awk '{ print $2 }'`
                                name=`grep -m 1 "${tid}$" ${arbeitsvz}/identifiedTaxClassifications.txt | cut -f3`         # Required grep to obtain the whole name as space was causing issue
                                #END_TAXINFO=$(date +%s)
                                #DIFF_TAXINFO=$(( $END_TAXINFO - $START_TAXINFO ))
                                #echo -ne "\nTIMING\t${DIFF_TAXINFO}\tObtaining taxonomic information that already existed"
                            else
                                echo -ne "\n[$(date)] Identified a new taxonomic ID, obtaining information and storing this information."
                                #START_TAX=$(date +%s)
                                get_species
                                FamSkTaxDetermination
                                #END_TAX=$(date +%s)
                                #DIFF_TAX=$(( $END_TAX - $START_TAX ))
                                #echo -ne "\nTIMING\t${DIFF_TAX}\tObtaining species and other taxonomic information"

                                #START_SCINAME=$(date +%s)
                                if [[ `grep -w "^${tid}" ${taxdir}/names.dmp | grep '\<scientific name\>'` != "" ]]          # if "scientific name" for taxid can be grepped and is not an empty string, then ...
                                    then
                                        name=`grep -w "^${tid}" ${taxdir}/names.dmp | grep '\<scientific name\>' | cut -f3 | tr -d [=\'=],[=\(=],[=\)=],[=\[=],[=\]=] | tr -s [=\/=] "." | tr -s "-" " " | tr -d "." | sed 's/[[:punct:]]/-/g'`    # grep $line in names.dmp, grep "scientific name" (if possible) and cut 3rd column for name
                                elif [[ `grep -w "^${tid}" ${taxdir}/names.dmp | grep -m 1 '\<synonym\>'` != "" ]]                                                                                                # else, ...
                                    then
                                        name=`grep -w "^${tid}" ${taxdir}/names.dmp | grep -m 1 '\<synonym\>' | cut -f3 | tr -d [=\'=],[=\(=],[=\)=],[=\[=],[=\]=] | tr -s [=\/=] "." | tr -s "-" " " | tr -d "." | sed 's/[[:punct:]]/-/g'`       # grep the first synonym you find as a name
                                    else
                                        name=`echo`                                                         # if no name is found write and empty line (important because otherwise the order is messed up)
                                fi
                                printf "${sktax}\t${famtax}\t${name}\t${tid}\n%.0s" >> ${arbeitsvz}/identifiedTaxClassifications.txt         # Save the taxonomic information retrieved in a file
                                #END_SCINAME=$(date +%s)
                                #DIFF_SCINAME=$(( $END_SCINAME - $START_SCINAME ))
                                #echo -ne "\nTIMING\t${DIFF_SCINAME}\tObtaining scientific name to use"
                        fi

                        # Get the number of reads that made up the scaffold
                        #START_READS=$(date +%s)
                        reads=`grep ${scaffold_name} reads_in_assigned_scaffolds.txt | cut -f 1-2,4 > reads_for_${scaffold_name}_scaffold.txt`      # Obtains the reads and their assembly scores
                        length=`wc -l < reads_for_${scaffold_name}_scaffold.txt`
                        #END_READS=$(date +%s)
                        #DIFF_READS=$(( $END_READS - $START_READS ))
                        #echo -ne "\nTIMING\t${DIFF_READS}\tObtaining reads for given scaffold"

                        # For each of the reads used within the scaffold, create a line and append assignment information
                        #START_DEF=$(date +%s)
                        if (( $length > 0 ))
                            then
                                method=Scaffold_Assembly
                                printf "${method}\t${sktax}\t${famtax}\t${tid}\t${name}\t${ident}\n%.0s" $(seq 1 $length) > ${scaffold_name}_info.txt
                                paste ${scaffold_name}_info.txt reads_for_${scaffold_name}_scaffold.txt > ${scaffold_name}_completeinfo.txt
                        fi
                        #END_DEF=$(date +%s)
                        #DIFF_DEF=$(( $END_DEF - $START_DEF ))
                        #echo -ne "\nTIMING\t${DIFF_DEF}\tStoring results"
                    done < $i &
            done; wait
        END_07=$(date +%s)
        DIFF_07=$(( $END_07 - $START_07 ))
        echo -ne "\nTIMING\t${DIFF_07}\tTaxonomic assignments"
        echo -ne "\n[`date`] Completed running taxonomic assignments."

        START_08=$(date +%s)
        cat *_completeinfo.txt > assembledReads.txt
        cat assembledReads.txt >> asignedAcc.txt
        cat non_scaffold_mapped_reads.txt | cut -f 1 > restAccScaffold.txt      # All non-mapped reads, which cannot therefore be assigned (continues into the remaining workflow)
        completed=`wc -l < asignedAcc.txt`
        remaining=`wc -l < restAccScaffold.txt`

        # Remove all unimportant files
        rm scaffoldBlastHits_part-*
        rm *_info.txt
        rm *_completeinfo.txt
        rm reads_for_*_scaffold.txt
        rm assembledReads.txt non_scaffold_mapped_reads.txt
        END_08=$(date +%s)
        DIFF_08=$(( $END_08 - $START_08 ))
        echo -ne "\nTIMING\t${DIFF_08}\tConsolidating and tidying up taxonomic assignments"

        echo -ne "\n[`date`] Completed taxonomic assignment of BLAST-identified reads."
        echo -ne "\n[`date`] ${completed} reads have been successfully assigned. ${remaining} reads remaining."
        echo -ne "$remaining\t$method\n" >> ${arbeitsvz}/report.txt                                                     # write restreads and method to the report.txt
fi