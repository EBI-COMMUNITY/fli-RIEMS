#!/bin/bash

#    ViRiems.sh - part of RIEMS - Reliable Information Extraction from MEtagenomic Sequence datasets
#    Copyright (C) 2009-2014  Matthias Scheuch, Dirk Hoeper
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




echo -ne "\nViRIEMS workflow was started.\nAll sequences will be first screened against a virus database. The whole RIEMS analysis will be started, subsequently."
method=ViRIEMS
    blastvek=(0 0 tax 0 0)                                                                                      # set blastvector for blast.sh              
    query=${arbeitsvz}/TrimmedReads.fasta                                                                       # set query for blast   
    referenz=${fasttaxidvirus}                                                                                            # set reference for blast 
    blastordner=${arbeitsvz}/MultiBlast/ViRiems                                                               # define blast-folder  
    ntblast=megablast                                                                                           # set further parameters for blastn (e.g. set blast to megablast)
    . ${installdir}/Blast.sh                                                                                    # call blast.sh
    getOutput

rm -r ${arbeitsvz}/MultiBlast/ViRiems/*
cut -f7 ${arbeitsvz}/ViRIEMS_asigned.txt > ${arbeitsvz}/MultiBlast/ViRiems/ViRIEMS.acc
${gsflxdir}/fnafile -i ${arbeitsvz}/MultiBlast/ViRiems/ViRIEMS.acc -o ${arbeitsvz}/MultiBlast/ViRiems/ViRIEMS.fna ${arbeitsvz}/TrimmedReads.fasta

method=ViRIEMS
    blastvek=(0 0 tax 0 0)                                                                                      # set blastvector for blast.sh              
    query=${arbeitsvz}/MultiBlast/ViRiems/ViRIEMS.fna                                                                       # set query for blast   
    referenz=${ntdb}                                                                                            # set reference for blast 
    blastordner=${arbeitsvz}/MultiBlast/ViRiems                                                               # define blast-folder  
    ntblast=megablast                                                                                           # set further parameters for blastn (e.g. set blast to megablast)
    . ${installdir}/Blast.sh                                                                                    # call blast.sh
    getOutput

cd ${arbeitsvz}

awk '$2 == 10239' ${arbeitsvz}/ViRIEMS_asigned.txt > ${arbeitsvz}/ViRIEMS_asigned_viral.txt
mv ${arbeitsvz}/ViRIEMS_asigned_viral.txt ${arbeitsvz}/ViRIEMS_asigned.txt
if [[ -s ${arbeitsvz}/ViRIEMS_asigned.txt ]]
    then
        cut -f 3 ${arbeitsvz}/ViRIEMS_asigned.txt > ${arbeitsvz}/ViRIEMS_famtax.txt                                                  # cut the 3rd column of asignedAcc.txt to get all famtax that were detected (necessary for resultprotocol) 
        sort -n ${arbeitsvz}/ViRIEMS_famtax.txt | uniq | sed '/^$/d' > ${arbeitsvz}/ViRIEMS_uniq_famtax.txt                             # sort, uniq and delete all empty lines in famtax.txt
        while read line                                                                                                 # while reading uniq_famtax.txt
            do
                famtax=$line                                                                                            # assign line to famtax
                    if [[ $famtax == NA ]]                                                                              # if famtax is NA (no family found), then ...
                        then 
                            echo "NA" >> ${arbeitsvz}/ViRIEMS_names.txt                                                         # print "NA" to names.txt, else ...
                        else
                            grep "^\<$famtax\>" ${taxdir}/names.dmp | grep '\<scientific name\>' | cut -f3 >> ${arbeitsvz}/ViRIEMS_names.txt    # grep famtax in names.dmp and grep "scientific name" (3rd column) and write it to names.txt
                    fi
            done < ${arbeitsvz}/ViRIEMS_uniq_famtax.txt
        paste ${arbeitsvz}/ViRIEMS_uniq_famtax.txt ${arbeitsvz}/ViRIEMS_names.txt > ${arbeitsvz}/ViRIEMS_famtax_names.txt                       # paste famtax and names to get a table for the resultprotocol
        
        export arbeitsvz
        
        echo -ne "\nResult tables for ViRIEMS will be created.\n"
        
        ${R}/Rscript --vanilla ${installdir}/ViRIEMS.R 12 &>>${arbeitsvz}/console.log
        
        blastv=`${blastdir}/blastn -version | tail -n1`                                                                   # get blastversion
        blastdbv=`${blastdir}/blastdbcmd -db /home/blast_db/ncbi/nt -info | grep "Date" | cut -f1 | sed 's/Date: //g'`              # get version of blast database
        texversion=`${latexdir}/tex --version | head -n 1`
        date3=`date`
        export blastdbv
        export blastv
        export texversion
        export projektname
        export date1
        export date3
        
        echo -ne "\nA PDF-file with ViRIEMS results will be created.\n"
        
        
        
        ${R}/Rscript -e "library(knitr); knit('${installdir}/ViRIEMS.Rnw')" 12 &>>${arbeitsvz}/console.log
        
        pdflatex ${arbeitsvz}/ViRIEMS.tex 12 &>>${arbeitsvz}/console.log
        pdflatex ${arbeitsvz}/ViRIEMS.tex 12 &>>${arbeitsvz}/console.log
fi

