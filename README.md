# fli-RIEMS

RIEMS - Reliable Information Extraction from Metagenomic Sequence datasets
    Copyright (C) 2009-2016  Ariane Belka, Maria Jenckel, Matthias Scheuch, Dirk Hoeper

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


   *** THANK YOU FOR USING RIEMS ON YOUR PC ***

 Please follow the instructions to complete the installation!



 EXPLANATION

 The software pipeline 'RIEMS' (Reliable Information Extraction of Metagenomic Sequence data sets)
 is a sorting tool written in Bash for classifying sequence reads according to taxa.
 


 IMPORTANT NOTES
 
 1. You are allowed to rename the folder 'RIEMS' before completing the installation.
 
 2. Do not change the names of any files in 'RIEMS'.
 
 3. Do not change the content of any files except during the installation.
 
 4. Do not add any files or folders to the present ones into the 'RIEMS' directory.
 
 5. 'RIEMS' requires the following software: 

        454 Genome Sequencer (version 3.0), 
        EMBOSS (version 6.3.1 or higher), 
        BLAST (version 2.4.0+ or higher),
        R (version 3.3.0 or higher) including Packages
            knitr
            ggplot2
            wordcloud
            xtable
            SciencesPo
            plyr
        Tex (version 3.141592 or higher) including Packages
            pdflscape
            colortbl
            color
            booktabs
            longtable
            array
            lastpage
            fancyhdr
            amssymb
            geometry
            inputenc
           
 6. 'RIEMS' requires the following databases: 
        --> database will be downloaded by the ncbiDBupdate.sh script
        BLAST nucleotide sequence database (using GIs and not Accessions),
        BLAST non-redundant protein sequence database (using GIs and not Accessions) and
        NCBI taxonomy database (using GIs and not Accessions)
                         
 7. 'RIEMS" requires the followring files:
        taxid_gi_nucl*
        --> the files can be created using the SplitDmp_VirGisTaxEx_Db1.sh tool

 8. Files to be analysed should not have any duplicate accession numbers.
    The assignments will be correct but the number of reads assigned to organisms
    will be imprecisely in the result protocol file.
    
 9. 'RIEMS' was evaluated under the software versions mentioned in 5. 
    The compatibility to other versions cannot be ensured. Please make sure when using
    BLAST versions other than 2.4.0+ that these put out the necessary result columns. 

 

 SET UP 'RIEMS'

 1. Move the folder 'RIEMS' to your desired directory and rename 'RIEMS' if desired.
 
 2. Open 'Config.txt' and adapt the paths-names (path-name/) of software and databases.

 3. Save and close 'Config.txt'.

 4. 'Riems.sh' has to read in the file 'Config.txt' before getting started. 
    For that purpose, in this file, the line displaying the config path (and the next two lines)
    must be replaced with the actual path of the 'path-name/Config.txt' file. 
    The line is situated in the first few lines of 'Riems.sh' file and marked with: 

    '# <- Path of config file'

    Do not change the '. ' at the beginning of this lines.
    
 5. Save and close the 'Riems.sh' file.

 6. SET UP completed!

 
 
 GETTING STARTED

 Enter 'sh RIEMS_installation_path/Riems.sh' without any quotes on your
 bash command line and follow the instructions

 or

 Enter 'sh RIEMS_installation_path/Riems.sh' and add the following parameters
 to adjust the analysis:

-j/i/t obligatory:  Specify sequence input file including full path
                    'RIEMS' will only accept fasta-, fastq- and sff-file formats
                    -i Illumina input data
                    -t Ion Torrent input data
                    -j unspecified input data
                    It is important to specify your input data because an adapter trimming is performed, 
                    if you have no information about your data or used a sequencing platform beside Illumina or Ion Torrent please
                    use the unspecified input (-j) option.
                    It is possible to add multiple sequence input files of multiple sequence platforms.
                    If using multiple sequence input files, please specify each data file separately.                                     

-o obligatory:      Enter target directory for the output (full path)                                                  

-p optional:        Project name, if not given project name will be named by input file name.
                    If multiple input files were selected the first one is used as project name.

-x optional:        Pre-Screening by given TaxID
                    multiple TaxID can be given but must be seperated by ",".
                    Please add TaxID on species level or below. TaxIDs above species level can not be used by RIEMS.

-r optional:        Specify if a resultreport should be created (n, y). Default is y

-f optional:        Further analysis (n, y), unclassified reads will be classified on protein level.
                    Default is n

-v optional:        ViRIEMS analysis (n, y), reads will be screened for viral sequences prior to RIEMS-Analysis.
                    Default is n
                    
Example 1:  Minimum of parameters needed

    sh path-name/RIEMS/Riems.sh -j path/filename.fna -o path/outputdir
    
    
Example 2:  Using multiple Illumina files for one analysis and add a project name

    sh path-name/RIEMS/Riems.sh -i path/filename1.fastq -i path/filename2.fastq -o path/outputdir -p ProjectX
    
    
Example 3:  Add multiple data files of different sequence platforms, use ViRIEMS and the further Analysis and add TaxIDs for Prescreenings        

    sh path-name/RIEMS/Riems.sh -j path/filename1.fna -i path/filename2.fastq -o path/outputdir -p ProjectX -x 9606,562 -f y -v y

    
    
Note:               Please to not use underlines "_". If you use them anyway they will be substituted to hyphens "-" 

 
 BUG FIXES

 28.08.2013 Input of several files enabled.

 29.08.2013 Bug fixed, caused by an incomplete EMBOSS link.

 25.09.2013 A failed Assembly is now able to be circumvented by an automatically restart.
  
 26.09.2013 Reads assigned to contigs in the assembly step are now listed under 'Megablast

            vs orgs' in the result protocol and not under 'Assembly' as before.

 09.10.2013 Initial sequence vector trimming implemented.

 11.10.2013 Implementing of a multithreading BLAST of partial no hit sequences.

 14.10.2013 Vector trimming implemented.

 14.11.2013 Partial Blast setting now available on command line settings.

 14.11.2013 Additional -n parameter added to mapping to avoid trimming.
 
 
 IMPROVEMENTS RIEMS 4.0
 
 - change output format to tab-seperated tables (resultprotocol-complete.txt, resultprotocol-compact.txt, resultprotocol-complete-pep.txt, resultprotocol-compact-pep.txt,)
 - Result report in PDF format (basic analysis and further analysis)
 - included Pre-Screening by TaxID (on species level)
 - specification of input by parameter (Illumina, IonTorrent, Unspecified)
 - adapted to larger datasets to run in appropriate time
  
   
   *** GOOD LUCK, HAVE FUN AND ENJOY THE RESULTS ***
