#!/bin/bash

anfang=`date`
zeit=`echo ${anfang} | tr -d " "`     # Datum fuer Verzeichnisnamen erzeugen

echo -e "\nProcessing of important NCBI database files. (duration ca. 30min)\n"

########################################################################## Dirks Eingaben ################################################################

arbeitsvz=/home/blast_db/taxLok.tmp
  # ist am Schluss auch gleichzeitig der Ordner mit allen relevanten Daten und ein Zugriffsort meines Metagenomanalyseprogramms
  
########################################################################### Definitionen #################################################################

array=(gi_taxid_nucl gi_taxid_prot taxid_gi_nucl childParentTaxid parentChildTaxid taxid_acc_nucl taxid_acc_prot)                                                         # array enthält prot und nucl .dmp Dateien
blast=/home/Software/ncbi-blast-2.4.0+/bin/                                                                                            # orginal Blastpfad
taxvz=/home/blast_db/taxonomy                                                                             # orginal Taxonomy-db pfad
data=TID-10239_Viridae_NuclGis.txt                                                                        # TID-10239_Viridae_Gis.txt
dat=TID-10239_Viridae.fna                                                                                 # TID-10239_Viridae.fna
virendmp=TID-10239_Viridae_gi-taxid-nucl.dmp                                                                
  # neue .dmp-Datei mit den Gi-Taxids nur für die Viren (ermöglicht beim Scannen schnelleren Zugriff)

# erste Vorbereitungen
if [ -d $arbeitsvz ] ; then                                                                               # |
  echo "There is already a Folder named '$arbeitsvz' in your current output directory."                   # | Abbruch des Programms, wenn der
  echo -e "Please delete or rename this Folder and start the analysis again.\n"                           # | Ordner bereits existiert.
  exit                                                                                                    # V
fi
for a in ${array[*]} ; do mkdir -p ${arbeitsvz}/$a ; done   																							# -p unterdrückt Fehlermeldungen
cd $arbeitsvz                                                                                             # Wechsel ins aktuelle Vz
for a in ${array[*]} ; do ln -s $taxvz/$a.dmp $a/$a.dmp 2>/dev/null ; done                                # Link von der Orginal .dmp-Dateien erzeugen

########################################################################### Funktionen ###################################################################

# Funktion zum Splitten der Orginal-.dmp-Datei in 100 Stücke zum effizienteren Durchsuchen
function OrgDmpDataSplitting()                                                                            # Funktionsname OrgDmpDataSplitting
{ grep "^$line" $a.dmp > $a/$a-${line}.dmp ;}                                                          # hole alle 12...... > gi_taxid_nucl-12.dmp
                                                         
# Schnellere Grep-Funktion zur Zuordnung der viralen-Taxids aus den gesplitteten .dmp-Dateien über Case  
function GiExtraction()                                                                                       
{
  case "${line}" in                                                                                       # line: jeweilige Viridae-Gi
    1*) case "${line}" in                                                                                 # beginnt diese mit der 1, schau nach der
          1)    echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-1.dmp`\n" >> tmp.dmp ;;
          10*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-10.dmp`\n" >> tmp.dmp ;;          # 2. Ziffer und geh dann in die 
          11*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-11.dmp`\n" >> tmp.dmp ;;          # entsprechende Datei und grep sie dort
          12*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-12.dmp`\n" >> tmp.dmp ;;
          13*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-13.dmp`\n" >> tmp.dmp ;;
          14*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-14.dmp`\n" >> tmp.dmp ;;
          15*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-15.dmp`\n" >> tmp.dmp ;;
          16*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-16.dmp`\n" >> tmp.dmp ;;
          17*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-17.dmp`\n" >> tmp.dmp ;;
          18*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-18.dmp`\n" >> tmp.dmp ;;
          19*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-19.dmp`\n" >> tmp.dmp ;;
        esac ;;
    2*) case "${line}" in
          2)    echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-2.dmp`\n" >> tmp.dmp ;;
          20*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-20.dmp`\n" >> tmp.dmp ;;           
          21*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-21.dmp`\n" >> tmp.dmp ;;          
          22*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-22.dmp`\n" >> tmp.dmp ;;
          23*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-23.dmp`\n" >> tmp.dmp ;;
          24*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-24.dmp`\n" >> tmp.dmp ;;
          25*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-25.dmp`\n" >> tmp.dmp ;;
          26*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-26.dmp`\n" >> tmp.dmp ;;
          27*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-27.dmp`\n" >> tmp.dmp ;;
          28*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-28.dmp`\n" >> tmp.dmp ;;
          29*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-29.dmp`\n" >> tmp.dmp ;;
        esac ;;
    3*) case "${line}" in
          3)    echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-3.dmp`\n" >> tmp.dmp ;;
          30*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-30.dmp`\n" >> tmp.dmp ;;           
          31*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-31.dmp`\n" >> tmp.dmp ;;          
          32*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-32.dmp`\n" >> tmp.dmp ;;
          33*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-33.dmp`\n" >> tmp.dmp ;;
          34*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-34.dmp`\n" >> tmp.dmp ;;
          35*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-35.dmp`\n" >> tmp.dmp ;;
          36*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-36.dmp`\n" >> tmp.dmp ;;
          37*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-37.dmp`\n" >> tmp.dmp ;;
          38*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-38.dmp`\n" >> tmp.dmp ;;
          39*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-39.dmp`\n" >> tmp.dmp ;;
        esac ;;
    4*) case "${line}" in
          4)    echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-4.dmp`\n" >> tmp.dmp ;;
          40*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-40.dmp`\n" >> tmp.dmp ;;           
          41*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-41.dmp`\n" >> tmp.dmp ;;          
          42*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-42.dmp`\n" >> tmp.dmp ;;
          43*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-43.dmp`\n" >> tmp.dmp ;;
          44*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-44.dmp`\n" >> tmp.dmp ;;
          45*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-45.dmp`\n" >> tmp.dmp ;;
          46*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-46.dmp`\n" >> tmp.dmp ;;
          47*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-47.dmp`\n" >> tmp.dmp ;;
          48*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-48.dmp`\n" >> tmp.dmp ;;
          49*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-49.dmp`\n" >> tmp.dmp ;;
        esac ;;
    5*) case "${line}" in
          5)    echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-5.dmp`\n" >> tmp.dmp ;;
          50*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-50.dmp`\n" >> tmp.dmp ;;           
          51*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-51.dmp`\n" >> tmp.dmp ;;          
          52*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-52.dmp`\n" >> tmp.dmp ;;
          53*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-53.dmp`\n" >> tmp.dmp ;;
          54*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-54.dmp`\n" >> tmp.dmp ;;
          55*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-55.dmp`\n" >> tmp.dmp ;;
          56*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-56.dmp`\n" >> tmp.dmp ;;
          57*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-57.dmp`\n" >> tmp.dmp ;;
          58*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-58.dmp`\n" >> tmp.dmp ;;
          59*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-59.dmp`\n" >> tmp.dmp ;;
        esac ;;
    6*) case "${line}" in
          6)    echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-6.dmp`\n" >> tmp.dmp ;;
          60*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-60.dmp`\n" >> tmp.dmp ;;           
          61*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-61.dmp`\n" >> tmp.dmp ;;          
          62*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-62.dmp`\n" >> tmp.dmp ;;
          63*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-63.dmp`\n" >> tmp.dmp ;;
          64*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-64.dmp`\n" >> tmp.dmp ;;
          65*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-65.dmp`\n" >> tmp.dmp ;;
          66*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-66.dmp`\n" >> tmp.dmp ;;
          67*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-67.dmp`\n" >> tmp.dmp ;;
          68*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-68.dmp`\n" >> tmp.dmp ;;
          69*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-69.dmp`\n" >> tmp.dmp ;;
        esac ;;
    7*) case "${line}" in
          7)    echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-7.dmp`\n" >> tmp.dmp ;;
          70*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-70.dmp`\n" >> tmp.dmp ;;           
          71*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-71.dmp`\n" >> tmp.dmp ;;          
          72*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-72.dmp`\n" >> tmp.dmp ;;
          73*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-73.dmp`\n" >> tmp.dmp ;;
          74*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-74.dmp`\n" >> tmp.dmp ;;
          75*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-75.dmp`\n" >> tmp.dmp ;;
          76*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-76.dmp`\n" >> tmp.dmp ;;
          77*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-77.dmp`\n" >> tmp.dmp ;;
          78*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-78.dmp`\n" >> tmp.dmp ;;
          79*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-79.dmp`\n" >> tmp.dmp ;;
        esac ;;
    8*) case "${line}" in
          8)    echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-8.dmp`\n" >> tmp.dmp ;;
          80*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-80.dmp`\n" >> tmp.dmp ;;           
          81*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-81.dmp`\n" >> tmp.dmp ;;          
          82*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-82.dmp`\n" >> tmp.dmp ;;
          83*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-83.dmp`\n" >> tmp.dmp ;;
          84*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-84.dmp`\n" >> tmp.dmp ;;
          85*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-85.dmp`\n" >> tmp.dmp ;;
          86*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-86.dmp`\n" >> tmp.dmp ;;
          87*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-87.dmp`\n" >> tmp.dmp ;;
          88*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-88.dmp`\n" >> tmp.dmp ;;
          89*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-89.dmp`\n" >> tmp.dmp ;;
        esac ;;
    9*) case "${line}" in
          9)    echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-9.dmp`\n" >> tmp.dmp ;;
          90*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-90.dmp`\n" >> tmp.dmp ;;           
          91*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-91.dmp`\n" >> tmp.dmp ;;          
          92*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-92.dmp`\n" >> tmp.dmp ;;
          93*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-93.dmp`\n" >> tmp.dmp ;;
          94*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-94.dmp`\n" >> tmp.dmp ;;
          95*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-95.dmp`\n" >> tmp.dmp ;;
          96*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-96.dmp`\n" >> tmp.dmp ;;
          97*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-97.dmp`\n" >> tmp.dmp ;;
          98*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-98.dmp`\n" >> tmp.dmp ;;
          99*)  echo -e "`grep "^\<$line\>" taxid_gi_nucl/taxid_gi_nucl-99.dmp`\n" >> tmp.dmp ;;
        esac ;;
  esac
}

# schnelle Funktion zur Extraktion der Sequencen aus der lokalen NCBInt-Datenbank
function GiDownload()                                                                                     # Funktionsname GiDownload
{ more $i | $blast/blastdbcmd -db /home/blast_db/ncbi/nt -entry_batch - -out ${i}.fna 2>/dev/null ;}      # Blast Sequenzextraktionsbefehl
  # jede i (part-xy) mit den darin enthaltenen GIs wird an die NCBI db übergeben und als Fasta Datei wieder ausgespuckt
  
####################################################################### Funktionsaufrufe #################################################################

# Vertauschen der Spalten GI und TaxID der Orginal gi_taxid_nucl.dmp Datei
echo -n "Writing gi_taxid_nucl.dmp to taxid_gi_nucl.dmp... "                                # Commandline Info
#awk '{print $1}' $taxvz/gi_taxid_nucl.dmp > gi_nucl.dmp                               # | GI- und TaxID Spalte
#echo -ne "\nCutting TaxID column from gi_taxid_nucl.dmp... "                          # Commandline Info
#awk '{print $2}' $taxvz/gi_taxid_nucl.dmp > taxid_nucl.dmp                            # | werden vertauscht und paste
#echo -ne "\nCombining columns to taxid_gi_nucl.dmp... "                               # Commandline Info
#paste taxid_nucl.dmp gi_nucl.dmp > taxid_gi_nucl.dmp                                  # V in taxid_gi_nucl.dmp 
#rm taxid_gi_nucl/taxid_gi_nucl.dmp ; rm gi_nucl.dmp ; rm taxid_nucl.dmp               # Entfernen der Einzelspalten
awk '{x="taxid_gi_nucl.dmp"} {print $2,$1 > x}' < $taxvz/gi_taxid_nucl.dmp
ln -s $arbeitsvz/taxid_gi_nucl.dmp $arbeitsvz/taxid_gi_nucl/taxid_gi_nucl.dmp         # neuen Link setzen 
cp ${taxvz}/gi_taxid_nucl.dmp ${arbeitsvz}/.
cp ${taxvz}/gi_taxid_prot.dmp ${arbeitsvz}/.
cp ${taxvz}/nucl_gb.accession2taxid  ${arbeitsvz}/.
cp ${taxvz}/prot.accession2taxid ${arbeitsvz}/.
cut -d "|" -f1 /home/blast_db/taxonomy/nodes.dmp | tr -d [:blank:] > childTaxid.tmp             # extract the first column holding the child-TAX-ID from nodes.dmp and remove any whitespace
cut -d "|" -f2 /home/blast_db/taxonomy/nodes.dmp | tr -d [:blank:] > parentTaxid.tmp            # extract the second column holding the parent-TAX-ID from nodes.dmp and remove any whitespace
awk 'BEGIN {FS="\t"} {x="taxid_acc_nucl.dmp"} {print $3 FS $1 > x}' < nucl_gb.accession2taxid                      # faster awk alternative to cut twice and paste subsequently
awk 'BEGIN {FS="\t"} {x="taxid_acc_prot.dmp"} {print $3 FS $1 > x}' < prot.accession2taxid                      # faster awk alternative to cut twice and paste subsequently

rm nucl_gb.accession2taxid; rm prot.accession2taxid
paste childTaxid.tmp parentTaxid.tmp | sort > childParentTaxid.dmp                              # fuse child-TAX-ID and parent-TAX-ID columns into a single file and sort ascending by child-TAX-ID
paste parentTaxid.tmp childTaxid.tmp | sort > parentChildTaxid.dmp                              # fuse parent-TAX-ID and child-TAX-ID columns into a single file and sort ascending by parent-TAX-ID
rm parentTaxid.tmp ; rm childTaxid.tmp


# Start der Splitfunktion der Orginal-.dmp -Datei für nucl und prot in 99 Einzeldateien zum schnelleren Zugriff
for a in ${array[*]} ; do                                                             # array enthält prot - und nucl .dmp Dateien
  echo -ne "\nSplitting of $a.dmp file... "                                           # Splitinfo zur entspr. .dmp Datei
  i=10                                                                                # Initialisierung i=10
  until (( $i == 100 )) ; do echo $i >> zahlen.txt ; ((i=$i+1)) ; done                # Zahlen von 10-99 werden in zahlen.txt geschrieben 
  while read line ; do OrgDmpDataSplitting $line & done < zahlen.txt ; wait         # Funktionsaufruf für Splitten der .dmp Datei in 89 Einzelteile
  i=1                                                                                 # Initialisierung i=1 für 1-9
  until (( $i == 10 )) ; do                                                          # im Prinzip äquivalent zur Splitfunktion, außer dass nur 
    grep "^\<$i\>" $a.dmp > $a/$a-$i.dmp &                                         # die GIs 1-9 verwendet werden, ohne darauffolgende Ziffern
    ((i=$i+1))                                                                        # nach jeder Gi Extraktion wird zahl um 1 erhöht  
  done ; wait ; rm zahlen.txt                                                        # Verteilung auf alle Kerne (max 9 in diesem Fall)
done                                                                                  # und fertig

rm childParentTaxid.dmp ; rm parentChildTaxid.dmp ; rm gi_taxid_nucl.dmp ; rm taxid_gi_nucl.dmp ; rm gi_taxid_prot.dmp
rm taxid_acc_nucl.dmp ; rm taxid_acc_prot.dmp

# Ermittlung aller viralen GIs
echo -ne "\nViral GI extraction (ca 20min)... "                                           # User Info
awk '$9~/9/ {print $1}' $taxvz/nodes.dmp > VirenTaxIDs-1.txt                              # wenn in Spalte 9 die ) steht, printe 1. Spalte der Datei
awk '$9~/9/ {print $3}' $taxvz/nodes.dmp > VirenTaxIDs-2.txt                              # wenn in Spalte 9 die ) steht, printe 3. Spalte der Datei
cat VirenTaxIDs-*.txt | sort | uniq > VirenTaxIDs-su.txt                                  # beide TaxID Spalten zusammenschreiben, ordnen und uniquen
(while read line ; do GiExtraction $line & done < VirenTaxIDs-su.txt ; wait)              # Funktionsaufruf
echo -ne "\nCutting viral GI column... " ; cut -f1 tmp.dmp | sed '/^$/d' > tmp_taxid.dmp  # |
echo -ne "\nCutting viral TaxID column... " ; cut -f2 tmp.dmp | sed '/^$/d' > tmp_gi.dmp  # | Gi & TaxID Spalten werden wieder vertauscht
echo -ne "\nCombining columns... " ; paste tmp_gi.dmp tmp_taxid.dmp > $virendmp           # V
rm VirenTaxIDs-*.txt ; rm tmp*.dmp                                                        # Entfernen der temporären TaxID Dateien 

# Extraktion aller viralen Sequenzen aus der lokalen NCBI Datenbank und Erstellen einer entsprechenden Blastdatenbank
echo -ne "\nSequence isolation... "                                                   # Startinfo
((p=`wc -l < $virendmp`/96))                                                          # Ermittlung der Zeilenanzahl für die Splitdateien
cut -f1 -d " " $virendmp | split -l $p - part-                                               # Split der .dmp Datei (nur Gi Spalte)
(for i in part-* ; do GiDownload $i & done ; wait)                                    # Start der Downloadfunktion
cat part-*.fna > $dat                                                                 # alle fasta Dateien werden in eine große überführt
rm part-*                                                                             # und die Teilstücke entfernt
echo -ne "\nCreate Blast database... "                                                # Info über Prozessbeendigung und die nächste Aktion
$blast/makeblastdb -dbtype nucl -in $dat -out ${dat}-DB 1>/dev/null                   # Erstellen einer Blastdatenbank aus der Vieren Fastadatei
rm TID-10239_Viridae.fna                  # Entfernen aller überflüssigen Dateien

if [ -d /home/blast_db/taxLok ]
    then
        mv /home/blast_db/taxLok /home/blast_db/taxLokBak_${zeit}
        cp -r ${arbeitsvz}/* /home/blast_db/taxLok
    else
        mkdir /home/blast_db/taxLok
        cp -r ${arbeitsvz}/* /home/blast_db/taxLok
fi

echo -e "finished!\n"                                                                 # & that's it!



 