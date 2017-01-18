#!/bin/bash


anfang=`date`   # Startzeit ermitteln und speichern
zeit=`echo ${anfang} | tr -d " "`     # Datum fuer Verzeichnisnamen erzeugen
HOST='ftp.ncbi.nih.gov'   # ftp-Server fuer den Download der Taxonomie-Datenbank spezifizieren
USER='anonymous'          # Benutzername fuer die Anmeldung am ftp-Server
PASSWD='your.email@server.xy'   # Kennwort fuer die Anmeldung am ftp-Server
dbs="nt nr taxdb"           # durch update_blastdb.pl zu aktualisierende Datenbanken; Auswahl: 16SMicrobial env_nr env_nt est est_human est_mouse est_others gss htgs human_genomic human_genomic_transcript mouse_genomic_transcript nr nt other_genomic pataa patnt pdbaa pdbnt refseq_genomic refseq_protein refseq_rna refseqgene sts swissprot taxdb tsa_nt vector wgs
#empfaenger="your.email@server.xy,your.email@server.xy,your.email@server.xy" 


mkdir ~/.tmp/ncbidb_${zeit}    # Verzeichnis fuer Download anlegen
cd ~/.tmp/ncbidb_${zeit}    # ins angelegte Verzeichnis wechseln
echo "ncbi" > md5ncbi        # Datei fuer den Vergleich der Dateien anlegen 
echo "download" > md5download      # Datei fuer den Vergleich der Dateien anlegen
while [[ `diff md5ncbi md5download` != "" ]]      # Download solange wiederholen, bis die Listen der md5-Summen identisch sind 
do 
  update_blastdb.pl ${dbs}       # Download der Datenbanken
  md5sum *.tar.gz | cut -d " " -f 1 | sort > md5download     # nach Abschluss des Downloads md5-Summen der heruntergeladenen Dateien berechnen und in Datei schreiben
  cat *md5 | cut -d " " -f 1 | sort > md5ncbi    # md5-Summen vom NCBI in eine Datei schreiben
done

for datei in *.tar.gz     # Schleife zum konsekutiven Entpacken der heruntergeladenen Archiv-Dateien
do 
    tar -zxvpf ${datei}      # Entpacken der Dateien (-f) ${datei} mittels zip (-z) unter Beibehaltung der Rechte (-p)
done

mkdir ~/.tmp/taxdb_${zeit}      # neues Verzeichnis fuer den Download der Taxonomie-Datenbank anlegen
cd ~/.tmp/taxdb_${zeit}        # ins Taxonomie-Verzeichnis wechseln
### im folgenden Teil werden die notwendigen Dateien vom NCBI-ftp-Server heruntergeladen; die Schritte sind im Einzelnen:
# 0. ftp-Aufruf mit den Parametern i (keine Abfragen bei multiplen Downloads), n (keine Abfrage der Authentifizierung durch den Server) und v (Verbose)
# 1. Verbindung zum ftp-Server herstellen
# 2. Anmeldung am ftp-Server mit Benutzername $USER und Passwort $PASSWD
# 3. Wechsel ins Quellverzeichnis
# 4. Download aller Dateien mit gz im Namen aus dem Quellverzeichnis
# 5. Schliessen der ftp-Verbindung 
ftp -inv ${HOST} <<END_SCRIPT
user $USER $PASSWD       
cd pub/taxonomy        
mget *.gz*  
quit     
END_SCRIPT


ftp -inv ${HOST} <<END_SCRIPT
user $USER $PASSWD       
cd pub/taxonomy/accession2taxid   
mget *.gz*  
quit     
END_SCRIPT

rm dead_*
rm nucl_gss*
rm nucl_est*
rm nucl_wgs*

for datei in *.tar.gz     # Schleife zum konsekutiven Entpacken der heruntergeladenen Archiv-Dateien
do 
    tar -zxvpf ${datei}      # Entpacken der Dateien (-f) ${datei} mittels zip (-z) unter Beibehaltung der Rechte (-p)
done
for datei in *gz     # Schleife zum konsekutiven Entpacken der heruntergeladenen Archiv-Dateien
do 
    gunzip -v ${datei}      # Entpacken der gz-Dateien
done


while [[ `pgrep "blast"` != "" ]]      # Pruefung auf laufende BLAST-Prozesse; falls ein BLAST laeuft, dann 5 min bis zur naechsten Pruefung warten
do
  echo date
  echo "verzoegert"
  sleep 5m
done





mv -f /home/blast_db/ncbi /home/blast_db/ncbi_bak_${zeit}       # nach Abschluss laufender BLAST-Prozesse alte DB in Backup verschieben
mv -f ~/.tmp/ncbidb_${zeit} /home/blast_db/ncbi        # neue DB in Ordner /home/blast_db/ncbi verschieben
mv -f /home/blast_db/taxonomy /home/blast_db/taxonomy_bak_${zeit}    # alte Taxonomie-Dateien in Backup verschieben
mv -f ~/.tmp/taxdb_${zeit} /home/blast_db/taxonomy     # neue Taxonomie-Dateien in /home/blast_db/Taxonomy verschieben



~/skript/SplitDmp_VirGisTaxEx_Db1.sh


schluss=`date`
#echo Die am ${anfang} gestartete Aktualisierung der Datenbanken ${dbs} vom NCBI wurde am ${schluss} abgeschlossen. Die zuvor vorhandenen Datenbank-Versionen wurden unter /home/blast_db/ncbi_bak_${zeit} und /home/blast_db/taxonomy_bak_${zeit} gesichert. > mailtext.txt 
# mail -s "Update of NCBI databases finished" your.email@server.xy < mailtext.txt       # Email-Benachrichtigung schicken -c ${empfaenger}

