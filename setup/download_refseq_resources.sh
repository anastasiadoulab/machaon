#!/bin/bash

retrieve_refseq_data() {
    if [ -f $1 ] && gzip -t $1; then
      echo "$1 : This file is already retrieved."
    else
      echo "Downloading : $1"
      ftp -in "ftp.ncbi.nlm.nih.gov" <<SCRIPT_END
    user anonymous -
    binary
    cd $2
    mget ${line}
SCRIPT_END
      if ! gzip -t $1; then
        echo "$1 : The file seems to be corrupt due to an erroneous transfer with the server. Please run this script again a while later to re-download any malformed files."
      fi
      sleep 10
    fi
}

echo "Downloading viral genomes..."
[ ! -d vir_fna ] && mkdir vir_fna
cd vir_fna
curl -v --silent ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/ 2>&1 | grep -o "viral.*.genomic.fna.gz" | while read -r line ; do
  retrieve_refseq_data ${line} "refseq/release/viral"
done
cd ..

echo "Downloading viral genome annotations..."
[ ! -d vir_gbff ] && mkdir vir_gbff
cd vir_gbff
curl -v --silent ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/ 2>&1 | grep -o "viral.*.genomic.gbff.gz" | while read -r line ; do
  retrieve_refseq_data ${line} "refseq/release/viral"
done
cd ..

echo "Downloading human transcripts..."
[ ! -d rna_fna ] && mkdir rna_fna
cd rna_fna
curl -v --silent ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/ 2>&1 | grep -o "human..*.rna.fna.gz" | while read -r line ; do
  retrieve_refseq_data ${line} "refseq/H_sapiens/mRNA_Prot"
done
cd ..


echo "Downloading human transcripts annotations..."
[ ! -d rna_gbff ] && mkdir rna_gbff
cd rna_gbff
curl -v --silent ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/ 2>&1 | grep -o "human..*.rna.gbff.gz" | while read -r line ; do
  retrieve_refseq_data ${line} "refseq/H_sapiens/mRNA_Prot"
done
cd ..
