#!/bin/bash

[ ! -d refseq ] && mkdir refseq

echo "Process viral genome data..."
cd vir_fna
gunzip *.gz
cat *.fna > viral.genomic.fna
mv viral.genomic.fna ../viral.genomic.fna
cd ..
rm -rf vir_fna
mv viral.genomic.fna refseq/viral.genomic.fna


echo "Process viral genome annotation data..."
cd vir_gbff
for f in *.gz; do
  t=${f%.gz}.chunk
  zgrep "ACCESSION\|CDS\(\s\)\{3\}\|5'UTR\(\s\)\{3\}\|3'UTR\(\s\)\{3\}\|GeneID:\|//$" $f| awk '$0 != prev{print; prev=$0;}' | awk '{$1=$1;print}' | tr "'" "-" > $t
done
cat *.chunk > virinfo.gbff
mv virinfo.gbff ../virinfo.gbff
cd ..
rm -rf vir_gbff
mv virinfo.gbff refseq/virinfo.gbff


echo "Process human transcript data..."
cd rna_fna
gunzip *.gz
cat *.fna > human.rna.fna
mv human.rna.fna ../human.rna.fna
cd ..
rm -rf rna_fna
mv human.rna.fna refseq/human.rna.fna


echo "Process human transcript annotation data..."
cd rna_gbff
for f in *.gz; do
  t=${f%.gz}.chunk
  zgrep "ACCESSION\|CDS" $f > $t
done
cat *.chunk > cdsinfo.gbff
mv cdsinfo.gbff ../cdsinfo.gbff
cd ..
rm -rf rna_gbff
mv cdsinfo.gbff refseq/cdsinfo.gbff

echo "RefSeq local resources are ready."
