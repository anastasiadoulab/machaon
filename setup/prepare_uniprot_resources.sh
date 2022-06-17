#!/bin/bash

[ ! -d refseq ] && mkdir refseq

echo 'Constructing RefSeq NM id mappings...'
zgrep -iP 'refseq_nt\tnm' idmapping.dat.gz > refseq_nm_map.tsv
echo 'Constructing RefSeq XM id mappings...'
zgrep -iP 'refseq_nt\txm' idmapping.dat.gz > refseq_xm_map.tsv
mv refseq_nm_map.tsv refseq/refseq_nm_map.tsv
mv refseq_xm_map.tsv refseq/refseq_xm_map.tsv
echo 'Id mappings are ready.'