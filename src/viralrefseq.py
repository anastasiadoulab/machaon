import os
import subprocess
import traceback
import re
from collections import defaultdict

from Bio.Seq import Seq


class ViralRefSeq:

    def __init__(self, refseq_id, gene_id):
        self.refseq_id = refseq_id
        self.gene_id = gene_id
        self.root_disk = ''
        self.refseq_dataset_path = ''
        self.refseq_annotation_path = ''
        self.sequence = ''
        self.genome_info = []

    def set_root_disk(self, root_disk):
        self.root_disk = root_disk
        self.refseq_dataset_path = os.path.join(self.root_disk, 'refseq', 'viral.genomic.fna')
        self.refseq_annotation_path = os.path.join(self.root_disk, 'refseq', 'virinfo.gbff')

    def fetch_viral_genome(self):
        # Fetch viral genome from RefSeq resources
        genomic_sequence = ''
        try:
            output = subprocess.run(' '.join(
                ['timeout', '80s', 'awk', '\'BEGIN{RS=">";FS="\\n"}NR>1{if', ''.join(['($1~/', self.refseq_id, '/)']),
                 'print ">"$0}\'',
                 ''.join(['\'', self.refseq_dataset_path, '\''])]), shell=True, capture_output=True)
            parts = output.stdout.decode("utf-8").split('\n')
            if (len(parts) > 1):
                genomic_sequence = ''.join(parts[1:])
        except:
            print(''.join([self.refseq_id, ' not found in refseq set']))
            print(traceback.format_exc())
        self.sequence = genomic_sequence

    def fetch_viral_genome_info(self):
        # Fetch viral genome annotations from RefSeq resources
        lines = []
        try:
            pattern = ''.join(['/', self.refseq_id, '/,/\/\//p'])
            output = subprocess.run(['timeout', '30s', 'sed', '-n', pattern, self.refseq_annotation_path], capture_output=True)
            lines = output.stdout.decode("utf-8").split('\n')
        except:
            print(''.join([self.refseq_id, ' not found in refseq annotations']))
            print(traceback.format_exc())
        self.genome_info = lines

    def parse_gene_data(self):
        # Parse all available data for a gene (sequence & annotations)
        data = defaultdict()
        if ('5-UTR' in self.genome_info[1]):
            data['5UTR'] = self.retrieve_sequence(self.genome_info[1].split()[1])
            if ('3-UTR' in self.genome_info[-3]):
                data['3UTR'] = self.retrieve_sequence(self.genome_info[-3].split()[1])
        parsed = self.fetch_gene_info(self.genome_info)
        data = {**parsed, **data}
        return data

    def retrieve_sequence(self, gene_info):
        # Retrieve gene sequence joining all its parts
        # according to the available annotations
        gene_info = gene_info.replace(')', '').split('(')
        complement = False
        multiple = False
        sequence = []
        ranges = []

        for part_index, part in enumerate(gene_info):
            if ('complement' in part):
                complement = True
            elif ('join' in part):
                multiple = True
                ranges = gene_info[part_index + 1].split(',')
            elif (multiple is False):
                ranges.append(part.replace('>', '').replace('<', ''))

        for segment in ranges:
            segment = segment.split('..')
            if (len(segment) > 1):
                sequence.append(self.sequence[int(segment[0]) - 1:int(segment[1])])
            else:
                sequence.append(self.sequence[int(segment[0])])

        if (complement is True):
            final_sequence = []
            for segment in sequence:
                sequence = Seq(segment)
                final_sequence.append(str(sequence.reverse_complement()))
            sequence = ''.join(final_sequence)
        else:
            sequence = ''.join(sequence)

        return sequence

    def fetch_gene_info(self, genome_info):
        # Find the annotations of a specific gene
        parsed = defaultdict()
        if (self.gene_id is not None):
            pattern = ''.join(['/db_xref="GeneID:', self.gene_id, '";'])
            result = re.search(''.join([pattern, '(.*)', pattern]), ';'.join(genome_info))
            if (result is not None):
                result = result.group(1)
                result = [x for x in result.split(';') if len(x) > 0]
                for entry in result:
                    if ('db_xref' not in entry):
                        parts = entry.split()
                        parsed[parts[0].replace('-', '')] = self.retrieve_sequence(parts[1])
        return parsed
