import os.path
from collections import defaultdict
from Bio import pairwise2, BiopythonDeprecationWarning
from Bio.pairwise2 import format_alignment
import traceback
import subprocess
from DSSPparser import parseDSSP
from src.pdbhandler import PDBHandler
import re
import pandas as pd
import numpy as np
from parasail import matrix_create
import parasail
import warnings
warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)
from Bio.SubsMat import MatrixInfo as matlist
from Bio import Align


class StructureAligner:

    def __init__(self, backend='biopython'):
        self.root_disk = ''
        self.verbose = False
        self.structure_data_dir = ''
        self.dssp_cache_dir = ''
        self.stride_cache_dir = ''
        self.pdb_dataset_path = ''
        self.dssp_executable_path = './mkdssp'
        self.stride_executable_path = './stride'
        self._max_sequence_size = 10000000
        self._dssp_notations = ['H', 'B', 'E', 'G', 'I', 'P', 'T', 'S', '.']

        # Group residues by their Hydropathy Index
        self._hydropathy = defaultdict()
        self._hydropathy['I'] = 0
        self._hydropathy['V'] = 0
        self._hydropathy['L'] = 1
        self._hydropathy['F'] = 2
        self._hydropathy['C'] = 2
        self._hydropathy['M'] = 2
        self._hydropathy['A'] = 2
        self._hydropathy['G'] = 3
        self._hydropathy['T'] = 3
        self._hydropathy['Y'] = 3
        self._hydropathy['W'] = 3
        self._hydropathy['P'] = 3
        self._hydropathy['S'] = 3
        self._hydropathy['H'] = 4
        self._hydropathy['N'] = 4
        self._hydropathy['Q'] = 4
        self._hydropathy['D'] = 4
        self._hydropathy['E'] = 4
        self._hydropathy['K'] = 4
        self._hydropathy['R'] = 5

        # 3-to-1 aminoacid letter code dictionary
        self._lettercode = defaultdict()
        self._lettercode['ILE'] = 'I'
        self._lettercode['VAL'] = 'V'
        self._lettercode['LEU'] = 'L'
        self._lettercode['PHE'] = 'F'
        self._lettercode['CYS'] = 'C'
        self._lettercode['CYX'] = 'C'
        self._lettercode['CYM'] = 'C'
        self._lettercode['MET'] = 'M'
        self._lettercode['ALA'] = 'A'
        self._lettercode['GLY'] = 'G'
        self._lettercode['THR'] = 'T'
        self._lettercode['TYR'] = 'Y'
        self._lettercode['TYM'] = 'Y'
        self._lettercode['TRP'] = 'W'
        self._lettercode['PRO'] = 'P'
        self._lettercode['SER'] = 'S'
        self._lettercode['HIS'] = 'H'
        self._lettercode['HIP'] = 'H'
        self._lettercode['HIE'] = 'H'
        self._lettercode['HID'] = 'H'
        self._lettercode['ASN'] = 'N'
        self._lettercode['GLN'] = 'Q'
        self._lettercode['ASP'] = 'D'
        self._lettercode['ASH'] = 'D'
        self._lettercode['GLU'] = 'E'
        self._lettercode['GLH'] = 'E'
        self._lettercode['LYS'] = 'K'
        self._lettercode['LYN'] = 'K'
        self._lettercode['ARG'] = 'R'
        self._lettercode['ARN'] = 'R'

        # Create mixed labels
        self._label_mix = []
        for notation in range(len(self._dssp_notations)):
            for hydropathy_index in range(7):
                self._label_mix.append((notation, hydropathy_index))

        self.backend = backend if backend != '' else 'biopython'
        if(backend == 'biopython'):
            self.perform_global_alignment = self.perform_global_alignment_b
            self.perform_local_alignment = self.perform_local_alignment_b
            self.get_reference_sequence_alignment = self.get_reference_sequence_alignment_b
        elif(backend == 'biopython_new'):
            self.perform_global_alignment = self.perform_global_alignment_bn
            self.perform_local_alignment = self.perform_local_alignment_bn
            self.get_reference_sequence_alignment = self.get_reference_sequence_alignment_bn
            self.get_alignment_info = self.get_alignment_info_bn
        elif (backend == 'parasail'):
            self.perform_global_alignment = self.perform_global_alignment_p
            self.perform_local_alignment = self.perform_local_alignment_p
            self.get_reference_sequence_alignment = self.get_reference_sequence_alignment_p
        else:
            raise Exception('Unknown backend.')

    def set_root_disk(self, root_disk):
        self.root_disk = root_disk
        self.structure_data_dir = os.path.sep.join([self.root_disk, 'prot_sec'])
        self.dssp_cache_dir = os.path.sep.join([self.root_disk, 'dssp_cache'])
        self.stride_cache_dir = os.path.sep.join([self.root_disk, 'stride_cache'])
        if (os.path.exists(self.structure_data_dir) is False):
            os.makedirs(self.structure_data_dir)
        if (os.path.exists(self.dssp_cache_dir) is False):
            os.makedirs(self.dssp_cache_dir)
        if (os.path.exists(self.stride_cache_dir) is False):
            os.makedirs(self.stride_cache_dir)

    def get_structure(self, structure_info, available_residues):
        pdb_id, chain_id = structure_info
        data_path = os.path.sep.join([self.structure_data_dir, ''.join([pdb_id, '_', chain_id, '.seq'])])

        # Read previously cached sequences
        if (os.path.exists(data_path)):
            with open(data_path, 'r', encoding='utf-8') as data_file:
                secondary_structure = data_file.readline().rstrip()
                primary_structure = data_file.readline().rstrip()
        else:
            # Extract 1D, 2D sequences from PDB
            secondary_structure, primary_structure = self.get_pdb_sequences(pdb_id, chain_id, available_residues)
            if (secondary_structure != ''):
                with open(data_path, 'w', encoding='utf-8') as data_file:
                    data_file.write('\n'.join([secondary_structure, primary_structure]))

        # If there was a failure on extracting 2D, extract only 1D
        if (secondary_structure == ''):
            primary_structure = self.get_pdb_primary_sequence(pdb_id, chain_id)

        if (secondary_structure == '' and self.verbose is True):
            print(''.join(['N/A 2D for ', pdb_id, '_', chain_id]))

        return secondary_structure, primary_structure

    # Extract 1D sequence from a PDB chain
    def get_pdb_primary_sequence(self, pdb_id, chain_id):
        sequence = []
        pdbhandler = PDBHandler(pdb_id)
        pdbhandler.root_disk = self.root_disk
        full_pdb_file_path = os.path.sep.join([self.pdb_dataset_path, ''.join([pdb_id, '.pdb'])])
        residue_names = pdbhandler.get_residues_names(full_pdb_file_path, chain_id)
        for name in residue_names:
            if(name in self._lettercode):
                sequence.append(self._lettercode[name])
            else:
                sequence.append('-')
        return ''.join(sequence)

    # Create a sequence based on hydropathy group label mappings
    def create_hydrophobicity_sequence(self, primary_sequence):
        labels = []
        for label_index, label in enumerate(primary_sequence):
            if ('-' == label):
                labels.append('-')
                continue
            labels.append(repr(6 if primary_sequence[label_index] not in self._hydropathy else self._hydropathy[
                primary_sequence[label_index]]))
        hydrophobicity_sequence = ''.join(labels)
        return hydrophobicity_sequence

    # Create a mixed representation by combining 1D, 2D and hydrophobicity sequences
    def create_mixed_sequence(self, primary_sequence, secondary_sequence):
        mixed_sequence = []
        labels = []
        for label_index, label in enumerate(secondary_sequence):
            if ('-' == label):
                mixed_sequence.append('-')
                labels.append('-')
                continue
            new_label = (self._dssp_notations.index(label.upper()),
                         6 if primary_sequence[label_index] not in self._hydropathy else self._hydropathy[
                             primary_sequence[label_index]])
            ascii_value = 33 + self._label_mix.index(new_label)
            # backslash
            if (ascii_value == 92):
                ascii_value = 125
            # alignment gap char
            if (ascii_value == 45):
                ascii_value = 126
            mixed_sequence.append(chr(ascii_value if ascii_value != 92 else 125))
            labels.append(repr(6 if primary_sequence[label_index] not in self._hydropathy else self._hydropathy[
                primary_sequence[label_index]]))
        hydrophobicity_sequence = ''.join(labels)
        return ''.join(mixed_sequence), hydrophobicity_sequence

    def read_ref_seq_id(self, uniprot_id):
        ref_seq_transcript_id = ''
        # Look for RefSeq ID in NM Refseq entries, checking also for isoform entries
        for i in range(5):
            uniprot_dataset_path = os.path.join(self.root_disk, 'refseq', 'refseq_nm_map.tsv')
            if (i == 0):
                isoform_uniprot_id = uniprot_id
            else:
                isoform_uniprot_id = ''.join([uniprot_id, '-', repr(i)])
            try:
                output = subprocess.run(['timeout', '30s', 'grep', '-i', isoform_uniprot_id, '-m', '1',
                                         uniprot_dataset_path], capture_output=True)
                parts = output.stdout.decode("utf-8").split('\t')
                if (len(parts) > 1):
                    ref_seq_transcript_id = parts[-1].strip()
            except:
                if (self.verbose is True):
                    print(''.join([uniprot_id, ' not found in idmapping (nm)']))
                    print(traceback.format_exc())

            # Look for RefSeq ID in XM Refseq entries (predicted sequences) if no sequence is found
            if (ref_seq_transcript_id == ''):
                try:
                    uniprot_dataset_path = os.path.join(self.root_disk, 'refseq', 'refseq_xm_map.tsv')
                    output = subprocess.run(['timeout', '30s', 'grep', '-i', isoform_uniprot_id, '-m', '1',
                                             uniprot_dataset_path], capture_output=True)
                    parts = output.stdout.decode("utf-8").split('\t')
                    if (len(parts) > 1):
                        ref_seq_transcript_id = parts[-1].strip()
                except:
                    if (self.verbose is True):
                        print(''.join([uniprot_id, ' not found in idmapping (xm)']))
                        print(traceback.format_exc())
            # Separate RefSeq ID from version number
            if (ref_seq_transcript_id != ''):
                parts = ref_seq_transcript_id.split('.')
                if (len(parts) > 1):
                    ref_seq_transcript_id = parts[0]
                break

        return ref_seq_transcript_id

    def fetch_transcript_sequence(self, ref_seq_transcript_id):
        refseq_dataset_path = os.path.join(self.root_disk, 'refseq', 'human.rna.fna')
        transcript_sequence = ''
        # Retrieve the transcript sequence belonging to the given RefSeq ID
        try:
            command = ['awk', '\'BEGIN{RS=">";FS="\\n"}NR>1{if', ''.join(['($1~/', ref_seq_transcript_id, '/)']),
                 'print ">"$0}\'', ''.join(['\'', refseq_dataset_path, '\''])]
            if (self.verbose is True):
                print(' '.join(command))
            output = subprocess.run(' '.join(['timeout', '80s'] + command), shell=True, capture_output=True)
            parts = output.stdout.decode("utf-8").split('\n')
            if (len(parts) > 1):
                transcript_sequence = ''.join(parts[1:])
        except:
            if (self.verbose is True):
                print(''.join([ref_seq_transcript_id, ' not found in refseq set']))
                print(traceback.format_exc())

        return transcript_sequence

    def fetch_transcript_cds(self, ref_seq_transcript_id):
        # Extract CDS region from RefSeq annotations
        refseq_annotation_path = os.path.join(self.root_disk, 'refseq', 'cdsinfo.gbff')
        cds_range = [-1, -1]
        try:
            command = ['grep', '-A', '1', ref_seq_transcript_id, refseq_annotation_path]
            if (self.verbose is True):
                print(' '.join(command))
            output = subprocess.run(['timeout', '30s'] + command, capture_output=True)
            lines = output.stdout.decode("utf-8").split('\n')
            if (len(lines) > 1 and 'CDS' in lines[1]):
                parts = lines[1].split()
                if (len(parts) > 1):
                    parts = parts[1].split('.')
                    if (len(parts) > 2):
                        cds_range = [int(parts[0]), int(parts[2])]
        except:
            if (self.verbose is True):
                print(''.join([ref_seq_transcript_id, ' not found in refseq annotations']))
                print(traceback.format_exc())
        return cds_range

    def fetch_transcript_parts(self, ref_seq_transcript_id):
        # Fetch untranslated and coding regions of a transcript sequence
        transcript_sequence = self.fetch_transcript_sequence(ref_seq_transcript_id)
        parts = []
        if (transcript_sequence != ''):
            cds_range = self.fetch_transcript_cds(ref_seq_transcript_id)
            if (cds_range[0] > 0):
                parts = [transcript_sequence[0:cds_range[0] - 1], transcript_sequence[cds_range[0] - 1:cds_range[1]],
                         transcript_sequence[cds_range[1]:], len(transcript_sequence)]
        return parts

    def get_pdb_sequences(self, pdb_id, chain_id, available_residues):
        # Extract 1D, 2D sequences of a PDB chain from its PDB data
        dssp_data_path = os.path.sep.join([self.dssp_cache_dir, ''.join([pdb_id, '.dssp'])])
        secondary_structure = []
        primary_structure = []
        traversed_residues = []
        if (len(available_residues['fullRange']) > 1):
            try:
                if (os.path.exists(dssp_data_path) is False):
                    self.call_dssp(pdb_id)
                if (os.path.exists(dssp_data_path)):
                    primary_structure, secondary_structure = self.parse_structures_by_dssp(dssp_data_path,
                                                             chain_id, available_residues['fullRange'])

                # If DSSP failed, call STRIDE
                if (len(secondary_structure) < 2):
                    primary_structure, secondary_structure = self.get_structures_by_stride(pdb_id, chain_id,
                                                                                           available_residues)
                    if(len(secondary_structure) == 0):
                        raise Exception('Failed to calculate 2D.')
                    else:
                        # assignment to trigger assertions
                        traversed_residues = [x for x in range(available_residues['fullRange'][0],
                                                               available_residues['fullRange'][-1] + 1)]

                if (len(primary_structure) < 2):
                    primary_structure = ['']
                    secondary_structure = ['']
                    traversed_residues = []
            except:
                if (self.verbose is True):
                    print(' | '.join([pdb_id, chain_id]))
                    print(traceback.format_exc())

            total_residue_range = []
            if (len(available_residues['fullRange']) > 1):
                total_residue_range = [x for x in range(available_residues['fullRange'][0],
                                                        available_residues['fullRange'][-1] + 1)]

            # Tests for extracted sequences
            if (len(traversed_residues) > 0):
                assert len(traversed_residues) == len(total_residue_range), ' | '.join([pdb_id, chain_id])
                assert np.sum(np.subtract(total_residue_range, traversed_residues)) == 0, ' | '.join([pdb_id, chain_id])
                assert len(secondary_structure) == len(total_residue_range), ' | '.join([pdb_id, chain_id])
                assert len(secondary_structure) == len(primary_structure), ' | '.join([pdb_id, chain_id])
        return ''.join(secondary_structure), ''.join(primary_structure)

    @staticmethod
    def parse_structures_by_dssp(dssp_data_path, chain_id, available_residues):
        secondary_structure = []
        primary_structure = []
        traversed_residues = []

        # Determine 2D sequence by the PDB file
        parser = parseDSSP(dssp_data_path)
        parser.parse()
        pddict = parser.dictTodataframe()
        pddict = pddict[pddict['chain'] == chain_id]
        previous = -1

        # Process residue ids
        pddict['inscode'] = pd.to_numeric(pddict['inscode'])
        pddict = pddict.sort_values(by='inscode')
        pddict = pddict[pddict['inscode'] > 0]

        # Add loops as dots in the sequence.
        # Missing residues are treated as gaps.
        for index, row in pddict.iterrows():
            current_residue_index = int(row['inscode'])
            if (repr(int(row['inscode'])) != repr(row['inscode'])):
                continue
            if (current_residue_index > available_residues[-1]):
                break
            if (previous == -1 and current_residue_index > available_residues[0]):
                previous = available_residues[0] - 1
            if (previous != -1 and current_residue_index - previous > 1):
                for resIndex in range(previous + 1, current_residue_index):
                    traversed_residues.append(resIndex)
                    primary_structure.append('-')
                    if (resIndex in available_residues):
                        secondary_structure.append('.')
                    else:
                        secondary_structure.append('-')
            if (current_residue_index in available_residues):
                traversed_residues.append(current_residue_index)
                primary_structure.append(row['aa'].strip().upper())
                struct = row['struct'].strip()
                secondary_structure.append(struct if len(struct) > 0 else '.')
            elif (previous != -1):
                traversed_residues.append(current_residue_index)
                if (traversed_residues[-1] <= available_residues[-1]):
                    primary_structure.append('-')
                    secondary_structure.append('-')
                else:
                    traversed_residues = traversed_residues[:-1]
            if (len(traversed_residues) > 0):
                previous = current_residue_index
        if(len(traversed_residues) > 0):
            if (traversed_residues[-1] < available_residues[-1]):
                for resIndex in range(traversed_residues[-1] + 1, available_residues[-1] + 1):
                    traversed_residues.append(resIndex)
                    primary_structure.append('-')
                    secondary_structure.append('-')
        return primary_structure, secondary_structure

    @staticmethod
    def get_alignment_content(alignment):
        # Find content of an alignment
        content = defaultdict(int)

        for i, (a, b) in enumerate(zip(alignment[0], alignment[2])):
            if a == b:
                content[a] += 1

        if(len(content) > 0):
            total = np.sum(list(content.values()))
            content = {key: round((content[key] / total) * 100, 2) for key in content}
        else:
            content = '-'

        return content

    @staticmethod
    def get_alignment_info(alignment, gap_included=True):
        # Find mismatches and matches of an alignment
        matches = alignment[1].count('|')
        mismatches = alignment[1].count('.') + alignment[1].count(':')
        gaps = alignment[1].count(' ')
        if gap_included is True:
            mismatches += gaps

        assert ((matches + mismatches) if gap_included is True else (matches + mismatches + gaps)) == max(
            len(alignment[0]),
            len(alignment[2]))
        return matches, mismatches, gaps

    @staticmethod
    def get_alignment_info_bn(alignment, gap_included=True):
        # Find mismatches and matches of an alignment
        matches = alignment[1].count('|')
        mismatches = alignment[1].count('.') + alignment[1].count(':')
        gaps = alignment[1].count('-')
        if gap_included is True:
            mismatches += gaps

        assert ((matches + mismatches) if gap_included is True else (matches + mismatches + gaps)) == max(
            len(alignment[0]),
            len(alignment[2]))
        return matches, mismatches, gaps

    @staticmethod
    def get_alignment_mask(alignment):
        # Construct a bit mask based on the matches of an alignment.
        # 1 for match or 0 for no match
        mask = [1 if a == b else 0 for i, (a, b) in enumerate(zip(alignment[0], alignment[1])) if a != '-']
        return mask

    def calculate_identity(self, alignment):
        # Calculate sequence identity (both for gaps included and excluded)
        matches, mismatches, gaps = self.get_alignment_info(alignment)
        identity = matches / (matches + mismatches)
        matches, mismatches, _ = self.get_alignment_info(alignment, False)
        total = (matches + mismatches)
        no_gap_identity = 0 if total == 0 else matches / total
        return identity, no_gap_identity, gaps

    @staticmethod
    def perform_global_alignment_b(data_type, reference_sequence, sequence):
        score = False
        alignment = False
        alignment_result = ''
        # EMBL-EBI EMBOSS Needle configurations are used
        if (data_type == 'PROTEIN'):
            # Use configuration for protein sequence alignments
            alignments = pairwise2.align.globalds(reference_sequence, sequence, matlist.blosum62, -10, -0.5,
                                                  penalize_end_gaps=False)
        else:
            alignments = pairwise2.align.globalms(reference_sequence, sequence, 5, -4, -10, -0.5,
                                                  penalize_end_gaps=False)
        if (len(alignments) > 0):
            alignment = alignments[0]
            alignment_result = format_alignment(*alignments[0]).split('\n')
            score = alignment[2]

        return score, alignment, alignment_result

    @staticmethod
    def perform_global_alignment_p(data_type, reference_sequence, sequence):
        # EMBL-EBI EMBOSS Needle configuration except for gap extension penalty
        # as parasail accepts only integers for scoring
        if (data_type == 'PROTEIN'):
            # Use configuration for protein sequence alignments
            scoring_matrix = parasail.blosum62
        else:
            scoring_matrix = matrix_create(''.join(set(''.join([reference_sequence, sequence]))), 5, -4)
        result = parasail.sg_qx_trace_striped_16(reference_sequence, sequence, 10, 1, scoring_matrix)
        alignment_traceback = result.get_traceback()
        score = False
        alignment = []
        alignment_result = []
        if len(alignment_traceback.comp) != 0:  # if there is an alignment
            cigar = result.get_cigar()
            alignment = [alignment_traceback.query, alignment_traceback.ref, result.score, cigar.beg_ref,
                         result.end_ref+1, cigar.beg_query, result.end_query+1]
            alignment_result = [alignment_traceback.ref, alignment_traceback.comp, alignment_traceback.query]
            score = result.score

        return score, alignment, alignment_result

    def align(self, reference_sequence, sequence, data_type, config={}):
        score = False
        alignment = False
        alignment_result = ''
        # Remove any gap characters present in the sequences (PDB data sometimes might be incomplete)
        if ('removeGaps' in config and config['removeGaps'] is True):
            reference_sequence = reference_sequence.replace('-', '')
            sequence = sequence.replace('-', '')
        # Convert contiguous gap regions in the sequences to single gaps
        if ('compressGaps' in config and config['compressGaps'] is True):
            if ('-' in reference_sequence):
                reference_sequence = re.sub('-{2,}', '-', reference_sequence)
            if ('-' in sequence):
                sequence = re.sub('-{2,}', '-', sequence)
        # Use a different character to represent gaps in the sequences
        # and avoid conflict with the '-' character that is used as an
        # alignment gap character in the alignment result
        if ('gapChars' in config and len(config['gapChars']) > 0):
            reference_sequence = reference_sequence.replace('-', config['gapChars'][0])
            sequence = sequence.replace('-', config['gapChars'][1])
        try:
            # Skip long sequences that lead to long running processes
            if (len(sequence) < self._max_sequence_size):
                score, alignment, alignment_result = self.perform_global_alignment(data_type,
                                                                                   reference_sequence, sequence)
            elif (self.verbose is True):
                raise Exception(''.join(['Sequence length: ', repr(len(sequence))]))
        except:
            if (self.verbose is True):
                print(traceback.format_exc())

        return score, alignment, alignment_result

    @staticmethod
    def get_reference_sequence_alignment_b(alignment_result):
        # Find the aligned range and the number of gaps of the reference sequence
        # in the alignment result
        reference_line = alignment_result[0].strip().split(' ')
        # If there are no position details in the alignment result, find the first aligned element
        if(len(reference_line) == 1):
            alignment_line = alignment_result[1]
            start = alignment_line.find('|')
            reference_line = reference_line[0]
        else:
            start = int(reference_line[0]) - 1
            reference_line = reference_line[1]
        end = start + len(reference_line.replace('-', ''))
        gaps = len(reference_line) - len(reference_line.replace('-', ''))
        aligned_range = [start, end]
        return aligned_range, gaps

    @staticmethod
    def get_reference_sequence_alignment_p(alignment_result):
        gaps = alignment_result[0].count('-')
        aligned_range = [alignment_result[3], alignment_result[4]]

        return aligned_range, gaps

    def local_align(self, reference_sequence, sequence, gap_chars=[], compress_gaps=False):
        # Convert contiguous gap regions in the sequences to single gaps
        if (compress_gaps is True):
            if ('-' in reference_sequence):
                reference_sequence = re.sub('-{2,}', '-', reference_sequence)
            if ('-' in sequence):
                sequence = re.sub('-{2,}', '-', sequence)
        # Use a different character to represent gaps in the sequences
        # and avoid conflict with the '-' character that is used as an
        # alignment gap character in the alignment result
        if (len(gap_chars) > 0):
            reference_sequence = reference_sequence.replace('-', gap_chars[0])
            sequence = sequence.replace('-', gap_chars[1])
        return self.perform_local_alignment(reference_sequence, sequence)

    @staticmethod
    def perform_local_alignment_b(reference_sequence, sequence):
        score = False
        alignment = False
        alignment_result = False
        # EMBL-EBI EMBOSS Water configuration
        alignments = pairwise2.align.localms(reference_sequence, sequence, 5, -4, -10, -0.5)
        if (len(alignments) > 0):
            alignment = alignments[0]
            alignment_result = format_alignment(*alignments[0]).split('\n')
            score = alignment[2]
        return score, alignment, alignment_result

    @staticmethod
    def perform_local_alignment_bn(reference_sequence, sequence):
        score = False
        alignment = False
        alignment_result = False
        # EMBL-EBI EMBOSS Water configuration
        aligner = Align.PairwiseAligner()
        aligner.mode = 'local'
        aligner.match_score = 5.0
        aligner.mismatch_score = -4.0
        aligner.open_gap_score = -10
        aligner.extend_gap_score = -0.5
        aligner.target_end_gap_score = 0.0
        aligner.query_end_gap_score = 0.0
        alignments = aligner.align(reference_sequence, sequence)
        if (len(alignments) > 0):
            alignment_path = alignments[0].path
            alignment = [reference_sequence, sequence, alignments[0].score,
                         alignment_path[0][0], alignment_path[0][1],
                         alignment_path[1][0], alignment_path[1][1]]
            alignment_result = alignments[0].format().split('\n')
            score = alignments[0].score
        return score, alignment, alignment_result

    @staticmethod
    def perform_local_alignment_p(reference_sequence, sequence):
        # EMBL-EBI EMBOSS Water configuration except for gap extension penalty
        # as parasail accepts only integers for scoring
        scoring_matrix = matrix_create(''.join(set(''.join([reference_sequence, sequence]))), 10, -4)
        result = parasail.sw_trace_striped_16(reference_sequence, sequence, 12, 2, scoring_matrix)
        alignment_traceback = result.get_traceback()
        score = False
        alignment = []
        alignment_result = []
        if len(alignment_traceback.comp) != 0:  # if there is an alignment
            end_ref = result.end_ref + 1
            end_query = result.end_query + 1
            # Parasail returns bounds of traversed parts of sequences
            beg_ref = result.end_ref + 1 - len(alignment_traceback.ref.replace('-', ''))
            beg_query = result.end_query + 1 - len(alignment_traceback.query.replace('-', ''))
            alignment = [reference_sequence, sequence, result.score, beg_ref, end_ref, beg_query, end_query]
            alignment_result = [alignment_traceback.query, alignment_traceback.comp,
                                alignment_traceback.ref, beg_query, end_query]
            score = result.score

        return score, alignment, alignment_result

    def call_dssp(self, pdb_id):
        # Determine protein secondary structure using DSSP
        result = False
        try:
            output = subprocess.run(
                ['timeout', '120s', self.dssp_executable_path, '--output-format=dssp',
                 os.path.sep.join([self.pdb_dataset_path, ''.join([pdb_id, '.pdb'])])], capture_output=True,
                cwd=os.getcwd())
            output = output.stdout.decode("utf-8")
            output = output.split('\n')
            if (len(output) > 2):
                with open(''.join([self.dssp_cache_dir, os.path.sep, pdb_id, '.dssp']), 'w',
                          encoding='utf-8') as dssp_file:
                    dssp_file.write('\n'.join(output))
                result = True
        except:
            print(traceback.format_exc())

        return result

    def call_stride(self, pdb_id, chain_id, output_dir):
        # Determine protein secondary structure using STRIDE
        result = False
        try:
            output = subprocess.run(
                ['timeout', '120s', self.stride_executable_path, ''.join(['-r', chain_id]), '-o',
                 os.path.sep.join([self.pdb_dataset_path, ''.join([pdb_id, '.pdb'])])], capture_output=True,
                cwd=os.getcwd())
            output = output.stdout.decode("utf-8")
            output = output.split('\n')
            if (len(output) > 2):
                with open(''.join([output_dir, os.path.sep, pdb_id, '_', chain_id, '.stride']), 'w',
                          encoding='utf-8') as stride_file:
                    stride_file.write('\n'.join(output))
                result = True
        except:
            print(traceback.format_exc())

        return result

    def get_structures_by_stride(self, pdb_id, chain_id, residue_range, output_dir=''):
        available_residues = residue_range['fullRange']
        if(output_dir == ''):
            output_dir = self.stride_cache_dir
        secondary_structure = []
        primary_structure = []
        sec_structure = []
        prim_structure = []
        output_filepath = ''.join([output_dir, os.path.sep, pdb_id, '_', chain_id, '.stride'])
        alternative_output_filepath = ''.join([output_dir, os.path.sep, pdb_id, '.stride'])
        result = True
        if(os.path.exists(output_filepath) is False):
            if (os.path.exists(alternative_output_filepath) is False):
                result = self.call_stride(pdb_id, chain_id, output_dir)
            else:
                output_filepath = alternative_output_filepath
        # If there was a previous run on this protein, read the cached result
        if(result is True):
            with open(output_filepath, 'r', encoding='utf-8') as stride_file:
                # Collect secondary structure output
                last_residue = -1
                for line in stride_file:
                    if(line[0:3] == 'SEQ'):
                        last_residue = int(line[60:65])
                        prim_structure.append(line[10:60])
                    if (line[0:3] == 'STR'):
                        sec_structure.append(line[10:60])
                sec_structure = ''.join(sec_structure)
                prim_structure = ''.join(prim_structure)
                # Add loops as dots in the sequence.
                sec_structure = sec_structure[residue_range['discardedStart']:last_residue].replace(' ', '.')
                prim_structure = prim_structure[residue_range['discardedStart']:last_residue]

                # Add gaps between missing residues
                for list_position, residue_number in enumerate(range(available_residues[0],
                                                                     available_residues[-1] + 1)):
                    if (residue_number in available_residues and list_position < len(sec_structure)):
                        secondary_structure.append(sec_structure[list_position])
                        primary_structure.append(prim_structure[list_position])
                    else:
                        secondary_structure.append('-')
                        primary_structure.append('-')

        return ''.join(primary_structure), ''.join(secondary_structure)

    @staticmethod
    def get_alignment_positions(alignment, sequence_b, total_residue_range, last_position=0):
        # alignment list variable containing:
        # [0] - sequence A
        # [1] - sequence B
        # [2] - score
        # [3] - sequence B alignment start
        # [4] - sequence B alignment end
        # [5] - sequence A alignment start (parasail only)
        # [6] - sequence A alignment end (parasail only)

        # Aligned range refers to the residues of the PDB sequence
        aligned_range = [len(alignment[1][:alignment[3]].replace('-', '')) + last_position,
                         len(alignment[1][:alignment[4]].replace('-', '')) + last_position]
        aligned_sequence = sequence_b[aligned_range[0]:aligned_range[1]]
        range_length = len(aligned_sequence)

        # Actual range refers to the range in the whole known sequence
        actual_range = total_residue_range[aligned_range[0]:aligned_range[1]]
        if(len(aligned_range) > 1):
            actual_range = [actual_range[0], actual_range[-1]]
        else:
            actual_range = []

        return aligned_range, actual_range, range_length

    def get_all_sequences(self, structure_info, available_residues):
        # Extract all the available sequences for a PDB file (1D, 2D, hydrophobicity, mixed)
        sequences = defaultdict()
        sequences['secondary'], sequences['primary'] = self.get_structure(structure_info, available_residues)
        if(len(sequences['secondary']) > 1):
            sequences['mixed'], sequences['hydrophobicity'] = self.create_mixed_sequence(sequences['primary'],
                                                                                         sequences['secondary'])
        else:
            if (len(sequences['primary']) > 1):
                sequences['hydrophobicity'] = self.create_hydrophobicity_sequence(sequences['primary'])
            else:
                sequences = defaultdict()
        return sequences
