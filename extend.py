import array
import pickle
import sys

import editdistance
from Bio import SeqIO

from cipher import AesCtrCipher


def sxor(s1, s2):
    sb1 = array.array('B', s1)
    sb2 = array.array('B', s2)
    result = array.array('B', [])
    for i in range(len(sb1)):
        result.append(sb1[i] ^ sb2[i])
    return result


reads_file = sys.argv[1]
client_results = sys.argv[2]
reference_file = sys.argv[3]
kmer_size = int(sys.argv[4])
max_edit_threshold = int(sys.argv[5])

kmerhash_position_map = pickle.load(open(reference_file + '.' + str(kmer_size) + '.rih', 'rb'))
seeds_file = reads_file + '.' + str(kmer_size) + '.seeds'
output = reads_file + '.' + str(kmer_size) + '.align'

key = 'DEMO_KEY'

read_cipher = AesCtrCipher(key)

if '.fq' in reads_file or '.fastq' in reads_file:
    read_sequences = SeqIO.parse(open(reads_file), 'fastq')
elif '.fa' in reads_file or '.fasta' in reads_file:
    read_sequences = SeqIO.parse(open(reads_file), 'fasta')
reads = []
read_headers = []
for i, sequence in enumerate(read_sequences):
    if i > 100000: break
    reads.append(sequence)
    read_headers.append(sequence.id)

reference_chromosomes = SeqIO.parse(open(reference_file), 'fasta')
reference = []
for sequence in reference_chromosomes:
    reference.append(str(sequence.seq))

seeds = []
with open(seeds_file, 'rb') as f:
    while True:
        i = f.read(37)
        if i:
            seeds.append(i)
        else:
            break

alignment_results = [{'score': 1000, 'location': -1, 'chr': 'chr0'}] * len(reads)

nunique = 0
extend_success = 0
total = 0

with open(client_results, 'rb') as f:
    while True:
        seed = f.read(37)
        if seed is None:
            break
        ref_kmer = f.read(16)
        if len(ref_kmer) != 16:
            break
        # Check if this is a match
        if ref_kmer != ('0' * 16) and len(ref_kmer) == 16:
            total += 1
            seed_second_part = read_cipher.decrypt(seed[10:])
            reference_position = sxor(seed_second_part[:6], ref_kmer[10:]).tostring().encode('hex')
            unique = int(reference_position[0:2], 16)
            chromosome_index = int(reference_position[2:4], 16)
            in_chromosome_position = int(reference_position[4:], 16)

            read_number = int(seed_second_part[6:9].encode('hex'), 16)
            in_read_position = int(seed_second_part[9].encode('hex'), 16)
            is_reverse = int(seed_second_part[10].encode('hex'), 16)

            # print seed_second_part, unique, chromosome_index, in_chromosome_position, read_number, in_read_position

            if is_reverse:
                corresponding_read = str(reads[read_number].reverse_complement().seq)
            else:
                corresponding_read = str(reads[read_number].seq)

            start = in_chromosome_position - in_read_position
            end = start + len(corresponding_read)

            aligned_str = reference[chromosome_index][start:end]

            distance = editdistance.eval(corresponding_read, aligned_str)
            if distance <= max_edit_threshold and alignment_results[read_number]['score'] > distance:
                alignment_results[read_number] = {'score': distance,
                                                  'location': start,
                                                  'chr': chromosome_index}
            if not unique:
                nunique += 1
                candidates = kmerhash_position_map[str(seed[:10])]
                for (chromosome_index, in_chromosome_position) in candidates:

                    start = in_chromosome_position - in_read_position
                    end = start + len(corresponding_read)

                    aligned_str = reference[chromosome_index][start:end]

                    distance = editdistance.eval(corresponding_read, aligned_str)
                    if distance <= max_edit_threshold and alignment_results[read_number]['score'] > distance:
                        alignment_results[read_number] = {'score': distance,
                                                          'location': start,
                                                          'chr': chromosome_index}
            if distance <= max_edit_threshold:
                extend_success += 1

with open(output, 'w') as g:
    aligned = 0
    naligned = 0
    for i, result in enumerate(alignment_results):
        g.write(read_headers[i] + '\n')
        g.write(str(reads[i].seq) + '\n')
        if result['score'] != 1000:
            aligned += 1
            g.write(reference[result['chr']][result['location']:result['location'] + len(reads[i])] + '\n')
            g.write(str(result['score']) + '-' + str(result['location']) + '-' + str(result['chr']) + '\n')
        else:
            naligned += 1
            g.write('-\n')
            g.write('-\n')
        g.write('+\n')

print aligned, naligned, nunique
