import copy
import hashlib
import pickle
import sys
import heapq

import array
from collections import defaultdict

from Bio import SeqIO


def heap_sort(iterable):
    h = []
    for value in iterable:
        heapq.heappush(h, value)
    return [heapq.heappop(h) for i in range(len(h))]


def is_valid_base(c):
    return not (c != 'A' and c != 'C' and c != 'G' and c != 'T' and c != 'N')


def index(filename, mer_size, key="DEMO_KEY"):
    key = bytearray(hashlib.sha1(key).digest())

    # Key hashed kmer is key.
    # Value is a list of positions where this kmer appears
    kmerhash_position_map = defaultdict(list)

    fasta_sequences = SeqIO.parse(open(filename), 'fasta')
    for seq_index, sequence in enumerate(fasta_sequences):
        name, data = sequence.id, str(sequence.seq)
        print 'Reading', name
        for i in range(mer_size, len(data)):
            if i % 1000000 == 0 or i == (len(data) - 1):
                print "Hashed", i, "kmers"
            kmer = data[i - mer_size:i]
            hash = bytearray(hashlib.sha1(kmer).digest())
            keyed_hash = bytearray([a ^ b for (a, b) in zip(hash, key)])
            kmerhash_position_map[str(keyed_hash)].append((seq_index, 0, i - mer_size))  # 0 for not reverse

        data = str(sequence.reverse_complement().seq)
        print 'Reading', name, 'in reverse'
        for i in range(mer_size, len(data)):
            if i % 1000000 == 0 or i == (len(data) - 1):
                print "Hashed", i, "kmers"
            kmer = data[i - mer_size:i]
            hash = bytearray(hashlib.sha1(kmer).digest())
            keyed_hash = bytearray([a ^ b for (a, b) in zip(hash, key)])
            kmerhash_position_map[str(keyed_hash)].append((seq_index, 1, i - mer_size))  # 1 for reverse

    # hashes = heap_sort(hashes)
    unique_indexed_hashes = []
    extend_map = {}

    for keyed_hash, positions in kmerhash_position_map.iteritems():
        unique = len(positions) == 1
        info_to_hide = bytearray(('%0.2X%0.2X%0.2X%0.6X' %
                                  (unique, positions[0][0], positions[0][1], positions[0][2])).decode('hex'))

        new_hash = bytearray(keyed_hash)
        for j in range(len(info_to_hide)):
            new_hash[j + 10] ^= info_to_hide[j]

        unique_indexed_hashes.append(new_hash[:16])
        del new_hash

        if not unique:
            del positions[0]
            extend_map[keyed_hash[:10]] = positions

    print len(unique_indexed_hashes), 'unique'

    with open(filename + '.' + str(mer_size) + '.uih', 'wb') as f:
        for hash in unique_indexed_hashes:
            f.write(hash)

    with open(filename + '.' + str(mer_size) + '.rih', 'wb') as f:
        pickle.dump(extend_map, f)


if __name__ == "__main__":
    index(sys.argv[1], int(sys.argv[2]))
