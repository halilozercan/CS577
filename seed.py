import hashlib
import sys
from random import shuffle

from Bio import SeqIO

from cipher import AesCtrCipher

def seed_read(read, mer_size, key):
    hashes = []
    seed_end_positions = range(mer_size, len(read), mer_size)
    if seed_end_positions[-1] != len(read) - 1:
        seed_end_positions.append(len(read) - 1)
    for i in seed_end_positions:
        kmer = read[i - mer_size:i]
        hash = bytearray(hashlib.sha1(kmer).digest())
        keyed_hash = bytearray([a ^ b for (a, b) in zip(hash, key)])
        hashes.append(keyed_hash[:16])

    return hashes, seed_end_positions


def seed(filename, mer_size, hash_key='DEMO_KEY', read_key='DEMO_KEY'):
    key = bytearray(hashlib.sha1(hash_key).digest())
    cipher = AesCtrCipher(read_key)
    if '.fq' in  filename or '.fastq' in filename:
        file_sequences = SeqIO.parse(open(filename), 'fastq')
    elif '.fa' in  filename or '.fasta' in filename:
        file_sequences = SeqIO.parse(open(filename), 'fasta')
    read_hashes = []
    for i, sequence in enumerate(file_sequences):
        if i > 100000: break
        if i % 10000 == 0:
            print i, "reads indexed as seeds"
        seeds, seed_end_positions = seed_read(str(sequence.seq), mer_size, key)
        for seed, seed_end_position in zip(seeds, seed_end_positions):
            try:
                seed.extend(('%0.6X%0.2X%0.2X' % (i, seed_end_position-mer_size, 0)).decode('hex'))
            except:
                print 'Problems', i, seed_end_position-mer_size
            part_to_encrypt = str(seed)[10:]
            encrypted = bytearray(cipher.encrypt(part_to_encrypt))
            seed[10:] = encrypted
            read_hashes.append(seed)

        seeds, seed_end_positions = seed_read(str(sequence.reverse_complement().seq), mer_size, key)
        for seed, seed_end_position in zip(seeds, seed_end_positions):
            try:
                seed.extend(('%0.6X%0.2X%0.2X' % (i, seed_end_position - mer_size, 1)).decode('hex'))
            except:
                print 'Problems', i, seed_end_position - mer_size
            part_to_encrypt = str(seed)[10:]
            encrypted = bytearray(cipher.encrypt(part_to_encrypt))
            seed[10:] = encrypted
            read_hashes.append(seed)

        # print cipher.decrypt(array.array('B', read_hashes[j]).tostring()[10:])

    shuffle(read_hashes)

    with open(filename + '.' + str(mer_size) + '.seeds', 'wb') as f:
        for hash in read_hashes:
            f.write(hash)


if __name__ == "__main__":
    seed(sys.argv[1], int(sys.argv[2]))
