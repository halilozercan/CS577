REFERENCE="../files/MG1655_REF.fasta"
KMER_SIZE="17"
READ_FILE="../files/MG1655_1.fastq"
MAX_DISTANCE=6
READ_COUNT=100000

Runtime() {
    res1=$(date +%s.%N)

    $1

    res2=$(date +%s.%N)
    dt=$(echo "$res2 - $res1" | bc)
    dd=$(echo "$dt/86400" | bc)
    dt2=$(echo "$dt-86400*$dd" | bc)
    dh=$(echo "$dt2/3600" | bc)
    dt3=$(echo "$dt2-3600*$dh" | bc)
    dm=$(echo "$dt3/60" | bc)
    ds=$(echo "$dt3-60*$dm" | bc)

    printf "%s: %02d:%02d:%02.4f\n" $2 $dh $dm $ds
}

Runtime "python index.py $REFERENCE $KMER_SIZE" Indexing
Runtime "python seed.py $READ_FILE $KMER_SIZE $READ_COUNT" Seed
Runtime "python client.py ${REFERENCE}.${KMER_SIZE}.uih ${READ_FILE}.${KMER_SIZE}.seeds ${READ_FILE}.${KMER_SIZE}.client" Client
Runtime "python extend.py $READ_FILE ${READ_FILE}.${KMER_SIZE}.client $REFERENCE $KMER_SIZE $READ_COUNT $MAX_DISTANCE" Extension
