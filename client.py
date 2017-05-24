import sys

reference_kmers = []
with open(sys.argv[1], 'rb') as f:
    while True:
        i = f.read(16)
        if i:
            reference_kmers.append(i)
        else:
            break

reference_kmers.sort()

seeds = []
with open(sys.argv[2], 'rb') as f:
    while True:
        i = f.read(36)
        if i:
            seeds.append(i)
        else:
            break

seeds.sort()
print len(seeds)


def binary_search(iterable, item, low, high, same_len=False):
    if low > high or low >= len(iterable):
        return -1
    a = iterable[low]
    if same_len:
        a = a[:len(item)]
    if low == high and a == item:
        return low
    elif low == high and a != item:
        return -1
    else:
        middle = (low + high) / 2
        a = iterable[middle]
        if same_len:
            a = a[:len(item)]
        if item <= a:
            return binary_search(iterable, item, low, middle, same_len=same_len)
        elif item > a:
            return binary_search(iterable, item, middle + 1, high, same_len=same_len)


def sequential_search(iterable, query, same_len=False):
    for index, item in enumerate(iterable):
        if same_len:
            if query == item[:len(query)]:
                return index
        else:
            if query == item:
                return index
    return -1


with open(sys.argv[3], 'w') as f:
    found = 0
    for i, seed in enumerate(seeds):
        if i % 1000000 == 0:
            print i, 'matching', found
        index = binary_search(reference_kmers, seed[:10], 0, len(reference_kmers), same_len=True)
        #index = sequential_search(reference_kmers, seed[:10], same_len=True)
        if index >= 0:
            f.write(seed)
            f.write(reference_kmers[index])
            found += 1
        else:
            f.write(seed)
            f.write("0"*16)
