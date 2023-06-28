
import cProfile   #calcul complexite
from Bio import SeqIO #biopython
import os


##########################################################################

#############################
######  DC3 Functions  ######
#############################

def DC3(S, verbose = False):
    '''
    Difference Cover size 3 is an algorithm that enables us
    to compute the suffix array of a string in a very efficient
    way.

    Input: String S (Genomic sequence : A,C,G and T)
           Verbose (True if you wish to visualize suffix array)

    Output: The suffix array I corresponding to the ordered list
            of the indexes of the first character of each
            suffixes.
    '''

    T = ACGTtoNum(S) # Translates String in a numeric sequence

    I = RecSort(T) # Recursive sort until we get final suffix array

    if(verbose):
        for i in I:
            print(S[i:])

    return I

def RecSort(T):

    #Appending three sentinel numbers (0 here)
    T.append(0)
    T.append(0)
    T.append(0)

    #Computing the P12 indexes
    size = len(T)
    P12 = P12construct(size)

    # Sorting triplets at index in P12 (returns sorted P12)
    index12 = radixTriplets(T, P12)

    # R12 is the list of ordered triplets
    R12 = R12construct(T, index12)

    # Computes the order of R12
    Order12 = Order(R12)

    # Test if we need to enter a new sorting loop
    if(len(R12) != Order12[-1]):

        # Computes new sequence to be sorted

        dicI12 = {index12[i]: i for i in range(len(index12))}

        newT = [0]*len(Order12)
        for i, num in enumerate(P12):
            newT[i] = Order12[dicI12.get(num)]

        # New sprting loop
        I = RecSort(newT)

        # Mapping indexes of I012 to higher order
        index12 = []
        for i in I:
            index12.append(P12[i])


    P0 = P0construct(size)

    # Constructs R0
    R0 = []
    dicI12 = {index12[i]: i for i in range(len(index12))}
    for num in P0:
        R0.append(T[num])
        R0.append(dicI12.get(num+1))

    # Computes index0 based on R0
    index0 = radix2(R0, [2*i for i in range(int(len(R0)/2))])
    index0 = [int(index0[i]*3/2) for i in range(len(index0))]

    return mergeIndexes(T, index12, index0)

def P12construct(size):
    '''Constructs P12 indexes based on the size of the array'''
    P1 = [i for i in range(size-2) if i%3 == 1]
    P2 = [i for i in range(size-2) if i%3 == 2]

    return P1 + P2

def P0construct(size):
    '''Constructs P12 indexes based on the size of the array'''
    P0 = [i for i in range(size-3) if i%3 == 0]

    return P0

def Order(T):
    n = 1
    orders = [1]
    for i in range(1,len(T)):
        if (T[i-1] == T[i]):
            orders.append(n)
        else:
            n += 1
            orders.append(n)
    return orders

def R12construct(T, idx):
    out = []
    for i in idx:
        out.append(triplet(T,i))
    return out


def triplet(T,i):
    return [T[i], T[i+1], T[i+2]]


def radix2(T, index):
    index = radix(T, index, 1)
    index = radix(T, index, 0)
    return index


def radixTriplets(T, index):
    index = radix(T, index, 2)
    index = radix(T, index, 1)
    index = radix(T, index, 0)
    return index


def radix(T, index, place):
    sortChar = [T[i+place] for i in index]
    counts = Counts(sortChar)
    sums = cumSum(counts)
    output = [0] * len(index)

    for i in index:
        num = T[i + place]
        output[sums[num]] = i
        sums[num] += 1

    return output

def Counts(T):
    '''Counts the number of occurences of each character'''

    c = [0] * (max(T)+1)
    for i in T:

        c[i] += 1
    return c


def compare(T, a, b, dic):

    if T[a] < T[b]:
        return True

    if T[a] > T[b]:
        return False

    if a % 3 != 0 and b % 3 != 0:
        return dic[a] < dic[b]

    return compare(T, a+1, b+1, dic)


def mergeIndexes(T, I12, I0):
    i = 0
    j = 0

    I012 = []
    I12dic = {I12[i] : i for i in range(len(I12))}

    while i < len(I12) and j < len(I0):
        if compare(T, I12[i], I0[j], I12dic):
            I012.append(I12[i])
            #print("added Val", I12[i])
            i += 1

        else:
            I012.append(I0[j])
            j += 1

    I012.extend(I12[i:])
    I012.extend(I0[j:])
    if(T[I012[0]] == 0):
        I012.pop(0)

    return I012


def cumSum(counts):
    '''Computes the cumulative sums based on counts'''
    res = [0] * len(counts)
    sums = 0
    for i, n in enumerate(counts):
        res[i] = sums
        sums+=n
    return res


def ACGTtoNum(S):
    dic = {"$" : 1,"A" : 2, "C" : 3, "G" : 4, "T" : 5}

    num_S = []
    for i in range(len(S)):
        num_S.append(dic[S[i]])

    return num_S

##########################################################################

###############################
######  Data management  ######
###############################

def getGenome(dataPath):
    genome=[]

    for seq_record in SeqIO.parse(dataPath, "fasta"):
        a={}
        a[seq_record.id]=str(seq_record.seq)[5:-2]
        genome.append(a)

    return genome

def storeChr(genome, number):

    file = "chr"+str(number)+".txt"
    f = open(file, "w")
    S=genome[number][list(genome[number].keys())[0]]
    S=S.upper()
    S=str(S)
    I = DC3(S)
    for i in I:
        f.write(str(i)+"\n")

    f.close()

def getI012(num):
    print(oui)



##########################################################################

############################
######  MAIN Program  ######
############################

genome=getGenome("GCF_000002765.5_GCA_000002765_genomic.fna")

for i in range(len(genome)):
    print(genome[i].keys())

S=genome[1]["NC_037280.1"]
S=S.upper()
s=str(S)

print("Number of chromosomes : ", len(genome))

#res1 = DC3(s)

for i in range(15):
    storeChr(genome, i)
    print("Done number ",i)
