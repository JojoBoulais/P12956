import math
import sys
with open('mutated_sequences.txt', 'r') as file:
    lines = file.readlines()

print(math.log(5))

epsilon = sys.float_info.epsilon

print(len(lines))

aa_freq = {}
shannon = {}
acides_amines = list("ARNDCEQGHILKMFPSTWYV")

for pos in range(0, len(lines[0])-1):

    aa_freq[str(pos+558)] = {}
    shannon[str(pos+558)] = {}

    for aa in acides_amines:
        aa_freq[str(pos + 558)][aa] = 0

    for alignment in lines:
        #print(list(alignment)[pos])
        aa_freq[str(pos + 558)][list(alignment)[pos]] += 1


    shannon[str(pos+558)] = -sum([ (int(aa_freq[str(pos + 558)][aa])/float(len(lines)))
                                   * ((math.log(int(aa_freq[str(pos + 558)][aa])+epsilon)/float(len(lines))))
                                    for aa in acides_amines])

print(aa_freq)
print(shannon)
    # for aa in acides_amines:
    #
    #     aa_freq[str(pos)][aa] = 0