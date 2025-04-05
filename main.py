import math
import random
from copy import deepcopy
import pandas as pd

df = pd.read_csv("msu081_Supplementary_Data/tableau_suplementaire_2.csv", sep=";")
PROTEIN = "EYSEEELKTHISKGTLGKFTVPMLKEACRAYGLKSGLKKQELLEALTKHFQD"
SAP_region = [558, 609]
N = 10**5

def calc_identity(mutated_sequence):

    sim = 0

    for i in range(len(mutated_sequence)):
        if mutated_sequence[i] == PROTEIN[i]:
            sim += 1

    return sim/len(mutated_sequence)

# Algorithme évolutif

mutated_sequences = []
random_aa = None

for _ in range(500):

    mutated_sequence = list(deepcopy(PROTEIN))
    while(calc_identity(mutated_sequence) > 0.50):

        s = "NA"
        k = 0

        # 1. Generate a random mutation
        random_position = random.randint(SAP_region[0], SAP_region[1])
        while (s == "NA" or float(s) < 0 or str(random_aa["mut_aa"]) == "*"):
            position = df[df["position"] == random_position]
            random_aa = position.iloc[random.randint(0, len(position)-1)]
            s = random_aa["fitness"]

        # 2. Calculate the probability of fixation
        s = float(s)
        pfix = (1 - math.e**(-2*s) ) / (1 - math.e**(-4*N*s))

        # 3. Determine if mutation is successful
        r = random.random()
        if r < pfix:

            mutated_sequence[random_position-SAP_region[0]] = str(random_aa["mut_aa"])

    mutated_sequences.append("".join(mutated_sequence))


# output sequences
outstream = open("mutated_sequences.txt", "w")
outstream.write("\n".join(mutated_sequences))
outstream.close()