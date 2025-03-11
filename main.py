import math
import random


# INIT
with open("P12956.fasta", "r") as f:
    lines = f.readlines()
    o_sequence = "".join([l.strip() for l in lines[1:]])

with open("deep_mutational_scan.tsv", "r") as f:
    lines = f.readlines()
    dms_labels = lines[0].strip().split("\t")
    dms = [l.strip().split("\t") for l in lines[1:]]

DEEPMUTDATA = dms
PROTEIN = o_sequence
N = 10**5
s = 5


# Algorithme évolutif


pfix = (1 - math.e**(-2*s) ) / (1 - math.e**(-4*N*s))
r = random.random()