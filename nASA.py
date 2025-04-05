import freesasa

structure = freesasa.Structure("1jeq.pdb")
result = freesasa.calc(structure)

A = result.residueAreas()["A"]
B = result.residueAreas()["B"]
SAP_region = [558, 609]

print(result.residueAreas())

for area in A.values():
    if int(area.residueNumber) > SAP_region[0] and int(area.residueNumber) < SAP_region[1]:
        print("-----------------------")
        print(area.residueNumber)
        print(area.residueType)
        print(area.relativeTotal) # This is the nASA (normalized ASA)


for area in B.values():
    if int(area.residueNumber) > SAP_region[0] and int(area.residueNumber) < SAP_region[1]:
        print("-----------------------")
        print(area.residueNumber)
        print(area.residueType)
        print(area.relativeTotal) # This is the nASA (normalized ASA)

print("-----------------------")