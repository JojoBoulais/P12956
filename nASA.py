import freesasa

structure = freesasa.Structure("1jeq.pdb")
result = freesasa.calc(structure)
print(result.residueAreas())


a = result.residueAreas()["A"]

for area in a.values():
    print(area.total)


