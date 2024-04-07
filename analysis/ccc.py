a={"atom_2": {"C": 2, "N": 0, "O": 0, "F": 0, "NO": 0, "OF": 0, "NF": 0, "CN": 1, "CO": 2, "CF": 0, "NOF": 0, "CNO": 0, "CNF": 0, "COF": 0, "CNOF": 0}, "atom_3": {"C": 3, "N": 0, "O": 0, "F": 0, "NO": 0, "OF": 0, "NF": 0, "CN": 1, "CO": 4, "CF": 0, "NOF": 0, "CNO": 1, "CNF": 0, "COF": 0, "CNOF": 0}, "atom_4": {"C": 7, "N": 0, "O": 0, "F": 0, "NO": 0, "OF": 0, "NF": 0, "CN": 5, "CO": 14, "CF": 0, "NOF": 0, "CNO": 5, "CNF": 0, "COF": 0, "CNOF": 0}, "atom_5": {"C": 17, "N": 0, "O": 0, "F": 0, "NO": 1, "OF": 0, "NF": 0, "CN": 19, "CO": 53, "CF": 1, "NOF": 0, "CNO": 37, "CNF": 0, "COF": 0, "CNOF": 0}, "atom_6": {"C": 56, "N": 0, "O": 0, "F": 0, "NO": 0, "OF": 0, "NF": 0, "CN": 95, "CO": 238, "CF": 2, "NOF": 0, "CNO": 207, "CNF": 2, "COF": 2, "CNOF": 0}, "atom_7": {"C": 196, "N": 0, "O": 0, "F": 0, "NO": 0, "OF": 0, "NF": 0, "CN": 407, "CO": 1183, "CF": 6, "NOF": 0, "CNO": 1283, "CNF": 18, "COF": 6, "CNOF": 0}, "atom_8": {"C": 832, "N": 0, "O": 0, "F": 0, "NO": 0, "OF": 0, "NF": 0, "CN": 2030, "CO": 6435, "CF": 20, "NOF": 0, "CNO": 8094, "CNF": 116, "COF": 38, "CNOF": 114}, "atom_9": {"C": 3781, "N": 0, "O": 0, "F": 0, "NO": 0, "OF": 0, "NF": 0, "CN": 10795, "CO": 37804, "CF": 0, "NOF": 0, "CNO": 53701, "CNF": 33, "COF": 13, "CNOF": 136}, "atom_10": {"C": 1282, "N": 0, "O": 0, "F": 0, "NO": 0, "OF": 0, "NF": 0, "CN": 20083, "CO": 14422, "CF": 3275, "NOF": 0, "CNO": 54855, "CNF": 11839, "COF": 11052, "CNOF": 15312}}
for i in a:
    tol=0
    aa=f'{i}'
    no=['O', 'F','OF','NF','NOF' ]
    for j in a[i]:
        #if j not in no:
        aa+=f'&{a[i][j]}'
        tol+=int(a[i][j])
    aa+=f'&\\textbf{{{tol}}}\\\\'
    print(aa)


ccc=''
for i in a['atom_2']:
    tol=0
    print(i)
    for j in a:
         tol+=int(a[j][i])
    ccc+=f'&\\textbf{{{tol}}}'
print(ccc)
