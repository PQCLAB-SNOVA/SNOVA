'''
Find parameters for irreducible S matrices
'''

import sys

q = 23

GF_q = GF(q)
x = GF_q.gen()
fetch_int = lambda x: x

params = {}

lset = [2, 3, 4, 5]

for L in lset:
    try:
        print(L)
        for i1 in range(0, q):
            for delta in range(0, 15):
                for i2 in range(0, q):
                    S1 = matrix(GF_q, L, L)
                    for i in range(L):
                        for j in range(i, L):
                            S1[i, j] = fetch_int((i1 + i + j) & delta)
                            S1[j, i] = S1[i, j]
                    S1[L - 1, L - 1] = fetch_int(i2)

                    if S1.charpoly().is_irreducible():
                        p = i1, delta, i2
                        params[p] = params.get(p, 0) + 1
    except Exception as exc:
        print(exc)

selection = []
for item in params.items():
    if item[1] == len(lset):
        selection.append(item[0])

if len(selection) > 0:
    print('Found', sorted(selection)[0], sorted(selection)[1])
else:
    print('No parameters found')
# quit()

# Check
print()
print('q=', q)

if q == 3:
    a = 2
    b = 1
    c = 1
elif q == 5:
    a = 3
    b = 2
    c = 1
elif q == 7:
    a = 4
    b = 6
    c = 1
elif q == 11:
    a = 0
    b = 3
    c = 6
elif q == 13:
    a = 2
    b = 11
    c = 3
elif q == 17:
    a = 1
    b = 11
    c = 10
elif q == 19:
    a = 1
    b = 3
    c = 15
elif q == 23:
    a = 1
    b = 11
    c = 22
elif q == 29:
    a = 3
    b = 12
    c = 11
elif q == 31:
    a = 2
    b = 5
    c = 8

for L in range(2, 6):
    if q == 7:
        if L == 3:
            c = 5
        else:
            c = 1

    print()
    print('L=', L)
    S1 = matrix(GF_q, L, L)
    for i in range(L):
        for j in range(L):
            S1[i, j] = fetch_int((a + i + j) & b)
    S1[L - 1, L - 1] = fetch_int(c)

    if not S1.charpoly().is_irreducible():
        print('Not irreducible')
        # raise Exception('Not irreducible')
    for l in range(1, 2):
        print(L, l)
        print(S1 ^ l)
