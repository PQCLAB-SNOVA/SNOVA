GF16 = GF(16, 'x')
x = GF16.gen()

fetch_int = lambda x: GF16._cache.fetch_int(x)

for L in range(2, 7):
    if L > 4:
        S1 = matrix(
            GF16, L, L, lambda i, j: fetch_int(19 - 2 * L)
            if i == L - 1 and j == L - 1 else fetch_int(abs(8 - (i + j))))
    else:
        S1 = matrix(GF16, L, L, lambda i, j: fetch_int(abs(8 - (i + j))))

    # print(S1)
    print(L, S1.charpoly('z'))
    print(L, S1.charpoly().is_irreducible())
