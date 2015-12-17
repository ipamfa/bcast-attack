from kannan import select_basis

# Test       
X = np.array([
    [3, -1, 5],
    [-5, 2, -1],
    [-3, 9, 2]
])
Y = np.array([
    [-11,  12,  -4],
    [3, -1, 5],
    [-5, 2, -1],
    [-3, 9, 2]
])
Z = np.array([
    [13, 15, -5,  8],
    [8, 6, 2, 6],
    [2, 4, 8, 4],
    [1, -1, -6, -8],
    [2, 6, -9, 6],
])

select_basis(5, Z)
