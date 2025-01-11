import numpy as np

def solve_pentadiagonal_v2(n, d, a, b, c, e, y):
    # Initialize variables
    psi = np.zeros(n + 1)
    sigma = np.zeros(n + 1)
    phi = np.zeros(n + 1)
    w = np.zeros(n + 1)
    rho = np.zeros(n + 1)

    # Step 3: Initialization
    psi[n] = d[n]
    sigma[n] = c[n] / psi[n]
    phi[n] = e[n] / psi[n]
    w[n] = y[n] / psi[n]

    # Step 4: Handle n-1 row
    rho[n - 1] = a[n - 1]
    psi[n - 1] = d[n - 1] - sigma[n] * rho[n - 1]
    sigma[n - 1] = (c[n - 1] - phi[n] * rho[n - 1]) / psi[n - 1]
    phi[n - 1] = e[n - 1] / psi[n - 1]
    w[n - 1] = (y[n - 1] - w[n] * rho[n - 1]) / psi[n - 1]

    # Step 5: Backward elimination for rows n-2 to 3
    for i in range(n - 2, 2, -1):
        rho[i] = a[i] - sigma[i + 2] * b[i]
        psi[i] = d[i] - phi[i + 2] * b[i] - sigma[i + 1] * rho[i]
        sigma[i] = (c[i] - phi[i + 1] * rho[i]) / psi[i]
        phi[i] = e[i] / psi[i]
        w[i] = (y[i] - w[i + 2] * b[i] - w[i + 1] * rho[i]) / psi[i]

    # Handle row 2
    rho[2] = a[2] - sigma[4] * b[2]
    psi[2] = d[2] - phi[4] * b[2] - sigma[3] * rho[2]
    sigma[2] = (c[2] - phi[3] * rho[2]) / psi[2]
    phi[2] = e[2] / psi[2]
    w[2] = (y[2] - w[4] * b[2] - w[3] * rho[2]) / psi[2]

    # Handle row 1
    rho[1] = a[1] - sigma[3] * b[1]
    psi[1] = d[1] - phi[3] * b[1] - sigma[2] * rho[1]
    w[1] = (y[1] - w[3] * b[1] - w[2] * rho[1]) / psi[1]

    # Step 6: Backward substitution to compute solution vector X
    x = np.zeros(n + 1)
    x[1] = w[1]
    x[2] = w[2] - sigma[2] * x[1]

    for i in range(3, n + 1):
        x[i] = w[i] - sigma[i] * x[i - 1] - phi[i] * x[i - 2]

    return x[1:]  # Return solution vector excluding the dummy 0th element


# Example usage
# n = 6  # Size of the matrix
# d = np.array([0, 4, 4, 4, 4, 4, 4])  # Main diagonal
# a = np.array([0, 1, 1, 1, 1, 1, 0])  # Subdiagonal 1
# b = np.array([0, 2, 2, 2, 2, 0, 0])  # Subdiagonal 2
# c = np.array([0, 0, 1, 1, 1, 1, 1])  # Superdiagonal 1
# e = np.array([0, 0, 0, 2, 2, 2, 2])  # Superdiagonal 2
# y = np.array([0, 5, 5, 5, 5, 5, 5])  # Right-hand side
#
# solution = solve_pentadiagonal_v2(n, d, a, b, c, e, y)




# pentadiagonal_matrix = create_toeplitz(1000, 5, -1, -2, 1, 2)
# print("5x5 portion of the pentadiagonal matrix:")
# print(pentadiagonal_matrix[:5, :5])
