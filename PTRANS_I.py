import numpy as np


def solve_pentadiagonal(n, d, a, b, c, e, f):

    # Initialize variables
    mu = np.zeros(n + 1)
    alpha = np.zeros(n + 1)
    beta = np.zeros(n + 1)
    gamma = np.zeros(n + 1)
    z = np.zeros(n + 1)

    # Step 3: Initialization
    mu[1] = d[1]
    alpha[1] = a[1] / mu[1]
    beta[1] = b[1] / mu[1]
    z[1] = f[1] / mu[1]

    # Step 4: Second row
    gamma[2] = c[2]
    mu[2] = d[2] - alpha[1] * gamma[2]
    alpha[2] = (a[2] - beta[1] * gamma[2]) / mu[2]
    beta[2] = b[2] / mu[2]
    z[2] = (f[2] - z[1] * gamma[2]) / mu[2]

    # Step 5: Forward elimination for rows 3 to n-2
    for i in range(3, n - 1):
        gamma[i] = c[i] - alpha[i-2] * e[i]
        mu[i] = d[i] - beta[i - 2] * e[i] - alpha[i - 1] * gamma[i]
        alpha[i] = (a[i] - beta[i - 1] * gamma[i]) / mu[i]
        beta[i] = b[i] / mu[i]
        z[i] = (f[i] - z[i - 2] * e[i] - z[i - 1] * gamma[i]) / mu[i]

    # Step 5: Handle the last two rows
    gamma[n - 1] = c[n - 1] - alpha[n - 3] * e[n - 1]
    mu[n - 1] = d[n - 1] - beta[n - 3] * e[n - 1] - alpha[n - 2] * gamma[n - 1]
    alpha[n - 1] = (a[n - 1] - beta[n - 2] * gamma[n - 1]) / mu[n - 1]
    z[n - 1] = (f[n - 1] - z[n - 3] * e[n - 1] - z[n - 2] * gamma[n - 1]) / mu[n - 1]

    gamma[n] = c[n] - alpha[n - 2] * e[n]
    mu[n] = d[n] - beta[n - 2] * e[n] - alpha[n - 1] * gamma[n]
    z[n] = (f[n] - z[n - 2] * e[n] - z[n - 1] * gamma[n]) / mu[n]

    # Step 6: Backward substitution
    x = np.zeros(n+1)

    x[n] = z[n]
    x[n - 1] = z[n - 1] - alpha[n - 1] * x[n]

    for i in range(n - 2, 0, -1):
        x[i] = z[i] - alpha[i] * x[i + 1] - beta[i] * x[i + 2]

    return x


# Example usage
# n = 6  # Size of the matrix
# d = np.array([0, 4, 4, 4, 4, 4, 4])  # Main diagonal
# a = np.array([0, 1, 1, 1, 1, 1, 0])  # Subdiagonal 1
# b = np.array([0, 2, 2, 2, 2, 0, 0])  # Subdiagonal 2
# c = np.array([0, 0, 1, 1, 1, 1, 1])  # Superdiagonal 1
# e = np.array([0, 0, 0, 2, 2, 2, 2])  # Superdiagonal 2
# f = np.array([0, 5, 5, 5, 5, 5, 5])  # Right-hand side
#
# solution = solve_pentadiagonal(n, d, a, b, c, e, f)
