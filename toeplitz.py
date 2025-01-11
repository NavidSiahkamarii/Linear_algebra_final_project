import numpy as np
from scipy.sparse import diags



def create_cyclic_matrix(n):
    J = np.zeros((n, n), dtype=int)

    for i in range(n):
        J[i, (i + 2) % n] = 1

    return J


def create_toeplitz(n, main_diag, sub_diag_1, sub_diag_2, super_diag_1, super_diag_2):
    main = np.full(n, main_diag)  # قطر اصلی
    sub1 = np.full(n - 1, sub_diag_1)  # زیرقطر اول
    sub2 = np.full(n - 2, sub_diag_2)  # زیرقطر دوم
    super1 = np.full(n - 1, super_diag_1)  # بالاقطر اول
    super2 = np.full(n - 2, super_diag_2)  # بالاقطر دوم

    diagonals = [main, sub1, super1, sub2, super2]
    offsets = [0, -1, 1, -2, 2]  # موقعیت قطرها
    toeplitz_matrix = diags(diagonals, offsets, shape=(n, n)).toarray()

    return toeplitz_matrix


def backward_substitution(A, b):

    n = A.shape[0]
    x = np.zeros(n)

    for i in range(n - 1, -1, -1):
        x[i] = (b[i] - np.dot(A[i, i + 1:], x[i + 1:])) / A[i, i]

    return x

def toeplitz_algorithm(n,main_diag, sub_diag_1, sub_diag_2, super_diag_1, super_diag_2, b):
    A_tilda = create_cyclic_matrix(n) @ create_toeplitz(n, main_diag, sub_diag_1, sub_diag_2, super_diag_1, super_diag_2)
    A11 = A_tilda[:-2, :-2]  # بلوک اصلی
    w = A_tilda[-2, :-2]  # بردار w^T
    s = A_tilda[-1, :-2]  # بردار s^T
    p = A_tilda[:-2, -2]  # ستون p
    r = A_tilda[:-2, -1]  # ستون r
    x = np.zeros(n)

    u = backward_substitution(A11, b[2:])
    v = backward_substitution(A11, p)
    z = backward_substitution(A11, r)

    B = np.array([
            [np.dot(w, v), np.dot(w, z)],
            [np.dot(s, v), np.dot(s, z)]
        ])

    c = np.array([
        [np.dot(w, u) - b[0]],
        [np.dot(s, u) - b[1]]])

    solution = np.linalg.solve(B, c)
    x[n-2], x[n-1] = solution.flatten()
    x[:-2] = u - x[n-2] * v - x[n-1] * z

    return x

# print(toeplitz_algorithm(6, 4, 1, 2, 1, 2,np.full(6, 5) ))