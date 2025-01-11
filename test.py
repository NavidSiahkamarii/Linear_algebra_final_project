import numpy as np
import time
from toeplitz import toeplitz_algorithm
from PTRANS_II import solve_pentadiagonal_v2
from PTRANS_I import solve_pentadiagonal

n = 500
main_diag = 0.5
sub_diag_1 = 7
sub_diag_2 = 3.5
super_diag_1 = 4
super_diag_2 = 1.8
b = np.full(n, 1)

d = np.full(n+1, main_diag)
a = np.full(n+1, sub_diag_1)
b = np.full(n+1, sub_diag_2)
c = np.full(n+1, super_diag_1)
e = np.full(n+1, super_diag_2)

temp_time = time.time()
solve_pentadiagonal(n, d, a, b, c, e, b)
print("time of PTRANS-I :",time.time()-temp_time)

temp_time = time.time()
solve_pentadiagonal_v2(n, d , a, b, c, e, b)
print("time of PTRANS-II :",time.time()-temp_time)

temp_time = time.time()
toeplitz_algorithm(n, main_diag, sub_diag_1, sub_diag_2, super_diag_1, super_diag_2, b)
print("time of toeplitz :",time.time()-temp_time)


