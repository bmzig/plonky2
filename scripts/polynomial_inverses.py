from sage.all import *

F = GF(2)
R = F['x']; (x,) = R._first_ngens(1)
K = F.extension(x**4 + x + 1, 'a')

print(K)

p = (x**2 + 1)
print(p)
q = p.inverse_mod(x**4 + x + 1)
print(q)

