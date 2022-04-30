import math
from sympy import GF, invert
import logging
import numpy as np
from sympy.abc import x
from sympy import ZZ, Poly, degree
from collections import Counter
import logging

log = logging.getLogger("mathutils")


def is_prime(n):
    for i in range(2, int(n ** 0.5) + 1):
        if n % i == 0:
            return False
    return True


def is_2_power(n):
    return n != 0 and (n & (n - 1) == 0)


def random_poly(length, d, neg_ones_diff=0):
    return Poly(np.random.permutation(
        np.concatenate((np.zeros(length - 2 * d - neg_ones_diff), np.ones(d), -np.ones(d + neg_ones_diff)))),
        x).set_domain(ZZ)


def invert_poly(f_poly, R_poly, p):
    inv_poly = None
    if is_prime(p):
        log.debug("Inverting as p={} is prime".format(p))
        inv_poly = invert(f_poly, R_poly, domain=GF(p))
    elif is_2_power(p):
        log.debug("Inverting as p={} is 2 power".format(p))
        inv_poly = invert(f_poly, R_poly, domain=GF(2))
        e = int(math.log(p, 2))
        for i in range(1, e):
            log.debug("Inversion({}): {}".format(i, inv_poly))
            inv_poly = ((2 * inv_poly - f_poly * inv_poly ** 2) % R_poly).trunc(p)
    else:
        raise Exception("Cannot invert polynomial in Z_{}".format(p))
    log.debug("Inversion: {}".format(inv_poly))
    return inv_poly


#test variable
n = 4
q = 41

f = Poly(x**3 -x + 1, x ,domain='ZZ')
#print(f.all_coeffs()[::-1])
#print(degree(f, gen=x))
#f = random_poly(n, n//3, neg_ones_diff=1)
R_poly = Poly(x**n - 1, x,domain='ZZ')
f_inv = invert_poly(f, R_poly, q)
#list_ = np.array(f_inv.all_coeffs()[::-1])

#f_poly = invert_poly(f, Poly(x**251+1, x, domain='ZZ'), 2521)
#log.info("f coeffs: {}".format(Counter(f_poly.coeffs())))
#f = Poly(-1 + x  + x**2,x,domain='ZZ')

#print(f"f : {f}")
print(f"f_inv di {q} : {f_inv}")
#print(f"coeffs of f_inv sebagai array : {list_}")
#print(f" hasil kali f dan f_inv : {((f * f_inv) % R_poly).trunc(q)}")
#print(type(list_))
