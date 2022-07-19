import time
import math
import random
from sympy import GF, invert
import logging
import numpy as np
from sympy.abc import x
from sympy import ZZ, Poly, degree
from sympy.polys.polyerrors import NotInvertible

def is_prime(n):
    for i in range(2, int(n ** 0.5) + 1):
        if n % i == 0:
            return False
    return True

def is_2_power(n):
    return n != 0 and (n & (n - 1) == 0)

def random_poly(length, d, neg_ones_diff=0):
    return Poly(np.random.permutation(
        np.concatenate((np.zeros(length - 2 * d - neg_ones_diff),
                        np.ones(d), -np.ones(d + neg_ones_diff)))),
                x).set_domain(ZZ)

def invert_poly(f_poly, R_poly, p):
    log = logging.getLogger("mathutils")
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

def shift_matriks(arr):
    F = np.zeros((len(arr), len(arr)))
    for i in range(len(arr)):
        for j in range(len(arr)):
            if i == 0:
                F[i,j] = arr[j]
            else:
                F[i,j] = F[(i-1)%len(arr)][(j-1)%len(arr)]
    return np.matrix(F)

def matrikskalimod(f,g,b): 
  result = g @ f 
  n = len(result)
  for i in range(n): 
    for j in range(n): 
      result[i,j] = np.remainder(result[i,j], b)
      if result[i,j] > b/2: 
        result[i,j] = result[i,j] - b  
  return result

def vektorkalimod(f,x,b): 
  result = x @ f
  n = result.size
  for i in range(0,n): 
    result[0,i] = np.remainder(result[0,i],b)
    if result[0,i] > b/2: 
      result[0,i] = result[0,i] - b  
  return result

def vektortambahmod(x,y,b): 
  result = np.add(x,y)
  n = result.size
  for i in range(n): 
    result[0,i] = np.remainder(result[0,i],b)
    if result[0,i] > b/2: 
      result[0,i] = result[0,i] - b  
  return result

def gen_keypriv(n,p,q):
    log = logging.getLogger("mathutils")
    R_poly = Poly(x**n - 1, x,domain='ZZ')
    tester = None
    g_poly = random_poly(n, int(math.sqrt(q)))
    deg_g = degree(g_poly,gen=x)
    g_arr = g_poly.all_coeffs()[::-1]
    if deg_g != n-1:
        for _ in range(n-1 - deg_g ):
            g_arr.append(int(0))
    g = shift_matriks(g_arr) 
    while tester is None:
        f_poly = random_poly(n, n // 3, neg_ones_diff=-1)
        try:
            f_p_poly = invert_poly(f_poly, R_poly, p)
            f_q_poly = invert_poly(f_poly, R_poly, q)
            tester =  f_p_poly * f_q_poly
        except NotInvertible as ex:
            log.debug(ex)
    if tester is None:
        raise Exception("Tidak bisa mengenerasi polinom f")
    deg_f = degree(f_poly,gen=x)
    f_arr = f_poly.all_coeffs()[::-1]
    if deg_f != n-1:
        for _ in range(n-1 - deg_f):
            f_arr.append(int(0))
    f = shift_matriks(f_arr)
    deg_f_p = degree(f_p_poly,gen=x)
    f_p_arr = f_p_poly.all_coeffs()[::-1]
    if deg_f_p != n-1:
        for _ in range(n-1 - deg_f_p):
            f_p_arr.append(int(0))
    f_p = shift_matriks(f_p_arr)
    deg_f_q = degree(f_q_poly,gen=x)
    f_q_arr = f_q_poly.all_coeffs()[::-1]
    if deg_f_q != n-1:
        for _ in range(n-1 - deg_f_q):
            f_q_arr.append(int(0))
    f_q = shift_matriks(f_q_arr)
    return f,g,f_p,f_q 

def gen_keypub(priv,q):
    f_q = priv[3]
    g = priv[1]
    h = matrikskalimod(f_q,g,q)
    return h 

def randomrgenerator(n, p):
  arr_r = []
  l = [-p,0,p]
  for i in range(n): 
    x = random.choice(l)
    arr_r.append(x)
  return np.array(arr_r)

def poly_r(n,p):
    r = randomrgenerator(n, p)
    return r 

def enkripsi(pub, m, r, q):
    size = m.shape[0]
    c = np.copy(m)
    for i in range(size): 
        c_temp = vektorkalimod(pub, r, q)
        c_cipherteks = vektortambahmod(m[i], c_temp, q)
        c[i] = c_cipherteks
    c = np.ravel(c)
    return c

def dekripsi(c,  priv, size, n, p, q):
    c_row_size = len(c) // n
    c = c.reshape(c_row_size, n)
    d = np.copy(c)
    for i in range(c_row_size):  
        m_temp1 = vektorkalimod(priv[0], c[i], q)
        m_plainteks = vektorkalimod(priv[2],m_temp1,p)
        d[i] = m_plainteks
    d = np.ravel(d)
    if size !=0:
        d = d[slice(-size)]
    d = bintostr(d)
    return d

def strToBinary(s):
    bin_conv = []
    arr_m = []
    byte_array = s.encode()
    binary_int = int.from_bytes(byte_array, "little")
    binary_string = bin(binary_int)
    for c in binary_string[2:]:
        bin_conv.append(c)
    for i in bin_conv: 
      for j in i:
        j = int(j)
        arr_m.append(j)
    return arr_m

def bintostr(result):
    arr = ''
    size = len(result)
    for i in range(size):
        int_ = int(result[i])
        str_ = str(int_)
        arr += str_
    binary_int = int(arr, 2)
    byte_number = binary_int.bit_length() + 7 // 8
    binary_array = binary_int.to_bytes(byte_number, byteorder="little")
    ascii_text = binary_array.decode()
    return ascii_text.replace('\x00','')

def custom_message_size(m, n):
    m_size_temp = len(m) % n
    if m_size_temp == 0:
        m_add_size = 0
    else:
        m_add_size = n - m_size_temp
    for _ in range(m_add_size):
        m.append(0)
    m_row_size = len(m) // n
    m_arr = np.array(m).reshape((m_row_size, n))
    return m_arr, m_add_size, m_row_size 

def readmessage(text_file):
    message='''''' 
    with open(text_file) as file:
        line = file.readlines()
    for text in line:
        message += text
    return message

def main(params, message):
    n,p,q = map(int, params)
    print(f"The public parameters are: n={n}, p={p}, q={q}")
    m, m_add_size = custom_message_size(strToBinary(message),n)[:2]
    start_time = time.time()
    keypriv = gen_keypriv(n,p,q) 
    keypub = gen_keypub(keypriv,q)
    r = poly_r(n, p)
    start_time_enkripsi = time.time()
    e = enkripsi(keypub,m,r,q)
    print(f"encryption result: \n {e}")
    print("--- encryption time : %s miliseconds ---\n" % 
        ((time.time() - start_time_enkripsi)*1000))
    start_time_dekripsi = time.time()
    d = dekripsi(e,  keypriv, m_add_size, n, p, q)
    print(f"decryption result:\n{d}")
    print("--- decryption time : %s miliseconds ---" % 
        ((time.time() - start_time_dekripsi)*1000))
    print("\n--- performance time (key generation + encryption + decryption time) : %s seconds---" % 
        (time.time() - start_time))

if __name__ == '__main__':
    params = list((251,3,2521)) 
    message = readmessage('./test_message.txt') 
    main(params, message)