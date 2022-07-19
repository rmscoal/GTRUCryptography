import time
import math
import random
import psutil
from sympy import GF, invert
import logging
import numpy as np
from sympy.abc import x
from sympy import ZZ, Poly, degree
from sympy.polys.polyerrors import NotInvertible
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
        np.concatenate((np.zeros(length - 2 * d - neg_ones_diff),
                        np.ones(d), -np.ones(d + neg_ones_diff)))),
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

def shift_matriks(arr):
    F = np.zeros((len(arr), len(arr)))
    for i in range(len(arr)):
        for j in range(len(arr)):
            if i == 0:
                F[i,j] = arr[j]
            else:
                F[i,j] = F[(i-1)%len(arr)][(j-1)%len(arr)]
    return np.matrix(F)

def invMatriks(inv_minor, f, b):
    size_f = len(f)
    size_minor = len(inv_minor)
    g_col = [1/f[0,0]]
    for i in range(size_minor): 
        g = 0
        for j in range(size_minor):
            g += inv_minor[i,j] * f[j+1,0]
        if i >= size_minor - 2:
            delta = inv_minor[i, size_minor - 2] * (inv_minor[i, size_minor-2] - 1)/2 * f[size_minor-1, size_minor-1]*f[size_minor-1, size_minor]
            delta += inv_minor[i, size_minor -1]*(inv_minor[i,size_minor-1] -1)/2 * f[size_minor,size_minor-1]*f[size_minor,size_minor]
            delta += inv_minor[i,size_minor-2] * inv_minor[i,size_minor-1]*f[size_minor-1,size_minor]*f[size_minor,size_minor-1]
            g += delta
        g  = -g_col[0]* g
        g_col.append(g)
    row = np.zeros((1,size_minor))
    inv = np.concatenate((row, inv_minor), axis=0)
    inv = np.concatenate((np.matrix(g_col).T, inv),axis=1)
    for i in range(size_f):
        for j in range(size_f):
            inv[i,j] = np.remainder(inv[i,j], b)
            if inv[i,j] > b/2:
                inv[i,j] = inv[i,j] - b
    return inv

def matriks_mul(f,g,b):
  mul = g @ f
  size = len(f)
  for i in range(size-2, size):
    delta  = g[i,size-2] * (g[i,size-2] -1) * 1/2 * f[size-2,size-2] * f[size-2, size-1]
    delta += g[i,size-1] * (g[i,size-1] -1) * 1/2 * f[size-1,size-2] * f[size-1,size-1]
    delta += g[i,size-2] * g[i,size-1] * f[size-2, size-1] * f[size-1,size-2]
    mul[i,0] += delta
    for i in range(size): 
        for j in range(size): 
            mul[i,j] = np.remainder(mul[i,j], b)
            if mul[i,j] > b/2: 
                mul[i,j] = mul[i,j] - b
  return mul

def vektor_mul(f,x,b):
    x  = np.matrix(x)
    mul = x @ f
    size = mul.size
    delta = x[0,size-2] * (x[0,size-2] -1) * 1/2 * f[size-2,size-2] * f[size-2,size-1]
    delta += x[0,size-1] * (x[0,size - 1]-1) *1/2 * f[size-1,size-2] * f[size-1,size-1]
    delta += x[0,size-2]*x[0,size-1]*f[size-2,size-1]*f[size-1,size-2]
    mul[0,0] += delta
    for i in range(0,size): 
        mul[0,i] = np.remainder(mul[0,i],b)
        if mul[0,i] > b/2: 
            mul[0,i] = mul[0,i] - b
    return mul

def vektor_add(x,y,b):
    x,y = np.matrix(x),np.matrix(y)
    add = np.add(x,y)
    size = add.size
    add[0,0] += x[0,size-1] * y[0,size-2]
    for i in range(0,size): 
        add[0,i] = np.remainder(add[0,i],b)
        if add[0,i] > b/2: 
            add[0,i] = add[0,i] - b
    return add

def gen_keypriv(n,p,q):
    R_poly = Poly(x**(n-3) - 1, x,domain='ZZ')
    g_poly = random_poly(n-3, int(math.sqrt(math.sqrt(q))))
    deg_g = degree(g_poly,gen=x)
    g_arr = g_poly.all_coeffs()[::-1]
    if deg_g != n-4:
        for i in range(n-4 - deg_g ):
            g_arr.append(int(0))
    little_g = shift_matriks(g_arr)
    g_2x2 = np.zeros((2,2))
    for i in range(2):
        for j in range(2):
            if i == 1 and j == 0:
                continue
            else:
                g_2x2[i,j] = np.random.choice([-1,1], p=[0.5,0.5])
    right_g = np.zeros((n-3,2))
    right_g_and_g2x2 = np.concatenate((right_g, g_2x2), axis=0)
    lower_g = np.zeros((2,n-3))
    lower_g_and_little_g =  np.concatenate((little_g, lower_g), axis=0)
    g_minor = np.concatenate((lower_g_and_little_g, right_g_and_g2x2), axis=1)
    left_g = []
    for i in range(n-1):
        if i >= n-3:
            left_g.append(np.random.choice([-1,0,1],p=[1/3,1/3,1/3]))
        else:
            left_g.append(int(0))
    g_minor_and_left_g = np.concatenate((np.matrix(left_g).T,g_minor),axis=1)
    top_g = []
    for i in range(n):
        if i == 0:
            top_g.append(g_minor[n-3,n-3] * g_minor[n-2,n-2])
        else:
            top_g.append(int(0))
    g = np.concatenate((np.matrix(top_g), g_minor_and_left_g), axis=0)
    #end of making matrix g
    tester = None
    while tester is None:
        f_poly = random_poly(n-3, (n-3) // 3 + 1, neg_ones_diff=-1)
        try:
            f_p_poly = invert_poly(f_poly, R_poly, p)
            f_q_poly = invert_poly(f_poly, R_poly, q)
            tester =  f_p_poly * f_q_poly
        except NotInvertible as ex:
            log.debug(ex)
    if tester is None:
        raise Exception("Unable to generate polynomial f")
    deg_f = degree(f_poly, gen=x)
    f_arr =  f_poly.all_coeffs()[::-1]
    if deg_f != n-4:
        for i in range(n - 4 - deg_f):
            f_arr.append(int(0))
    little_f = shift_matriks(f_arr)
    f_2x2 = np.zeros((2,2))
    for i in range(2):
        for j in range(2):
            if i == 1 and j == 0:
                continue
            else:
                f_2x2[i,j] = np.random.choice([-1,1], p=[0.5,0.5])
    right_f = np.zeros((n-3,2))
    right_f_and_f2x2 = np.concatenate((right_f, f_2x2), axis=0)
    lower_f = np.zeros((2,n-3))
    lower_f_and_little_f =  np.concatenate((little_f, lower_f), axis=0)
    f_minor = np.concatenate((lower_f_and_little_f, right_f_and_f2x2), axis=1)
    left_f = []
    for i in range(n-1):
        left_f.append(np.random.choice([-1,0,1],p=[1/3,1/3,1/3]))
    f_minor_and_left_f = np.concatenate((np.matrix(left_f).T,f_minor),axis=1)
    top_f = []
    for i in range(n):
        if i == 0:
            top_f.append(f_minor[n-3,n-3] * f_minor[n-2,n-2])
        else:
            top_f.append(int(0))
    f = np.concatenate((np.matrix(top_f), f_minor_and_left_f), axis=0)
    #end of the making of f
    #f_p_poly and f_q_poly  alr made
    deg_f_p = degree(f_p_poly,gen=x)
    f_p_arr = f_p_poly.all_coeffs()[::-1]
    if deg_f_p != n-4:
        for i in range(n-4 - deg_f_p):
            f_p_arr.append(int(0))
    little_f_p = shift_matriks(f_p_arr)
    deg_f_q = degree(f_q_poly,gen=x)
    f_q_arr = f_q_poly.all_coeffs()[::-1]
    if deg_f_q != n-4:
        for i in range(n-4 - deg_f_q):
            f_q_arr.append(int(0))
    little_f_q = shift_matriks(f_q_arr)

    f_p_2x2 = np.remainder(np.linalg.inv(f_2x2), p)
    f_q_2x2 = np.remainder(np.linalg.inv(f_2x2), q)

    right_f_p = np.zeros((n-3,2))
    right_f_q = np.zeros((n-3,2))

    right_f_p_and_f_p_2x2 = np.concatenate((right_f_p, f_p_2x2), axis=0)
    right_f_q_and_f_q_2x2 = np.concatenate((right_f_q, f_q_2x2), axis=0)

    lower_f_p = np.zeros((2,n-3))
    lower_f_p_and_little_f_p =  np.concatenate((little_f_p, lower_f_p), axis=0)

    lower_f_q = np.zeros((2,n-3))
    lower_f_q_and_little_f_q =  np.concatenate((little_f_q, lower_f_q), axis=0)

    f_minor_p = np.concatenate((lower_f_p_and_little_f_p, right_f_p_and_f_p_2x2 ), axis=1)
    f_minor_q = np.concatenate((lower_f_q_and_little_f_q, right_f_q_and_f_q_2x2), axis=1)

    f_p = invMatriks(f_minor_p, f, p)
    f_q = invMatriks(f_minor_q, f,q)
    return f, g, f_p, f_q

def gen_keypub(keypriv, q):
    f_q = keypriv[3]
    g = keypriv[1]
    h = matriks_mul(f_q,g,q)
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
        c_temp = vektor_mul(pub, r, q)
        c_cipherteks = vektor_add(m[i], c_temp, q)
        c[i] = c_cipherteks
    c = np.ravel(c)
    return c

def dekripsi(c,  priv, size, n, p, q):
    c_row_size = len(c) // n
    c = c.reshape(c_row_size, n)
    d = np.copy(c)
    for i in range(c_row_size):  
        m_temp1 = vektor_mul(priv[0], c[i], q)
        m_plainteks = vektor_mul(priv[2],m_temp1,p)
        d[i] = m_plainteks
    d = np.ravel(d)
    if size != 0:
        d = d[slice(-size)]
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
    for i in range(m_add_size):
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