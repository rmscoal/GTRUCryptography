'''
___________________________________________

Aplikasi GTRU Menggunakan Z^n Versi 2
___________________________________________


____________________________________________

In this version, we use the subset E_t(G)
which is the lower triangular matrix of
entries in the set {-1,0,1} and with no 0's
in the main diagonal. Previous matrix inversion
is done by using the numpy.linalg.inv method.
In this version, we use the Sympy library
to calculate the matrix inversion. After some
testing, the numpy version has better speed yet
cannot afford large size matrices. On the other
hand, sympy is able to deliver inversion of large
size matrices yet slower speed is achieved.

Using this method, the largest size matrices that
have been tried is 251 x 251, which is the recom-
mended size for NTRU.
_____________________________________________

'''


#import modul
from sympy.matrices import Matrix
from sympy import *
import numpy as np
import math
import time
import random

'''
PEMBANGUNAN MATRIKS, VEKTOR, DAN OPERASI MOD DI MATRIKS
'''


def lower_matriks(n):
    f = Matrix.zeros(n)
    for i in range(n):
        for j in range(n):
            if i == j:
                f[i, j] = np.random.choice([-1, 1], p=[1/2, 1/2])
            elif i > j:
                f[i, j] = np.random.choice([-1, 0, 1], p=[1/3, 1/3, 1/3])
    return f


def randomrgenerator(n, p):
  arr_r = []
  l = [-p, 0, p]
  for i in range(n):
    x = random.choice(l)
    arr_r.append(x)

  return Matrix(1, n, np.array(arr_r))


def matrikskalimod(f, g, b):
    result = Matrix(Mod(g @ f, b))
    n = shape(result)[0]
    for i in range(n):
        for j in range(n):
            if result[i, j] > b/2:
                result[i, j] -= b
    return result


def vektorkalimod(f, x, b):
    result = Matrix(Mod(x @ f, b))
    n = shape(result)[1]
    for i in range(n):
        if result[0, 1] > b/2:
            result[0, i] -= b
    return result


def mod(mat, q):
    n, m = shape(mat)
    mat = Matrix(mat)
    for i in range(n):
        for j in range(m):
            mat[i, j] = mat[i, j] % q
            if mat[i, j] > q/2:
                mat[i, j] -= q
    return mat


'''
PEMBANGUNAN KUNCI
'''


def poly_r(n, p):
    r = randomrgenerator(n, p)
    return r


def transversalchecker(f, g, m, r, q):
    size = shape(m)[1]
    for i in range(size):
        result_1 = Matrix(m.row(i) * f)
        result_2 = Matrix(r * g)
        result = Matrix(result_1 + result_2)
        n = shape(result)[1]
        for i in range(n):
            if result[0, i] < -(q)/2 or result[0, i] > q/2:
              return False
              break
            else:
                return True


def gen_keypriv(n, p, q):
    f = lower_matriks(n)
    g = lower_matriks(n)
    f_q = Mod(f.inv(), q)
    f_p = Mod(f.inv(), p)

    return f, g, f_p, f_q


def gen_keypub(keypriv):
    f_q = keypriv[3]
    g = keypriv[1]
    h = matrikskalimod(f_q, g, q)

    return h


def enkripsi(pub, m, r, q):
    size = shape(m)[0]
    n = shape(m)[1]
    c = Matrix(m.copy())
    for i in range(size):
        c_temp = r * pub
        c_cipherteks = m.row(i) + c_temp
        c = Matrix(c.row_insert(i, c_cipherteks))
    return c[:size, :].reshape(1, len(c[:size, :]))


def dekripsi(c,  priv, size, n):
    c_row_size = len(c) // n
    c = c.reshape(c_row_size, n)
    d = c.copy()
    for i in range(c_row_size):
        m_temp1 = c.row(i) * priv[0]
        m_plainteks = Mod(m_temp1 * priv[2], p)
        d = Matrix(d.row_insert(i, m_plainteks))
    des = np.array(d[:c_row_size, :].reshape(1, len(d[:c_row_size, :])))
    des = des[0, slice(-size)]
    des = bintostr(des)
    return des


'''
PENYESUAIAN MESSAGE
'''


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
    return ascii_text.replace('\x00', '')


def custom_message_size(m, n):
    m_size_temp = len(m) % n
    m_add_size = n - m_size_temp
    for i in range(m_add_size):
        m.append(0)
    m_row_size = len(m) // n
    m_arr = Matrix(m_row_size, n, np.array(m))
    return m_arr, m_add_size, m_row_size


'''
JALANKAN PROGRAM
'''

n = 251
p = 3
q = 2048

start_time = time.time()
keypriv = gen_keypriv(n, p, q)
r = poly_r(n, p)

message = '''Lorem ipsum dolor sit amet,
consectetur adipiscing elit,
sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.
Ut enim ad minim veniam,
quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo
consequat.
Duis aute irure dolor in reprehenderit in voluptate velit esse cillum
dolore eu fugiat nulla pariatur.
Excepteur sint occaecat cupidatat non proident,
sunt in culpa qui officia deserunt mollit anim id est laborum.
Saya disini adalah seorang mahasiswa di Universitas Gadjah Mada jurusan Matematika!
'''
m = custom_message_size(strToBinary(message), n)[0]
m_add_size = custom_message_size(strToBinary(message), n)[1]


f, g = keypriv[0], keypriv[1]
tries = 5
while tries > 0 and transversalchecker(f, g, m, r, q) == False:
    f, g = gen_keypriv(n, p, q)[0], gen_keypriv(n, p, q)[1]
    tries -= 1
    print(f"after {5-tries} the key is built")

keypub = gen_keypub(keypriv)

start_time_enkripsi = time.time()
e = enkripsi(keypub, m, r, q)
print(f"cipherteks: \n {mod(e,q)}")
print("--- waktu enkripsi : %s milidetik ---\n" %
      ((time.time() - start_time_enkripsi)*1000))

start_time_dekripsi = time.time()
d = dekripsi(e,  keypriv, m_add_size, n)
print(f"hasil dekripsi: {d}")
print("--- waktu dekripsi : %s milidetik ---" %
      ((time.time() - start_time_dekripsi)*1000))

print("\n--- waktu total (temasuk pembangunan kunci) : %s seconds ---" %
      (time.time() - start_time))
