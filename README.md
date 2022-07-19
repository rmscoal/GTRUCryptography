# GTRUCryptography

## What is this project
This project is a work that I'm doing on my **undergraduate thesis**. The title of my undergraduate thesis "*GTRU Public-key Cryptosystem and its Application using $\mathbb{Z}^n$ and Poly-$\mathbb{Z}$ Group*". The purpose of this project is to provide the Python implementation of the discussed cryptosystem. This program is to simulate key building, encryption and decryption process of the cryptosystem. In the thesis, I strictly discussed about the performance GTRU has and later on be compared with RSA public-key cryptosystem. Moreover, I also discussed on how efficient is this cryptosystem implementation.

## What is Cryptography
To put it short, cryptography is a field in mathematics that study on how to secure sending message over the noisy and dirtyh internet. Cryptography is based on mathematical problems such as the difficulty on factorizing large numbers, shortest vector problem, and many more. Cryptosystem is an algorithm made from simple cryptography. A cryptosystem can be said as the by product of cryptography. In the real world,  a cryptosystem should be able to  prioritize its security and it is hard to break, private keys is hard to find or plaintext is hardly recover from ciphertext. On top of that, a cryptosystem should also be able to deliver a great performance, as in speed, and is efficient. If so, then that cryptosystem is ideal.

## What is GTRU Public-Key Cryptosystem
Before going further, let me explain about NTRU public-key cryptosystem. 

### NTRU Public-Key Cryptosystem
NTRU (*Nth degree Truncated polynomial Ring Unit*) is a public-key cryptosystem working on the ring $\mathcal{R} = \mathbb{Z}[x]/\langle x^N - 1 \rangle$. The operations are addition and multiplication (or cyclic convolution product) in the ring. Another operation on NTRU is the modulo reduction of $p$ and $q$ where $\gcd(p,q) = 1$ and $p << q$. NTRU was created by Hoffstein, Jeffrey, and Silverman in 1996 and was a lattice based problems. It is a cryptosystem post-quantum era. As Hoffstein claimed, NTRU features reasonably short, easily created keys, high speed, and low memory requirements.

### GTRU Public-Key Cryptosystem
GTRU (*a Group-based NTRU-like*) public-key cryptosystem is a cryptosystem similar to NTRU, as its name suggests. GTRU is made by generalizing NTRU into the group structure. The idea comes from looking cyclic convolution products as a vector-matrix product in the group $\mathbb{Z}^n$. With this then, it is generalized further where the vector-matrix product is seen as the image of an endomorphism of $\mathbb{Z}^n$. With this moving forward, NTRU is brought into the group structure by using normal subgroups $P$ and $Q$ in $G$ and natural homomorphism 
$$\pi_P : G \rightarrow G/P$$ 
as the replacement for modulo reduction of $p$ and 
$$\pi_Q : G \rightarrow G/Q$$ 
for modulo reduction of $q$. Furthermore, GTRU make uses of the functions 
$$\pi_N : G \rightarrow G/N, \pi_N(g) = gN$$
$$\rho_{T_N} : G/N \rightarrow G, \rho_{T_N}(gN) = g_{T_N} \in gN \cap T_N,$$
$$\bar{\pi}_N : E(G)_N \rightarrow E(G/N), \bar{\pi}_N(\mathfrak{f})(gN) = \mathfrak{f}(g)N.$$

The complete algorithm can be seen as follows.

**Public Parameter:** Given a group $G$, normal subgroups $P$ and $Q$ in $G$, $T_P$ transversal of $P$ and $T_Q$ transversal of $Q$, the sets $\mathcal{L}_f, \mathcal{L}_g \subseteq E(G)$, and $\mathcal{L}_m, \mathcal{L}_r \subseteq G$.

**Key Generation:** Choose $\mathfrak{f} \in \mathcal{L}_f$ and $\mathfrak{g} \in \mathcal{L}_g$ such that there exist $\mathfrak{f}_P$ and $\mathfrak{f}_Q$ that follows 
$$\bar{\pi}_P(\mathfrak{f}_P \circ \mathfrak{f}) \circ \pi_P = \pi_P,$$ 
$$\bar{\pi}_Q(\mathfrak{f} \circ \mathfrak{f}_Q) \circ \pi_Q = \pi_Q.$$ 
Calculate $\mathfrak{h} = \bar{\pi}_Q(\mathfrak{f}_Q \circ \mathfrak{g})$. The public key is $\mathfrak{h}$ and the private key is the pair $(\mathfrak{f}, \mathfrak{f}_P)$.

**Encryption:** To encrypt a message $m \in \mathcal{L}_m$, randomly choose $r \in \mathcal{L}_r$ and calculate 
$$c=\pi_Q(m) \cdot \mathfrak{h} \circ \pi_Q(r)$$

**Decryption:** To decrypt the ciphertext $c \in G/Q$, calculate 
$$m = \rho_{T_P} \circ \bar{\pi}_P(\mathfrak{f}_P) \circ \pi_P \circ \rho_{T_Q} \circ \bar{\pi}_Q(f)(c).$$ 

Now, not all groups can be used for GTRU cryptosystem. There are surely a number of requirements a group must have. One example is a group $G$ must at least have two non-trivial normal subgroups. With all this said, this program is a GTRU simulation program using the groups $\mathbb{Z}^n$ under vector addition and a special poly-$\mathbb{Z}$ group 
$$G_n = \mathbb{Z}^{n-1} \rtimes_\phi \mathbb{Z}$$ 
where $\phi: \mathbb{Z} \rightarrow Aut(\mathbb{Z}^{n-1})$ defined by 
$$
\phi (a) = \begin{pmatrix}
1 & 0 &\dots& 0 \\
0 & 1 & \dots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
a & 0 & \dots & 1 
\end{pmatrix}, a \in \mathbb{Z}.
$$

## Repository Explanation
The folder ***GTRU-G*** contains the program for GTRU using the group $\mathbb{Z}^n$. Meanwhile, the folder ***GTRU-Gn*** contains the program for GTRU using the poly-$\mathbb{Z}$ group $G_n$ as mentioned. In each folder, you will see a brief explanationn about what the codes are doing.