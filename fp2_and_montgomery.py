import random
import numpy
from math import sqrt
import math

def ModExp(x,y,n): 
  # Input: (x,y,n) with x,y,n integers and n>=2
  # Output: z = (x**y) mod n
  q = x
  z = 1
  while y>0:
    if y&1 == 1:
        z = (z*q) % n
    q = q*q % n
    y = y >> 1
  return z


def fpinv(a):
    # Input: prime residue class a
    # Output: the inverse of a mod 127
    x = ModExp(a,125,127) % 127
    return x
    

def fp2add(r,s):
  # Addition of Fp^2 Elements, with p =127
  # Elements are represented with x + iy, where x and y are in Fp
  # Add the numbers r = a + ib and s = c + id
    (a,b) = r
    (c,d) = s
    x = (a + c) % 127
    y = (b + d) % 127
    z = (x,y)
    return z

def fp2sub(r,s):
  # Subtraction of Fp^2 Elements, with p =127
  # Elements are represented with x + iy, where x and y are in Fp
  # Subtract the number c + id from a + ib
    (a,b) = r
    (c,d) = s
    x = (a - c) % 127
    y = (b - d) % 127
    z = (x,y)
    return z

def fp2mul(r,s):
  # Mulitplication of Fp^2 Elements, with p =127
  # Elements are represented with x + iy, where x and y are in Fp
  # Multiply the numbers a + ib and c + id
    (a,b) = r
    (c,d) = s
    x = (a*c - b*d) % 127
    y = (a*d + b*c) % 127
    z = (x,y)
    return z

def fp2inv(r):
  # Inversion of Fp^2 Element, with p =127
  # Elements are represented with x + iy, where x and y are in Fp
  # Find invers of the number a + ib
    (a,b) = r
    c = (a*a + b*b) % 127
    d = fpinv(c)
    x = (a*d) % 127
    y = (- b*d) % 127
    z = (x,y)
    return z

def fp2sqrt(s):
  # find the square root of a + ib in Fp^2
  # brute force, possible since p = 127
  # if a + ib is a quadratic non residue return 0
    for x in range(0,127):
        for y in range(0,127):
            z = (x,y)
            if fp2mul(z,z)==s:
                return z
    return(0,0)

def randompoint(b,a):
  # find random point on the elliptic curve E over Fp^2 given
  # in Montgomery form with coefficients b,a by
  # B y^2 = x^3 + A x^2 + x
    t = 0
    while t == 0:
        c = random.randint(1, 127)
        d = random.randint(1, 127)
        x = (c,d)
        x2 = fp2mul(x,x)
        x3 = fp2mul(x2,x)
        xx3 = fp2add(x,x3)
        ax2 = fp2mul(a,x2)
        f = fp2add(xx3,ax2)
        B1 = fp2inv(b)
        fB = fp2mul(f,B1)
        y = fp2sqrt(fB)
        if y != (0,0):
            t = 1
            return (f,x3,x,y)

   
def xdbl(X,Z,A,C):
  # optimized SIKE code
  # Input: projective XZ-coords of point P and projective curve coefficients A,C
  # Output: projective XZ-coords of [2]P
  # precompute coeff A24plus und C24
    C2 = fp2add(C,C)
    C24 = fp2add(C2,C2)
    A24plus = fp2add(A,C2)
  # doubling
    t0 = fp2sub(X, Z);                          # t0 = X1-Z1
    t1 = fp2add(X, Z);                          # t1 = X1+Z1
    t0 = fp2mul(t0, t0);                        # t0 = (X1-Z1)^2 
    t1 = fp2mul(t1, t1);                        # t1 = (X1+Z1)^2 
    Z2 = fp2mul(C24, t0);                       # Z2 = C24*(X1-Z1)^2   
    X2 = fp2mul(t1,Z2);                         # X2 = C24*(X1-Z1)^2*(X1+Z1)^2
    t1 = fp2sub(t1, t0);                        # t1 = (X1+Z1)^2-(X1-Z1)^2 
    t0 = fp2mul(A24plus, t1);                   # t0 = A24plus*[(X1+Z1)^2-(X1-Z1)^2]
    Z2 = fp2add(Z2, t0);                        # Z2 = A24plus*[(X1+Z1)^2-(X1-Z1)^2] + C24*(X1-Z1)^2
    Z2 = fp2mul(Z2, t1);                        # Z2 = [A24plus*[(X1+Z1)^2-(X1-Z1)^2] + C24*(X1-Z1)^2]*[(X1+Z1)^2-(X1-Z1)^2]
    return (X2,Z2)

def xdble(X,Z,A,C,e):
  # iterate xdbl e times to find [2^e]P
  # Input: projective XZ-coords of point P and projective curve coefficients A,C
  # Output: projective XZ-coords of [2^e]P
    X2e = X
    Z2e = Z
    for i in range(0,e):
        (X2e,Z2e) = xdbl(X2e,Z2e,A,C)
    return (X2e,Z2e)

def xtrl(X,Z,A,C):
  # optimized SIKE code
  # Input: projective XZ-coords of point P and projective curve coefficients A,C
  # Output: projective XZ-coords of [3]P
  # precompute coeff A24plus und A24minus
    C2 = fp2add(C,C)
    A24plus = fp2add(A,C2)
    A24minus = fp2sub(A,C2)
  # tripling   
    t0 = fp2sub(X, Z);                          # t0 = X-Z 
    t2 = fp2mul(t0, t0);                        # t2 = (X-Z)^2           
    t1 = fp2add(X, Z);                          # t1 = X+Z 
    t3 = fp2mul(t1, t1);                        # t3 = (X+Z)^2
    t4 = fp2add(t0, t1);                        # t4 = 2*X
    t0 = fp2sub(t1, t0);                        # t0 = 2*Z 
    t1 = fp2mul(t4, t4);                        # t1 = 4*X^2
    t1 = fp2sub(t1, t3);                        # t1 = 4*X^2 - (X+Z)^2 
    t1 = fp2sub(t1, t2);                        # t1 = 4*X^2 - (X+Z)^2 - (X-Z)^2
    t5 = fp2mul(t3, A24plus);                   # t5 = A24plus*(X+Z)^2 
    t3 = fp2mul(t3, t5);                        # t3 = A24plus*(X+Z)^3
    t6 = fp2mul(A24minus, t2);                  # t6 = A24minus*(X-Z)^2
    t2 = fp2mul(t2, t6);                        # t2 = A24minus*(X-Z)^3
    t3 = fp2sub(t2, t3);                        # t3 = A24minus*(X-Z)^3 - coeff*(X+Z)^3
    t2 = fp2sub(t5, t6);                        # t2 = A24plus*(X+Z)^2 - A24minus*(X-Z)^2
    t1 = fp2mul(t1, t2);                        # t1 = [4*X^2 - (X+Z)^2 - (X-Z)^2]*[A24plus*(X+Z)^2 - A24minus*(X-Z)^2]
    t2 = fp2add(t3, t1);                        # t2 = [4*X^2 - (X+Z)^2 - (X-Z)^2]*[A24plus*(X+Z)^2 - A24minus*(X-Z)^2] + A24minus*(X-Z)^3 - coeff*(X+Z)^3
    t2 = fp2mul(t2, t2);                        # t2 = t2^2
    X3 = fp2mul(t4, t2);                        # X3 = 2*X*t2
    t1 = fp2sub(t3, t1);                        # t1 = A24minus*(X-Z)^3 - A24plus*(X+Z)^3 - [4*X^2 - (X+Z)^2 - (X-Z)^2]*[A24plus*(X+Z)^2 - A24minus*(X-Z)^2]
    t1 = fp2mul(t1, t1);                        # t1 = t1^2
    Z3 = fp2mul(t0, t1);                        # Z3 = 2*Z*t1
    return (X3,Z3)

def xtrle(X,Z,A,C,e):
  # iterate xtrl e times to find [3^e]P
  # Input: projective XZ-coords of point P and projective curve coefficients A,C
  # Output: projective XZ-coords of [3^e]P
    X3e = X
    Z3e = Z
    for i in range(0,e):
        (X3e,Z3e) = xtrl(X3e,Z3e,A,C)
    return (X3e,Z3e)

def xDBLADD(XP,ZP,XQ,ZQ,xPQ,A,C):
  #  Simultaneous doubling and differential addition.
  #  Input: projective Montgomery points P=(XP:ZP) and Q=(XQ:ZQ) such that xP=XP/ZP and xQ=XQ/ZQ, affine difference xPQ=x(P-Q) and Montgomery curve constant A24=(A+2)/4.
  #  Output: projective Montgomery points P <- 2*P = (X2P:Z2P) such that x(2P)=X2P/Z2P, and Q <- P+Q = (XQP:ZQP) such that = x(Q+P)=XQP/ZQP. 
    #C2 = fp2add(C,C)
    #A24plus = fp2add(A,C2)
    #A24plus = fp2mul(A24plus,(32,0))
    c1 = fp2inv(C)
    a = fp2mul(c1,A)                  # unproj
    A24plus = fp2add(a,(2,0))         # minus 2
    A24plus = fp2mul(A24plus,(32,0))  # durch 4
    #C2 = fp2add(C,C)
    #A24plus = fp2add(A,C2)   
    t0 = fp2add(XP, ZP);                         # t0 = XP+ZP
    t1 = fp2sub(XP, ZP);                         # t1 = XP-ZP
    XP = fp2mul(t0,t0);                          # XP = (XP+ZP)^2
    t2 = fp2sub(XQ, ZQ);                         # t2 = XQ-ZQ
    #fp2correction(t2);
    XQ = fp2add(XQ, ZQ);                       # XQ = XQ+ZQ
    t0 = fp2mul(t0, t2);                       # t0 = (XP+ZP)*(XQ-ZQ)
    ZP = fp2mul(t1, t1);                       # ZP = (XP-ZP)^2
    t1 = fp2mul(t1, XQ);                       # t1 = (XP-ZP)*(XQ+ZQ)
    t2 = fp2sub(XP, ZP);                       # t2 = (XP+ZP)^2-(XP-ZP)^2
    XP = fp2mul(XP, ZP);                       # XP = (XP+ZP)^2*(XP-ZP)^2
    XQ = fp2mul(t2, A24plus);                      # XQ = A24*[(XP+ZP)^2-(XP-ZP)^2]
    ZQ = fp2sub(t0, t1);                       # ZQ = (XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)
    ZP = fp2add(XQ, ZP);                       # ZP = A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2
    XQ = fp2add(t0, t1);                       # XQ = (XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)
    ZP = fp2mul(ZP, t2);                       # ZP = [A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2]*[(XP+ZP)^2-(XP-ZP)^2]
    ZQ = fp2mul(ZQ, ZQ);                       # ZQ = [(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2
    XQ = fp2mul(XQ, XQ);                       # XQ = [(XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)]^2
    ZQ = fp2mul(ZQ, xPQ);                      # ZQ = xPQ*[(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2
# XP and ZP are projective coords of 2P
# XQ and ZQ are projective coords of P+Q
    return (XP,ZP,XQ,ZQ)

def xsev(X,Z,A,C):
  # computing the [7]-multiple of a point with xDBLADD
  # Input: projective XZ-coords of point P and projective curve coefficients A,C
  # Output: projective XZ-coords of [7]P
    (X2,Z2) = xdbl(X,Z,A,C)
    z1 = fp2inv(Z)
    x = fp2mul(X,z1)
    (X4,Z4,X3,Z3) = xDBLADD(X2,Z2,X,Z,x,A,C)
    (X8,Z8,X7,Z7) = xDBLADD(X4,Z4,X3,Z3,x,A,C)
    return (X7,Z7)

def get3isog(XP,ZP):
  # optimized SIKE code
  # Computes the corresponding 3-isogeny of a projective Montgomery point (X3:Z3) of order 3.
  # Input:  projective point of order three P = (X3:Z3).
  # Output: the 3-isogenous Montgomery curve with projective coefficient A/C. 
    co0 = (0,0)
    co1 = (0,0)
    co0 = fp2sub(XP,ZP);                        # coeff0 = X-Z
    t0 = fp2mul(co0, co0);                      # t0 = (X-Z)^2
    co1 = fp2add(XP, ZP);                       # coeff1 = X+Z
    t1 = fp2mul(co1, co1);                      # t1 = (X+Z)^2
    t2 = fp2add(t0, t1);                        # t2 = (X+Z)^2 + (X-Z)^2
    t3 = fp2add(co0, co1);                      # t3 = 2*X
    t3 = fp2mul(t3, t3);                        # t3 = 4*X^2
    t3 = fp2sub(t3, t2);                        # t3 = 4*X^2 - (X+Z)^2 - (X-Z)^2 
    t2 = fp2add(t1, t3);                        # t2 = 4*X^2 - (X-Z)^2 
    t3 = fp2add(t3, t0);                        # t3 = 4*X^2 - (X+Z)^2
    t4 = fp2add(t0, t3);                        # t4 = 4*X^2 - (X+Z)^2 + (X-Z)^2 
    t4 = fp2add(t4, t4);                        # t4 = 2(4*X^2 - (X+Z)^2 + (X-Z)^2) 
    t4 = fp2add(t1, t4);                        # t4 = 8*X^2 - (X+Z)^2 + 2*(X-Z)^2
    A24minus = fp2mul(t2, t4);                  # A24minus = [4*X^2 - (X-Z)^2]*[8*X^2 - (X+Z)^2 + 2*(X-Z)^2]
    t4 = fp2add(t1, t2);                        # t4 = 4*X^2 + (X+Z)^2 - (X-Z)^2
    t4 = fp2add(t4, t4);                        # t4 = 2(4*X^2 + (X+Z)^2 - (X-Z)^2) 
    t4 = fp2add(t0, t4);                        # t4 = 8*X^2 + 2*(X+Z)^2 - (X-Z)^2
    A24plus = fp2mul(t3, t4);                   # A24plus = [4*X^2 - (X+Z)^2]*[8*X^2 + 2*(X+Z)^2 - (X-Z)^2]
    # modify
    A2 = fp2add(A24plus,A24minus)
    C4 = fp2sub(A24plus,A24minus)
    i = fp2inv((4,0))
    k = fp2inv((2,0))
    A = fp2mul(A2,k)
    C = fp2mul(C4,i)
    j = j_inv(A,C)
    return (j,A,C,co0,co1)
    #return (A24minus,A24plus,co0,co1)

def eval3isog(XQ,ZQ,co0,co1):
  # Computes the 3-isogeny R=phi(X:Z), given projective point (X3:Z3) of order 3 on a Montgomery curve and 
  # a point P with 2 coefficients in coeff (computed in the function get_3_isog()).
  # Inputs: projective points P = (X3:Z3) and Q = (X:Z).
  # Output: the projective point Q <- phi(Q) = (X3:Z3). 
    t0 = fp2add(XQ, ZQ);                     # t0 = X+Z
    t1 = fp2sub(XQ, ZQ);                     # t1 = X-Z
    t0 = fp2mul(t0, co0);                    # t0 = coeff0*(X+Z)
    t1 = fp2mul(t1, co1);                    # t1 = coeff1*(X-Z)
    t2 = fp2add(t0, t1);                     # t2 = coeff0*(X+Z) + coeff1*(X-Z)
    t0 = fp2sub(t1, t0);                     # t0 = coeff1*(X-Z) - coeff0*(X+Z)
    t2 = fp2mul(t2, t2);                     # t2 = [coeff0*(X+Z) + coeff1*(X-Z)]^2
    t0 = fp2mul(t0, t0);                     # t0 = [coeff1*(X-Z) - coeff0*(X+Z)]^2
    phiXQ = fp2mul(XQ, t2);                  # X3final = X*[coeff0*(X+Z) + coeff1*(X-Z)]^2        
    phiZQ = fp2mul(ZQ, t0);                  # Z3final = Z*[coeff1*(X-Z) - coeff0*(X+Z)]^2
    return (phiXQ,phiZQ)

def get4isog(XP,ZP):
  # Computes the corresponding 4-isogeny of a projective Montgomery point (X4:Z4) of order 4.
  # Input:  projective point of order four P = (X4:Z4).
  # Output: the 4-isogenous Montgomery curve with projective coefficients A+2C/4C and the 3 coefficients 
  #         that are used to evaluate the isogeny at a point in eval_4_isog().
    
    co1 = fp2sub(XP, ZP);                   # coeff[1] = X4-Z4
    co2 = fp2add(XP, ZP);                   # coeff[2] = X4+Z4
    co0 = fp2mul(ZP, ZP);                   # coeff[0] = Z4^2
    co0 = fp2add(co0, co0);                 # coeff[0] = 2*Z4^2
    C24 = fp2mul(co0, co0);                 # C24 = 4*Z4^4
    co0 = fp2add(co0, co0);                 # coeff[0] = 4*Z4^2
    A24plus = fp2mul(XP,XP);                # A24plus = X4^2
    A24plus = fp2add(A24plus, A24plus);     # A24plus = 2*X4^2
    A24plus = fp2mul(A24plus, A24plus);     # A24plus = 4*X4^4
  # modify
    i = fp2inv((4,0))
    C = fp2mul(C24,i)
    A = fp2sub(A24plus,C)
    A = fp2sub(A,C)
    j = j_inv(A,C)
    return (j,A,C,co0,co1,co2)
    #return(A24plus,C24,co0,co1,co2)


def eval4isog(XP,ZP,co0,co1,co2):
# Evaluates the isogeny at the point (X:Z) in the domain of the isogeny, given a 4-isogeny phi defined 
# by the 3 coefficients in coeff (computed in the function get_4_isog()).
# Inputs: the coefficients defining the isogeny, and the projective point P = (X:Z).
# Output: the projective point P = phi(P) = (X:Z) in the codomain. 
    
    t0 = fp2add(XP, ZP);                    # t0 = X+Z
    t1 = fp2sub(XP, ZP);                    # t1 = X-Z
    phiXP = fp2mul(t0, co1);                # X = (X+Z)*coeff[1]
    phiZP = fp2mul(t1, co2);                # Z = (X-Z)*coeff[2]
    t0 = fp2mul(t0, t1);                    # t0 = (X+Z)*(X-Z)
    t0 = fp2mul(t0, co0);                   # t0 = coeff[0]*(X+Z)*(X-Z)
    t1 = fp2add(phiXP, phiZP);              # t1 = (X-Z)*coeff[2] + (X+Z)*coeff[1]
    phiZP = fp2sub(phiXP, phiZP);           # Z = (X-Z)*coeff[2] - (X+Z)*coeff[1]
    t1 = fp2mul(t1, t1);                    # t1 = [(X-Z)*coeff[2] + (X+Z)*coeff[1]]^2
    phiZP = fp2mul(phiZP,phiZP);            # Z = [(X-Z)*coeff[2] - (X+Z)*coeff[1]]^2
    phiXP = fp2add(t1, t0);                 # X = coeff[0]*(X+Z)*(X-Z) + [(X-Z)*coeff[2] + (X+Z)*coeff[1]]^2
    t0 = fp2sub(phiZP, t0);                  # t0 = [(X-Z)*coeff[2] - (X+Z)*coeff[1]]^2 - coeff[0]*(X+Z)*(X-Z)
    phiXP = fp2mul(phiXP, t1);               # Xfinal
    phiZP = fp2mul(phiZP, t0);               # Zfinal
    return (phiXP,phiZP)

def kernel7(X,Z,A,C):
  # for degree 7 Iso we need the coords of P,2P and 3P
  (X2,Z2) = xdbl(X,Z,A,C)
  (X3,Z3) = xtrl(X,Z,A,C)
  return (X,Z,X2,Z2,X3,Z3)

def crisscross(a,b,c,d):
  # write criss cross function as suggested in Costello & Hisil
  t1 = fp2mul(a, d);
  t2 = fp2mul(b, c);
  return (fp2add(t1,t2),fp2sub(t1,t2))

def get7isog(X1,Z1,X2,Z2,X3,Z3,A,C):
  # calculate the Montgomery curve coefficient with Theorem 1 in Costello & Hisil (only okay if degree of isogeny is very low)
  # get affine x coord
  t0 = fp2inv(Z1)
  x1 = fp2mul(X1,t0)
  ix1 = fp2inv(x1)
  t0 = fp2inv(Z2)
  x2 = fp2mul(X2,t0)
  ix2 = fp2inv(x2)
  t0 = fp2inv(Z3)
  x3 = fp2mul(X3,t0)
  ix3 = fp2inv(x3)
  # get affine curve
  t0 = fp2inv(C)
  a = fp2mul(A,t0)
  # sigma, sigmat, pi
  t1 = fp2add(x1,x2)
  sigma = fp2add(x3,t1)
  t1 = fp2add(ix1,ix2)
  sigmat = fp2add(ix3,t1)
  t1 = fp2mul(x1,x2)
  t1 = fp2mul(t1,x3)
  pisq = fp2mul(t1,t1)
  # get a'
  t2 = fp2sub(sigmat,sigma)
  t3 = fp2mul(t2,(6,0))
  t4 = fp2add(t3,a)
  a = fp2mul(t4,pisq)
  return a
  
def eval7isog(X1,Z1,X2,Z2,X3,Z3,XQ,ZQ):
  # Input points are the three kernel elements P,2P,3P and the point Q that needs to be evaluated under the 7 isogeny
  # with kernel generated by P
  # the calculation follows the general formulae from Costello & Hisil in Chapter 5
  # uses crisscross for simplicity and needs the proj-coords of the P,2P and 3P
  # preparation for crisscross
  hXQ = fp2add(XQ , ZQ) 
  hZQ = fp2sub(XQ , ZQ)
  hX1 = fp2add(X1 , Z1) 
  hZ1 = fp2sub(X1 , Z1)
  hX2 = fp2add(X2 , Z2) 
  hZ2 = fp2sub(X2 , Z2)
  hX3 = fp2add(X3 , Z3) 
  hZ3 = fp2sub(X3 , Z3)
  # start with P1
  (phiX,phiZ) = crisscross(hX1,hZ1,hXQ,hZQ) # Algorithm 1
  # multiplication for P2
  (t0, t1) = crisscross(hX2,hZ2,hXQ,hZQ) # Algorithm 1
  phiX =  fp2mul(t0,phiX)
  phiZ =  fp2mul(t1,phiZ)
  # multiplication for P3
  (t0, t1) = crisscross(hX3,hZ3,hXQ,hZQ) # Algorithm 1
  phiX =  fp2mul(t0,phiX)
  phiZ =  fp2mul(t1,phiZ)
  # square
  phiX = fp2mul(phiX,phiX) # 2S
  phiZ = fp2mul(phiZ,phiZ)
  # multiply with initial coords
  phiX = fp2mul(XQ ,phiX) # 2M
  phiZ = fp2mul(ZQ ,phiZ)
  return (phiX,phiZ)

  
def j_inv(A, C):
# Computes the j-invariant of a Montgomery curve with projective constant.
# Input: A,C in GF(p^2).
# Output: j=256*(A^2-3*C^2)^3/(C^4*(A^2-4*C^2)), which is the j-invariant of the Montgomery curve B*y^2=x^3+(A/C)*x^2+x
# or (equivalently) j-invariant of B'*y^2=C*x^3+A*x^2+C*x.
    jinv = fp2mul(A,A);                        # jinv = A^2        
    t1 = fp2mul(C,C);                          # t1 = C^2
    t0 = fp2add(t1, t1);                       # t0 = t1+t1 = 2C^2
    t0 = fp2sub(jinv, t0);                     # t0 = jinv-t0 = A^2 -2C^2
    t0 = fp2sub(t0, t1);                       # t0 = t0-t1 = A^2 -3C^2
    jinv = fp2sub(t0, t1);                     # jinv = t0-t1 = A^2 - 4C^2
    t1 = fp2mul(t1, t1);                       # t1 = t1^2 = C^4
    jinv = fp2mul(jinv, t1);                   # jinv = jinv*t1 = (A^2 - 4C^2)C^4
    t0 = fp2add(t0, t0);                       # t0 = t0+t0 
    t0 = fp2add(t0, t0);                       # t0 = t0+t0 = 4(A^2 -3C^2)
    t1 = fp2mul(t0, t0);                       # t1 = t0^2 = 16(A^2 -3C^2)^2
    t0 = fp2mul(t0, t1);                       # t0 = t0*t1 = 64(A^2 -3C^2)^3
    t0 = fp2add(t0, t0);                       # t0 = t0+t0
    t0 = fp2add(t0, t0);                       # t0 = t0+t0 = 256(A^2 -3C^2)^3
    jinv = fp2inv(jinv);                       # jinv = 1/jinv      # equals the formulae we gave in chapter 3.2 for proj Montgomery curves
    jinv = fp2mul(jinv, t0);                   # jinv = t0*jinv
    return jinv

def affinsub(x1,y1,x2,y2,B):
  # Input: two point in affine coordinates on a Montgomery curve
  #        P = (x1,y1) Q =(x2,y2)
  # Output: the x-coordinate of the difference P-Q = (x4,y4)
  # formula from Montgomery "Speeding the Pollard"
    t0 = fp2sub(x1,x2)                      # t0 = x1 - x2
    t0 = fp2mul(t0,t0)                      # t0 = (x1 - x2)^2
    t1 = fp2mul(x2,y1)                      # t1 = x2 * y1
    t2 = fp2mul(x1,y2)
    t2 = fp2add(t1,t2)                      # t2 = (x2y1 + x1y2)
    t2 = fp2mul(t2,t2)                      # t2 = (x2y1 + x1y2)^2
    t3 = fp2mul(x1,x2)
    t3 = fp2mul(t0,t3)                      # t3 = (x1 - x2)^2 x1x2
    t3 = fp2inv(t3)                         # inv
    x4 = fp2mul(t2,B)
    x4 = fp2mul(t3,x4)                      # x4 = B(x2y1 + x1y2)^2 / ((x1 - x2)^2 x1x2)
    return x4

#### achtung in 3 point ladder nur affine coord einsetzen momentan
def threeptladder(XP,ZP,XQ,ZQ,XPQ,t,A,C):
  # 3 point ladder as in De Feo, Jao, Plut
  # Input: P,Q and x-coord of P-Q on the elliptic curve with proj curve coeff A,C
  # Output: P + [t]Q
  l = len(bin(t))-2
  AX = (0,0)
  AZ = (0,0)
  BX = XQ
  BZ = ZQ
  CX = XP
  CZ = ZP
  for i in range(0,l):
    x = t >> l-1-i
    if x&1 == 0:
      print(0)
      #(AX,AZ) = xDBL(AX,AZ,A,C)
      (a,b,BX,BZ) = xDBLADD(AX,AZ,BX,BZ,XQ,A,C)
      (AX,AZ,CX,CZ) = xDBLADD(AX,AZ,CX,CZ,XP,A,C)
    else:
      print(1)
      if AZ == (0,0):
        (AX,AZ) = (BX,BZ)
      else:
        (a,b,AX,AZ) = xDBLADD(AX,AZ,BX,BZ,XQ,A,C)
      (BX,BZ,CX,CZ) = xDBLADD(BX,BZ,CX,CZ,XPQ,A,C)
  return (CX,CZ)

def recova(XP,XQ,XPQ):
# formula from Costello & Hisil
# recover the curve coefficient of the image curve from the x-coordinates of the points P, Q and P-Q
# Chapter 4 equation (16)
  t0 = fp2mul(XP,XQ)
  t1 = fp2mul(XP,XPQ)
  t2 = fp2mul(XQ,XPQ)
  t3 = fp2mul(t0,XPQ)
  a = fp2sub((1,0),t0)
  a = fp2sub(a,t1)
  a = fp2sub(a,t2)
  a = fp2mul(a,a)
  x = fp2mul((4,0),t3)
  x = fp2inv(x)
  a = fp2mul(a,x)
  a = fp2sub(a,XP)
  a = fp2sub(a,XQ)
  a = fp2sub(a,XPQ)
  return a


def norm(X,Z):
  # normalize projective coordinates X,Z to x=X/Z or projective curve coefficients A,C to a=A/C
  a = fp2inv(Z)
  x = fp2mul(X,a)
  return x


