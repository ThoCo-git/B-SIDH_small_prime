## B-SIDH_small_prime
This is an assembly of functions to perform field arithmetic in a finite field F_(p^2), elliptic curve arithmetic and degree 3, 4, and 7 isogeny computation in projective XZ Montgomery coordinates.
It is implemented for p=127, but can easily be adapted to other small prime numbers congruent to 3 mod 4.
It provides all necessary functions to perform a walk in the supersingular isogeny graph, as described in Costello's work "B-SIDH: supersingular isogeny Diffie-Hellman using twisted torsion"

### F_(p^2) field arithmetic
We choose F_(p^2) = F_p(i), so the quadratic extension of F_p is realized by adjoing a square root of -1.
The implemented arithmetic operations are by no means optimal, the focus was to keep the code as short and simple as possible.
The purpose of the code is to give non-experts in this field a help to comprehend the basic steps of the B-SIDH protocol, so we assume that only small primes are considered, like p=127 in this example.
Functions:  - fp2add
            - fp2sub
            - fp2mul
            - fp2inv (uses fpinv and ModExp)
            - fp2sqrt (brute-force, only applicable for small primes)

### Elliptic curve arithmetic
The point operations that are necessary for the protocol depend on the factors of p+1 and p-1. In our example p=127 doublings, triplings and multiplication by 7 maps are sufficient. The functions work with projectivized XZ Montgomery arithmetic.
Functions:  - xdbl, xdble ([2] and [2^e], taken from the Microsoft SIDH library)
            - xtrl, xtrle ([3] and [3^e], taken from the Microsoft SIDH library)
            - xsev ([7], selfmade by using xDBLADD)
            - xDBLADD (double and add, taken from the Microsoft SIDH library)
            - randompoint (get a random point on the elliptic curve E with Montgomery coefficients b,a)
            - j_inv (get the j-invariant of a curve from the projective curve coefficients A and C)

### Isogeny computation
Similar to the elliptic curve arithmetic, the necessary isogeny equations depend on p+1 and p-1. Since the only only factors are 2,3 and 7 we only need isogeny formulas of degree 2,3 and 7. To save some iterations we replace degree 2 by degree 4 isogenies.
Functions:  - get3isog (get image curve, taken from the Microsoft SIDH library)
            - eval3isog (evaluate point, taken from the Microsoft SIDH library)
            - get4isog (get image curve, taken from the Microsoft SIDH library)
            - eval4isog (evaluate point, taken from the Microsoft SIDH library)
            - get7isog (get image curve, pseudo code from Costello & Hisil)
            - eval7isog (evaluate point, pseudo code from Costello & Hisil)
