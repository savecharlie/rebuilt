"""The closed forms as printed in arXiv:2609.21772 (Prop 2.3 and Theorem 3.1),
evaluated in exact rationals, so a typo in a printed coefficient would show up
as a non-integer or as a disagreement with the paper's own Table."""
from fractions import Fraction as F
from math import factorial as f

def P3(n):  # Yu & Zhang, earlier paper, eq (3); n >= 2
    return (F(f(2*n-2), 3*f(n-1)*2**n)*(8*n**6+88*n**5+366*n**4+640*n**3+325*n**2-155*n-114)
            - F(1, 3)*f(n+1)*F(2)**(n-4)*(48*n**3+367*n**2+959*n+840))

def P4(n):  # Theorem 3.1; n >= 2
    a = F(f(2*n-4), 189*f(n-2)*2**(n-1))*(504*n**9+10836*n**8+91414*n**7+362607*n**6
        +568813*n**5-326256*n**4-1950557*n**3-1378566*n**2+440523*n+367416)
    b = F(1, 3)*f(n+1)*F(2)**(n-4)*(32*n**5+558*n**4+3901*n**3+13523*n**2+23008*n+15264)
    return a - b

TABLE4 = {1: 109, 2: 3881, 3: 113424, 4: 3212190, 5: 92486235, 6: 2759501745,
          7: 86034598650, 8: 2812696803000, 9: 96530594780025, 10: 3477318517324575}

if __name__ == '__main__':
    for n in range(2, 11):
        v = P4(n)
        print(n, v, 'integer' if v.denominator == 1 else 'NOT INTEGER',
              'matches table' if v == TABLE4[n] else 'DIFFERS from table')
