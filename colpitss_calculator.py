from math import pi
import cmath
import math
from tabulate import tabulate
from sympy import symbols, Eq, solve
from si_prefix import si_format

VO = 1.0
beta = 1/10
Vf = beta * VO

F = 1E6
w = 2*pi*F
L = 100E-6
XL = w*L
ZL = XL*1j
Ct = 1/(w*w*L)

RL = round(abs(ZL),3)
x, y = symbols('c1 c2')
eq1 = Eq((beta*x) - y )
eq2 = Eq((Ct*x) + (Ct*y) + (x*y) )
temp = solve((eq1,eq2), (x, y)) 
c1,c2 = temp[0]
c1 = float(abs(c1))
c2 = float(abs(c2))
ZC1 = 1/(c1*w*1j)
ZC2 = 1/(c2*w*1j)

RC1 = abs(ZC1)
RC2 = abs(ZC2)



zf = cmath.polar(1/((1/ZC1)+(1/(ZL+ZC2))))
zO = cmath.polar(1/((1/ZC2)+(1/(ZL+ZC1))))

Zf = (si_format(round(zf[0]),precision=2), math.degrees(zf[1]))
ZO = (si_format(round(zO[0]),precision=2), math.degrees(zO[1]))

Rf = 1/((1/RC1)+(1/(RL+RC2)))
RO = 1/((1/RC2)+(1/(RL+RC1)))

ZT = ZL + ZC1 + ZC2

If = abs(Vf/Rf)
IO = abs(VO/RO)

print(tabulate([['Beta', beta],
                ['Freq', si_format(F,precision=2)],
                ['Omega',w],
                ['L',si_format(L,precision=2)],
                ['C1',si_format(c1,precision=2)],
                ['C2',si_format(c2,precision=2)],
                ['Ct',si_format(Ct,precision=2)],
                ['ZL',cmath.polar(ZL)],
                ['ZC1',cmath.polar(ZC1)],
                ['ZC2',cmath.polar(ZC2)],
                ['ZT',cmath.polar(ZT)],
                ['Zf',Zf],
                ['ZO',ZO],
                ['RL',RL],
                ['RC1',si_format(RC1,precision=2)],
                ['RC2',si_format(RC2,precision=2)],
                ['Rf',si_format(Rf,precision=2)],
                ['RO',si_format(RO,precision=2)],
                ['Vf',si_format(Vf,precision=2)],
                ['VO',si_format(VO,precision=2)],
                ['If',si_format(If,precision=2)],
                ['IO',si_format(IO,precision=2)]],
                headers=['Components', 'Value']))

