from sympy import solve 
from sympy import symbols
import numpy as np
import matplotlib.pyplot as plt
x_l, y_l, x_r, y_r = symbols('x_l y_l x_r y_r')
a, b, c, d = symbols('a b c d')

s = solve([-y_l+a*x_l**3 + b*x_l**2 + c*x_l +d, 
       -y_r+a*x_r**3 + b*x_r**2 + c*x_r +d,
       3*a*x_l**2 + 2*b*x_l + c,
       3*a*x_r**2 + 2*b*x_r + c], [a,b,c,d], dict=True)
coeffs = s[0]
print("a", coeffs[a])
print("b", coeffs[b])
print("c", coeffs[c])
print("d", coeffs[d])
def smoothstep(x__l,x__r,y__l,y__r, x):
    if x<1:
        return -1
    if x >= 1 and x < 2:
        
        A = coeffs[a].subs({x_l: x__l, x_r: x__r, y_l: y__l,y_r: y__r})
        B = coeffs[b].subs({x_l: x__l, x_r: x__r, y_l: y__l,y_r: y__r})
        C = coeffs[c].subs({x_l: x__l, x_r: x__r, y_l: y__l,y_r: y__r})
        D = coeffs[d].subs({x_l: x__l, x_r: x__r, y_l: y__l,y_r: y__r})
        P = A*x**3 + B*x**2 + C*x + D
        
        return P
     
    return 1
print(smoothstep(1,2,-1,1,1.5))
x = np.linspace(0,3,1000)
y = []
for i in range(len(x)):
    y.append(smoothstep(1,2,-1,1,x[i]))
plt.plot(x, y)
plt.show()