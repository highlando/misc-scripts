# this script symbolically computes the P2 basis 
# on an arbitrary triangle
# following the approach and notation of 
# http://www.mathematik.uni-dortmund.de/~kuzmin/cfdintro/lecture7.pdf

# then one can set the edge points of the triangle
# in order to evaluate the basis functions at the 
# center of gravity (cog)

# it shows, in particular, that in general 
# the nodal functions are not zero at the cog
# and that the P3-exact Newton-Cotes quadrature in 2D
# does not give a diagonal mass matrix

# the use of the lower order Newton-Cotes formulae
## cf. http://www.mathematik.uni-dortmund.de/~kuzmin/cfdintro/lecture5.pdf

# the P2-exact formula gives a singular mass matrix
# as does the P1-exact formula 

# in addition the use of the P1-exact formula introduces an
# consistency error of order 2

import sympy as sp

x, y = sp.symbols('x,y')
x1, x2, x3 = sp.symbols('x1, x2, x3')
y1, y2, y3 = sp.symbols('y1, y2, y3')

# ===> set the edge points here <===
subsdic = {x1:0,
        y1:0,
        x2:2,
        y2:0,
        x3:1,
        y3:1}

A = sp.Matrix([
        [1, x1, y1],
        [1, x2, y2],
        [1, x3, y3]]
        )
A1 = sp.Matrix([
        [1, x, y],
        [1, x2, y2],
        [1, x3, y3]]
        )
A2 = sp.Matrix([
        [1, x1, y1],
        [1, x, y],
        [1, x3, y3]]
        )
A3 = sp.Matrix([
        [1, x1, y1],
        [1, x2, y2],
        [1, x, y]]
        )

dA, dA1, dA2, dA3 = A.det(), A1.det(), A2.det(), A3.det() 

l1, l2, l3 = sp.Abs(dA1/dA), sp.Abs(dA2/dA), sp.Abs(dA3/dA)

phi1, phi2, phi3 = l1*(2*l1-1),  l2*(2*l2-1),  l3*(2*l3-1) 
phi12, phi13, phi23 = 4*l1*l2, 4*l1*l3, 4*l2*l3

#center of gravity
cogx = sp.Rational(1,3)*(x1+x2+x3)
cogy = sp.Rational(1,3)*(y1+y2+y3)

sbcog = {x:cogx, y:cogy}

sbnodesl = [{x:x1,y:y1},
        {x:x2,y:y2},
        {x:x3,y:y3},
        {x:0.5*(x1+x2),y:0.5*(y1+y2)},
        {x:0.5*(x1+x3),y:0.5*(y1+y3)},
        {x:0.5*(x2+x3),y:0.5*(y2+y3)}]

sbcogeva = {x:cogx.subs(subsdic), y:cogy.subs(subsdic)}


print 'The nodes of the triangle:', subsdic
print 'The center of gravity:', sbcogeva

phiL = [ phi1, phi2, phi3, phi12, phi13, phi23] 

for phi in phiL:
    print 'Value at cog: ', phi.subs({x:cogx, y:cogy}).subs(subsdic)
    print 'Values at the nodes:'
    for sbnod in sbnodesl:
        print phi.subs(sbnod).subs(subsdic)




