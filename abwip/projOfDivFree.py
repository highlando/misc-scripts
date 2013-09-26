# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# On the $L^2$-FEM-Projection of Solenoidal Functions
# ---
# 
# An example that illustrates that the $L^2$ projection of a divergence free velocity field onto a 
# 
# *  stable 
# *  mixed FEM
# 
# function space is **not necessarily divergence free**!
# 
# We use a uniform triangularization of the unit-square and the socalled *Taylor-Hood* elements.
# 
# The implementation is in python and uses the `dolfin` interface (Version 1.2.0) to the open source FEM package [FEniCs](http://fenicsproject.org/).

# <codecell>

import sympy as smp
import numpy as np
from dolfin import *

import dolfin_to_nparrays as dtp
# <markdowncell>

# We set up a divergence free velocity field `sol_u` on the unit square. Also, `sol_u` is zero at the boundary.

# <codecell>

x, y = smp.symbols('x[0], x[1]')
u_x = x*x*(1-x)*(1-x)*2*y*(1-y)*(2*y-1)
u_y = y*y*(1-y)*(1-y)*2*x*(1-x)*(1-2*x)
print 'div u = {}'.format( smp.simplify(smp.diff(u_x,x) + smp.diff(u_y,y)) )

# <markdowncell>

# Turn the symbolic expression into a Fenics expression.

# <codecell>

from sympy.printing import ccode
sol_u = Expression((ccode(u_x), ccode(u_y)))

# <markdowncell>

# Then, we set up the $P_2-P_1$ (Taylor-Hood) discretization of $H_0^1-L^2$ on the unit square and project `sol_u` onto it. 

# <codecell>

# function spaces
t0, tE = 0, 1
DT = 1e-1
N = 12
nu = 1

mesh = UnitSquareMesh(N, N)
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)

# No-slip boundary condition for velocity
noslip = Constant((0.0, 0.0))
bc0 = DirichletBC(V, noslip, 'on_boundary')

# Stokes Matrices
stms = get_stokessysmats( V, Q, nu)

M, A, J = stmsc['M'], stmsc['A'], stmsc['J'],
# remove the boundary nodes
stmsc = condense_sysmatsbybcs(stms, bc0)[0]

biga = sps.vstack([
                    sps.hstack([M+DT*A, J.T])
                    sps.hstack([J, sps.csr_matrix((Np,Np))])
                   ])

# projection
u_fun = project(sol_u, V, solver_type = "lu", bcs=bc0)

# <markdowncell>

# The divergence in the $P_2-P_1$ FEM space is given as the divergence of a $P_2$-function tested against all basis functions of $P_1$. 
# 
# In our case we compute `div(u_fun)` and assemble the vector of the divergence by testing against all `q` in `Q`.
# 
# Finally, we compute the norm (in `Q`) of the divergence of `u_fun` and divide it by the norm (in `V`) of `u_fun`, i.e. the relative residual in the discrete continuity equation. 

# <codecell>

q = TestFunction(Q)

disc_div_u = assemble( q*div(u_fun) *dx)

rel_div_res = norm(disc_div_u)/norm(u_fun)
print '\n And the relative residual in the \n discrete divergence for the projected u is:\n\n {}\n'.format(rel_div_res)

# <markdowncell>

# Questions or comments? Please visit me at [www.janheiland.de](http://www.janheiland.de).

