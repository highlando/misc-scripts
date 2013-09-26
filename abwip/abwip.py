from dolfin import *
import numpy as np
import scipy.sparse as sps
import matplotlib.pyplot as plt

import dolfin_to_nparrays as dtn
import linsolv_utils

import data_output_utils as dou

parameters.linear_algebra_backend = 'uBLAS'

def time_int_params(Nts):
    t0 = 0.0
    tE = 1.0
    dt = (tE - t0)/Nts
    tip = dict(t0 = t0,
            tE = tE,
            dt = dt, 
            Nts = Nts,
            Residuals = NseResiduals(), 
            ParaviewOutput = True, 
            nu = 1e-3,
            )

    return tip

def abetterworld(N = 20, Nts = 4):

    tip = time_int_params(Nts)
    prp = problem_params(N)
    prp['fv'].om = 1
    prp['p_sol'].om = 1
    prp['v_sol'].om = 1

###
## start with the Stokes problem for initialization
###

    stokesmats = dtn.get_stokessysmats(femp['V'], femp['Q'],
                                         tip['nu'])
    rhsd_vf = dtn.setget_rhs(femp['V'], femp['Q'], 
                            femp['fv'], femp['fp'], t=0)

    # remove the freedom in the pressure 
    stokesmats['J'] = stokesmats['J'][:-1,:][:,:]
    stokesmats['JT'] = stokesmats['JT'][:,:-1][:,:]
    rhsd_vf['fp'] = rhsd_vf['fp'][:-1,:]

    # reduce the matrices by resolving the BCs
    (stokesmatsc, 
            rhsd_stbc, 
            invinds, 
            bcinds, 
            bcvals) = dtn.condense_sysmatsbybcs(stokesmats,
                                                femp['bc0'])

    # casting some parameters 
    NV, NP = len(femp['invinds']), stokesmats['J']/shape[0]
    DT, INVINDS = tip['dt'], femp['invinds']
    M, A, J = stokesmatsc['M'], stokesmatsc['A'], stokesmatsc['J'] 


###
## Compute the time-dependent flow
###

    inivalvec = np.zeros((NV+NP,1))

    for newtk in range(1, tip['nnewtsteps']+1):

        biga = sps.vstack([
                    sps.hstack([M+DT*A, J.T])
                    sps.hstack([J, sps.csr_matrix((Np,Np))])
                        ])

        v_old = inivalvec
        for t in np.arange(tip['t0'], tip['tE'], DT):

            fvpn = dtn.setget_rhs(femp['V'], femp['Q'], femp['fv'], femp['fp'], t=t)

            cur_rhs = np.vstack([fvpn['fv'][INVINDS,:],
                                np.zeros((NP,1))])

            ret = krypy.linsys.minres(IterA, Iterrhs, 
                    x0=vp_old, tol=TolCor*TsP.linatol,
                    M=MInv)
            vp_old = ret['xk'] 

            v_old = vp[:NV,]

        tip['norm_nwtnupd'].append(norm_nwtnupd)

    print tip['norm_nwtnupd']

def exact_stokes_sol():
    import sympy as smp
    from sympy import diff, sin, cos, pi
    from sympy.printing import ccode

    x, y, t, nu, om = smp.symbols('x[0],x[1],t,nu,om')

    ft = smp.sin(om*t)
    u1 = ft*x*x*(1-x)*(1-x)*2*y*(1-y)*(2*y-1)
    u2 = ft*y*y*(1-y)*(1-y)*2*x*(1-x)*(1-2*x)
    p = ft*x*(1-x)*y*(1-y)

    # Stokes case
    rhs1 = smp.simplify(diff(u1,t) - nu*(diff(u1,x,x) + diff(u1,y,y)) + diff(p,x))
    rhs2 = smp.simplify(diff(u2,t) - nu*(diff(u2,x,x) + diff(u2,y,y)) + diff(p,y))

    sol_p = Expression(ccode(p))
    sol_v = Expression((ccode(u1), ccode(u2)))
    fv = Expression((ccode(rhs1), ccode(rhs2))) 
    fp = Expression(ccode(rhs3))

    return sol_v, sol_p, fv, fp

def problem_params(N):
    """dictionary for the fem items of the (unit) driven cavity

    """
    mesh = UnitSquareMesh(N, N)
    V = VectorFunctionSpace(mesh, "CG", 2)
    Q = FunctionSpace(mesh, "CG", 1)

    # No-slip boundary condition for velocity
    noslip = Constant((0.0, 0.0))
    bc0 = DirichletBC(V, noslip, 'on_boundary')

    sol_v, sol_p, fv, fp = exact_stokes_sol(): 

    dfems = dict(mesh = mesh,
            V = V,
            Q = Q,
            bc0 = bc0,
            fv = fv,
            fp = fp,
            sol_v = sol_v,
            sol_p = sol_p)

    return dfems


class NseResiduals(object):
    def __init__(self):
        self.ContiRes = []
        self.VelEr = []
        self.PEr = []


if __name__ == '__main__':
    abwip()
