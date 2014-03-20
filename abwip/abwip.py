from dolfin import *
import krypy
import numpy as np
import scipy.sparse as sps
import matplotlib.pyplot as plt

import dolfin_to_nparrays as dtn
import time_int_schemes as tis


parameters.linear_algebra_backend = 'uBLAS'

def time_int_params():
    t0 = 0.0
    tE = 1.0
    tip = dict(t0 = t0,
            tE = tE,
            Residuals = NseResiduals(),
            nu = 1
            )

    return tip

def abetterworld(N = 20, Nts = [16, 32, 64, 128]):

    tip = time_int_params()
    prp = problem_params(N)

    fv, fp, sol_v, sol_p = prp['fv'], prp['fp'], prp['sol_v'], prp['sol_p']

    # fixing some parameters
    # fv.om, v_sol.om, p_sol.om = 1, 1, 1
    # fv.nu, v_sol.nu, p_sol.nu = 1, 1, 1

###
## start with the Stokes problem for initialization
###

    stokesmats = dtn.get_stokessysmats(prp['V'], prp['Q'],
                                         tip['nu'])
    rhsd_vf = dtn.setget_rhs(prp['V'], prp['Q'], 
                            prp['fv'], prp['fp'], t=0)

    # remove the freedom in the pressure 
    stokesmats['J'] = stokesmats['J'][:-1,:][:,:]
    stokesmats['JT'] = stokesmats['JT'][:,:-1][:,:]
    rhsd_vf['fp'] = rhsd_vf['fp'][:-1,:]

    # reduce the matrices by resolving the BCs
    (stokesmatsc, 
            rhsd_stbc, 
            INVINDS, 
            bcinds, 
            bcvals) = dtn.condense_sysmatsbybcs(stokesmats,
                                                prp['bc0'])

    # casting some parameters 
    NV, NP = len(INVINDS), stokesmats['J'].shape[0]
    M, A, J = stokesmatsc['M'], stokesmatsc['A'], stokesmatsc['J'] 


###
## Compute the time-dependent flow
###

    inivalvec = np.zeros((NV+NP,1))

    for nts in Nts:
        DT = (1.0*tip['t0'] - tip['tE'])/nts

        biga = sps.vstack([
                    sps.hstack([M+DT*A, J.T]),
                    sps.hstack([J, sps.csr_matrix((NP, NP))])
                        ])

        vp_old = inivalvec

        ContiRes, VelEr, PEr = [], [], []

        for tcur in np.linspace(tip['t0']+DT,tip['tE'],nts):

            fvpn = dtn.setget_rhs(prp['V'], prp['Q'], fv, fp, t=tcur)

            cur_rhs = np.vstack([fvpn['fv'][INVINDS,:],
                                np.zeros((NP,1))])

            vp_old = krypy.linsys.minres(biga, cur_rhs,
                    x0=vp_old, maxiter=100)['xk'] 

            vc = vp_old[:NV,]
            pc = vp_old[NV:,]

            v, p = tis.expand_vp_dolfunc(tip, vp=None, vc=vc, pc=pc)

            v_sol.t, p_sol.t = tcur

            # the errors  
            ContiRes.append(tis.comp_cont_error(v,fp,tip['Q']))
            VelEr.append(errornorm(vCur,v))
            PEr.append(errornorm(pCur,p))

        tip['Residuals'].ContiRes.append(ContiRes)
        tip['Residuals'].VelEr.append(VelEr)
        tip['Residuals'].PEr.append(PEr)

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
    rhs3 = smp.simplify(diff(u1,x) + diff(u2,y))

    sol_p = Expression(ccode(p), om=1, t=0)
    sol_v = Expression((ccode(u1), ccode(u2)), om=1, t=0)
    fv = Expression((ccode(rhs1), ccode(rhs2)), om=1, nu=1, t=0) 
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
    bc0 = [DirichletBC(V, noslip, 'on_boundary')]

    sol_v, sol_p, fv, fp = exact_stokes_sol()

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
    abetterworld()
