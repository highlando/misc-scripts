import dolfin
import scipy.sparse as sps
import scipy.io

from dolfin import dx, grad, inner

# set the backend
dolfin.parameters.linear_algebra_backend = "uBLAS"


def mat_dolfin2sparse(A):
    rows, cols, values = A.data()
    return sps.csr_matrix((values, cols, rows))

# example usage
mesh = dolfin.UnitSquareMesh(10, 10)
V = dolfin.VectorFunctionSpace(mesh, 'CG', 1)
v = dolfin.TestFunction(V)
u = dolfin.TrialFunction(V)

# define the form and assemble
lap_form = inner(grad(u), grad(v))*dx
lap_ass = dolfin.assemble(lap_form)

# convert to csr matrix
lap_mat = mat_dolfin2sparse(lap_ass)

scipy.io.savemat('matfile', dict(A=lap_mat))
# gives the file matfile.mat with lap_mat saved as A
