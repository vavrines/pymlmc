#%%
from firedrake import *
from numpy import ceil, exp

class DarcyScalarProblem(object):
    # l: current level
    # M: refinement factor
    def __init__(self, l, M, cost = None):
        self.init_problem(l, M)
        self.cost = cost

    def init_problem(self, l, M):
        N = int(ceil(4*M**l))
        mesh = UnitSquareMesh(N, N)

        V = FunctionSpace(mesh, "CG", 1)
        self.V = V
        u = Function(V)
        v = TestFunction(V)
        x,y = SpatialCoordinate(mesh)
        f = Function(V).interpolate(sin(2.0*pi*x))

        K0 = Constant(1.) # Deterministic field
        self.Kw = Constant(0.0) # Random scalar
        K = K0 + self.Kw # Diffusion coefficient

        bcs = [DirichletBC(V, Constant(1.0), [1]), DirichletBC(V, Constant(0.0), [2])]

        a = inner(K*grad(u), grad(v))*dx
        L = f*v*dx
        F = a-L

        self.u = u
        problem = NonlinearVariationalProblem(F, self.u, bcs)
        self.solver = NonlinearVariationalSolver(problem, solver_parameters = {"ksp_type": "preonly", "pc_type": "lu"})

    def evaluate(self, sample):
        self.Kw.assign(exp(sample))
        self.u.assign(Function(self.V))
        self.solver.solve()
        return norm(self.u)**2

#%%
from pymlmc import mlmc_test, mlmc_plot, mlmc_fn
from numpy import sqrt
from numpy.random import randn
from matplotlib import pyplot as plt

# %%
N0 = 10 # initial samples on coarse levels
Lmin = 2  # minimum refinement level
Lmax = 7 # maximum refinement level
M = 2 # refinement factor
N = 200 # samples for convergence tests
L = 7 # levels for convergence tests
Eps = [0.000002, 0.000005, 0.00001, 0.00002, 0.00005, 0.0001]#, 0.0002, 0.0005]
# %%
sig = 1.0 # standard deviation 
def sampler(N, l):
    sample = sqrt(sig)*randn(N)
    return (sample, sample)
l_range = range(Lmax+1)
# %%
problems = [DarcyScalarProblem(l, M) for l in l_range]
def darcyscalar_l(l, N):
    return mlmc_fn(l, N, problems, sampler = sampler)
# %%
