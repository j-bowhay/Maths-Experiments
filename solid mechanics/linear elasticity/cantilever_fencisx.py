from dolfinx import mesh, fem, plot, io, default_scalar_type
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
import ufl
import numpy as np

# Parameters
L = 1
W = 0.2

rho = 2700  # density
E = 70e9  # Young's modulus
nu = 0.3  # Poisson ratio

# Lame parameters
mu = E/(2*(1 + nu))
lambda_ = (E*nu)/((1 + nu)*(1 - 2*nu))

g = 9.81  # Gravity

# Mesh

domain = mesh.create_box(MPI.COMM_WORLD, [np.array([0, 0, 0]), np.array([L, W, W])],
                         [20, 6, 6], cell_type=mesh.CellType.hexahedron)
V = fem.functionspace(domain, ("Lagrange", 1, (domain.geometry.dim, )))

# Boundary conditions
# Clamped at x = 0

def clamped_boundary(x):
    return np.isclose(x[0], 0)

fdim = domain.topology.dim - 1
boundary_facets = mesh.locate_entities_boundary(domain, fdim, clamped_boundary)

u_D = np.array([0, 0, 0], dtype=default_scalar_type)
bc = fem.dirichletbc(u_D, fem.locate_dofs_topological(V, fdim, boundary_facets), V)

# Traction free

T = fem.Constant(domain, default_scalar_type((0, 0, 0)))

ds = ufl.Measure("ds", domain=domain)

# Variational problem

def epsilon(u):
    return ufl.sym(ufl.grad(u))


def sigma(u):
    e = epsilon(u)
    return lambda_*ufl.tr(e)*ufl.Identity(len(u)) + 2*mu*e

# Function spaces
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

# Body force
f = fem.Constant(domain, default_scalar_type((0, 0, -rho * g)))

# Bilinear form
a = ufl.inner(sigma(u), epsilon(v))*ufl.dx

# Linear form
L = ufl.dot(f, v)*ufl.dx + ufl.dot(T, v)*ds

problem = LinearProblem(a, L, bcs=[bc],
                        petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
uh = problem.solve()

with io.XDMFFile(domain.comm, "deformation.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    uh.name = "Deformation"
    xdmf.write_function(uh)