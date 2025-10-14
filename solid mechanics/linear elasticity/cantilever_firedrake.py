import firedrake as fd

# Parameters
L = 1
W = 0.2

rho = fd.Constant(2700)  # density
E = 70e9  # Young's modulus
nu = 0.3  # Poisson ratio

# Lame parameters
mu = fd.Constant(E/(2*(1 + nu)))
lambda_ = fd.Constant((E*nu)/((1 + nu)*(1 - 2*nu)))

g = fd.Constant(9.81)  # Gravity

# Mesh

domain = fd.BoxMesh(20, 6, 6, L, W, W)
V = fd.VectorFunctionSpace(domain, "Lagrange", 1)

# Boundary conditions
# Clamped at x = 0

bc = fd.DirichletBC(V, fd.Constant((0, 0, 0)), 1)

# Variational problem

def epsilon(u):
    return fd.sym(fd.grad(u))

def sigma(u):
    e = epsilon(u)
    return lambda_*fd.tr(e)*fd.Identity(domain.geometric_dimension()) + 2*mu*e

# Function spaces
u = fd.TrialFunction(V)
v = fd.TestFunction(V)

# Body force
f = fd.Constant((0, 0, -rho * g))

# Bilinear form
a = fd.inner(sigma(u), epsilon(v))*fd.dx

# Linear form
L = fd.dot(f, v)*fd.dx

uh = fd.Function(V, name="Displacement")
fd.solve(a == L, uh, bcs=bc,
         solver_parameters={"ksp_type": "preonly", "pc_type": "lu",
                            "ksp_monitor": None})

fd.VTKFile("output.pvd").write(uh)