from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
import timeit

start_time = timeit.default_timer()

# Set unit
GPa = 10**9
MPa = 10**6

# Import mesh, subdomains, and boundaries from files
mesh = Mesh("data/vessel.xml")
subdomains = MeshFunction("size_t", mesh, "data/vessel_physical_region.xml")
boundaries = MeshFunction("size_t", mesh, "data/vessel_facet_region.xml")

# Young's modulus and Poisson ratio
E = Constant(201.0*GPa)
nu = Constant(0.306)

# Lame constants for constitutive relation
mu = E/2/(1+nu)
lmbda = E*nu/(1+nu)/(1-2*nu)

# Define funtion space and basis functions
V = VectorFunctionSpace(mesh, 'Lagrange', degree=2)
u = TrialFunction(V)
v = TestFunction(V)

# Define Dirichlet boundary conditions
bcs = [DirichletBC(V, Constant((0., 0., 0.)), boundaries, 1),  # fully fixed (u, v, w = 0)
       DirichletBC(V.sub(0), Constant(0.0), boundaries, 2),    # displacement u = 0
       DirichletBC(V.sub(1), Constant(0.0), boundaries, 3)]    # displacement v = 0

# Define strain and stress
def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)

def sigma(u):
    dim = u.geometric_dimension()  # domain dimension
    return lmbda*tr(epsilon(u))*Identity(dim) + 2.0*mu*epsilon(u)

# Define variational problem
n = FacetNormal(mesh)
f = -40.0*MPa

dx = Measure("dx")(subdomain_data=subdomains)
ds = Measure("ds")(subdomain_data=boundaries)

a = inner(sigma(u), epsilon(v))*dx
L = inner(f*n, v)*ds(4)

# Compute solution
u = Function(V, name="displacement")
solve(a == L, u, bcs)

# Compute stresses
s = sigma(u) - (1./3)*tr(sigma(u))*Identity(3)
von_Mises = sqrt(3./2*inner(s, s))
V = FunctionSpace(mesh, 'Lagrange', 2)
von_Mises = project(von_Mises, V)
von_Mises.rename("von_Mises", "Von Mises")

file_results = XDMFFile("vessel_results.xdmf")
file_results.parameters["flush_output"] = True
file_results.parameters["functions_share_mesh"] = True
file_results.write(u, 0.)
file_results.write(von_Mises, 0.)

stop_time = timeit.default_timer()
run_time = round(stop_time - start_time, 2)

print(f'Finish solving in {run_time}s.')
