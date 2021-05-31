from dolfin import *
from rbnics import *
from matplotlib import pyplot as plt
from tools.clear_offline_data import Cleaner

from matplotlib import rc
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)

@EIM()
class Gaussian(EllipticCoerciveProblem):

    # Default initialization of members
    def __init__(self, V, **kwargs):
        # Call the standard initialization
        EllipticCoerciveProblem.__init__(self, V, **kwargs)
        # ... and also store FEniCS data structures for assembly
        assert "subdomains" in kwargs
        assert "boundaries" in kwargs
        self.subdomains, self.boundaries = kwargs["subdomains"], kwargs["boundaries"]
        self.u = TrialFunction(V)
        self.v = TestFunction(V)
        self.dx = Measure("dx")(subdomain_data=subdomains)
        self.g = ParametrizedExpression(
            self, "exp(- 2 * pow(x[0] - mu[0], 2) - 2 * pow(x[1] - mu[1], 2))", mu=(0., 0.),
            element=V.ufl_element())
        # note that we cannot use self.mu in the initialization of self.f, because self.mu has not been initialized yet

    # Return custom problem name
    def name(self):
        return "GaussianEIM"

    # Return the alpha_lower bound.
    def get_stability_factor_lower_bound(self):
        return 1.

    # Return theta multiplicative terms of the affine expansion of the problem.
    def compute_theta(self, term):
        if term == "a":
            return (1.,)
        elif term == "f":
            return (1.,)
        else:
            raise ValueError("Invalid term for compute_theta().")

    # Return forms resulting from the discretization of the affine expansion of the problem operators.
    def assemble_operator(self, term):
        v = self.v
        dx = self.dx
        if term == "a":
            u = self.u
            a0 = inner(grad(u), grad(v)) * dx
            return (a0,)
        elif term == "f":
            g = self.g
            f0 = g * v * dx
            return (f0,)
        elif term == "dirichlet_bc":
            bc0 = [DirichletBC(self.V, Constant(0.0), self.boundaries, 1),
                   DirichletBC(self.V, Constant(0.0), self.boundaries, 2),
                   DirichletBC(self.V, Constant(0.0), self.boundaries, 3)]
            return (bc0,)
        elif term == "inner_product":
            u = self.u
            x0 = inner(grad(u), grad(v)) * dx
            return (x0,)
        else:
            raise ValueError("Invalid term for assemble_operator().")


# 1. Read the mesh for this problem
mesh = Mesh("data/gaussian.xml")
subdomains = MeshFunction("size_t", mesh, "data/gaussian_physical_region.xml")
boundaries = MeshFunction("size_t", mesh, "data/gaussian_facet_region.xml")

# 2. Create Finite Element space for this problem
V = FunctionSpace(mesh, "Lagrange", 1)

# 3. Allocate an object for UnitCell class
problem = Gaussian(V, subdomains=subdomains, boundaries=boundaries)
mu_range = [(-1.0, 1.0), (-1.0, 1.0)]
problem.set_mu_range(mu_range)

# 4. Prepare reduction with a ReducedBasis method
reduction_method = ReducedBasis(problem)
reduction_method.set_Nmax(50, EIM=20)
reduction_method.set_tolerance(1e-4, EIM=1e-3)

# 5. Perform the offline phase
solve_stage = "online"
reduction_method.initialize_training_set(50, EIM=50)
# Delete old offline data
if solve_stage == "offline":
    clean = Cleaner(problem.name())
    names = clean.list_names(clean.path)
    if problem.name() in names:
        clean.delete_offline_folder()
elif solve_stage == "online":
    pass
else:
    print("Invalid name for solving stage!")
reduced_problem = reduction_method.offline()

# 6. Perform an online solve
online_mu = (0.7, 0.7)
reduced_problem.set_mu(online_mu)
reduced_solution = reduced_problem.solve()

# 7. Visualize results
P = plot(reduced_solution, reduced_problem=reduced_problem)
P.set_cmap('jet')
plt.ylabel(r'$x_2$', fontsize=16)
plt.xlabel(r'$x_1$', fontsize=16)
plt.colorbar(P)
plt.savefig('rb_solution.pdf')
plt.show()
