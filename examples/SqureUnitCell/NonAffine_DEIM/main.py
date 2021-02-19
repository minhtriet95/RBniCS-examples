from dolfin import *
from rbnics import *
from subfiles.clear_offline_data import Cleaner
import matplotlib.pyplot as plt

from matplotlib import rc
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)

@DEIM(basis_generation="POD")
@PullBackFormsToReferenceDomain()
@ShapeParametrization(
    ("(mu[0]/0.5)*x[0]", "(mu[0]/0.5)*x[1]"),  # subdomain 1
    ("(-4.0*(2.0*mu[0] - 1.0)*(sqrt(x[0]**2 + x[1]**2) - 0.75) + 1.0)*x[0]",
     "(-4.0*(2.0*mu[0] - 1.0)*(sqrt(x[0]**2 + x[1]**2) - 0.75) + 1.0)*x[1]"),  # subdomain 2
    ("x[0]", "x[1]"),  # subdomain 3
)
class UnitCell(EllipticCoerciveProblem):

    # Default initialization of members
    # @generate_function_space_for_stability_factor
    def __init__(self, V, **kwargs):
        # Call the standard initialization
        EllipticCoerciveProblem.__init__(self, V, **kwargs)
        # ... and also store FEniCS data structures for assembly
        assert "subdomains" in kwargs
        assert "boundaries" in kwargs
        self.subdomains, self.boundaries = kwargs["subdomains"], kwargs["boundaries"]
        self.u = TrialFunction(V)
        self.v = TestFunction(V)
        self.dx = Measure("dx")(subdomain_data=self.subdomains)
        self.ds = Measure("ds")(subdomain_data=self.boundaries)

    # Return custom problem name
    def name(self):
        return "UnitCell"

    # Return the alpha lower bound.
    def get_stability_factor_lower_bound(self):
        return 1.0

    # Return theta multiplicative terms of the affine expansion of the problem.
    # @compute_theta_for_stability_factor
    def compute_theta(self, term):
        mu = self.mu
        if term == "a":
            theta_a0 = 10.0
            theta_a1 = 1.0
            theta_a2 = 1.0
            return (theta_a0, theta_a1, theta_a2)
        elif term == "f":
            theta_f0 = 1.0
            return (theta_f0, )
        else:
            raise ValueError("Invalid term for compute_theta().")

    # Return forms resulting from the discretization of the affine expansion of the problem operators.
    # @assemble_operator_for_stability_factor
    def assemble_operator(self, term):
        v = self.v
        dx = self.dx
        if term == "a":
            u = self.u
            a0 = inner(grad(u), grad(v)) * dx(1)
            a1 = inner(grad(u), grad(v)) * dx(2)
            a2 = inner(grad(u), grad(v)) * dx(3)
            return (a0, a1, a2)
        elif term == "f":
            ds = self.ds
            f0 = v * ds(1)
            return (f0, )
        elif term == "dirichlet_bc":
            bc0 = [DirichletBC(self.V, Constant(0.0), self.boundaries, 3)]
            return (bc0,)
        elif term == "inner_product":
            u = self.u
            x0 = inner(grad(u), grad(v)) * dx
            return (x0,)
        else:
            raise ValueError("Invalid term for assemble_operator().")


# 1. Read the mesh for this problem
mesh = Mesh("data/thermal_block.xml")
subdomains = MeshFunction(
    "size_t", mesh, "data/thermal_block_physical_region.xml")
boundaries = MeshFunction(
    "size_t", mesh, "data/thermal_block_facet_region.xml")

# 2. Create Finite Element space for this problem
V = FunctionSpace(mesh, "Lagrange", 1)

# 3. Allocate an object for UnitCell class
problem = UnitCell(V, subdomains=subdomains, boundaries=boundaries)
mu_range = [(0.3, 0.6)]
problem.set_mu_range(mu_range)

# 4. Prepare reduction with a ReducedBasis method
reduction_method = ReducedBasis(problem)
reduction_method.set_Nmax(50, DEIM=20)
reduction_method.set_tolerance(1e-4, DEIM=1e-3)

# 5. Perform the offline phase
solve_stage = "online"
reduction_method.initialize_training_set(50, DEIM=25)
# Delete old offline data
if solve_stage == "offline":
    clean = Cleaner(problem.name())
    clean.delete_folder()
elif solve_stage == "online":
    pass
else:
    print("Invalid name for solving stage!")
reduced_problem = reduction_method.offline()

# 6. Perform an online solve
online_mu = (0.6, )
reduced_problem.set_mu(online_mu)
reduced_solution = reduced_problem.solve()
reduced_problem.export_solution(filename="online_solution")

# 7. Visualize results
P = plot(reduced_solution, reduced_problem=reduced_problem, levels=256)

P.set_cmap('jet')
plt.ylabel(r'$x_2$', fontsize=16)
plt.xlabel(r'$x_1$', fontsize=16)
plt.colorbar(P)
plt.savefig('rb_solution.pdf')
plt.show()


# 8. Perform an error analysis
# reduction_method.initialize_testing_set(100, EIM=100)
# reduction_method.error_analysis()

# reduction_method.initialize_testing_set(100, EIM=100)
# reduction_method.speedup_analysis()