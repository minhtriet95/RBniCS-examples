from dolfin import *
from rbnics import *
import matplotlib.pyplot as plt

from matplotlib import rc
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)

@SCM()
class ThermalBlock(EllipticCoerciveProblem):

    # Default initialization of members
    @generate_function_space_for_stability_factor
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
        # Customize eigen solver parameters (remove when not using SCM)
        self._eigen_solver_parameters.update({
            "bounding_box_minimum": {
                "problem_type": "gen_hermitian", "spectral_transform": "shift-and-invert",
                "spectral_shift": 1.e-5, "linear_solver": "mumps"
            },
            "bounding_box_maximum": {
                "problem_type": "gen_hermitian", "spectral_transform": "shift-and-invert",
                "spectral_shift": 1.e5, "linear_solver": "mumps"
            },
            "stability_factor": {
                "problem_type": "gen_hermitian", "spectral_transform": "shift-and-invert",
                "spectral_shift": 1.e-5, "linear_solver": "mumps"
            }
        })

    # Return custom problem name
    def name(self):
        return "ThermalBlock"

    # Return theta multiplicative terms of the affine expansion of the problem.
    @compute_theta_for_stability_factor
    def compute_theta(self, term):
        mu = self.mu
        if term == "a":
            theta_a0 = mu[0]
            theta_a1 = mu[1]
            theta_a2 = mu[2]
            theta_a3 = mu[3]
            theta_a4 = mu[4]
            theta_a5 = mu[5]
            theta_a6 = mu[6]
            theta_a7 = mu[7]
            theta_a8 = 1.
            return (theta_a0, theta_a1, theta_a2, theta_a3, theta_a4, theta_a5, theta_a6, theta_a7, theta_a8)
        elif term == "f":
            theta_f0 = 1.0
            return (theta_f0, )
        else:
            raise ValueError("Invalid term for compute_theta().")

    # Return forms resulting from the discretization of the affine expansion of the problem operators.
    @assemble_operator_for_stability_factor
    def assemble_operator(self, term):
        v = self.v
        dx = self.dx
        if term == "a":
            u = self.u
            a0 = inner(grad(u), grad(v)) * dx(1)
            a1 = inner(grad(u), grad(v)) * dx(2)
            a2 = inner(grad(u), grad(v)) * dx(3)
            a3 = inner(grad(u), grad(v)) * dx(4)
            a4 = inner(grad(u), grad(v)) * dx(5)
            a5 = inner(grad(u), grad(v)) * dx(6)
            a6 = inner(grad(u), grad(v)) * dx(7)
            a7 = inner(grad(u), grad(v)) * dx(8)
            a8 = inner(grad(u), grad(v)) * dx(9)
            return (a0, a1, a2, a3, a4, a5, a6, a7, a8)
        elif term == "f":
            ds = self.ds
            f0 = v * ds(1)
            return (f0, )
        elif term == "dirichlet_bc":
            bc0 = [DirichletBC(self.V, Constant(0.0), self.boundaries, 5)]
            return (bc0,)
        elif term == "inner_product":
            u = self.u
            x0 = inner(grad(u), grad(v)) * dx
            return (x0,)
        else:
            raise ValueError("Invalid term for assemble_operator().")


mesh = Mesh("data/elastic_block.xml")
subdomains = MeshFunction("size_t", mesh, "data/elastic_block_physical_region.xml")
boundaries = MeshFunction("size_t", mesh, "data/elastic_block_facet_region.xml")

V = FunctionSpace(mesh, "Lagrange", 1)

problem = ThermalBlock(V, subdomains=subdomains, boundaries=boundaries)
mu_range = [
    (0.1, 10.0),
    (0.1, 10.0),
    (0.1, 10.0),
    (0.1, 10.0),
    (0.1, 10.0),
    (0.1, 10.0),
    (0.1, 10.0),
    (0.1, 10.0),
]
problem.set_mu_range(mu_range)

reduction_method = ReducedBasis(problem)
reduction_method.set_Nmax(50, SCM=20)
reduction_method.set_tolerance(1e-5, SCM=7.5e-1)

reduction_method.initialize_training_set(100, SCM=500)
reduced_problem = reduction_method.offline()

online_mu = (1.0, 1.0, 1.0, 1.0, 5.0, 1.0, 1.0, 1.0)
reduced_problem.set_mu(online_mu)
reduced_solution = reduced_problem.solve()
P = plot(reduced_solution, reduced_problem=reduced_problem)

# Define default LaTeX fonts
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

P.set_cmap('jet')
plt.ylabel(r'$x_2$', fontsize=14)
plt.xlabel(r'$x_1$', fontsize=14)
plt.colorbar(P)
plt.savefig('rb_solution.pdf')
plt.show()

# reduction_method.initialize_testing_set(100)
# reduction_method.error_analysis()

# reduction_method.initialize_testing_set(100)
# reduction_method.speedup_analysis()
