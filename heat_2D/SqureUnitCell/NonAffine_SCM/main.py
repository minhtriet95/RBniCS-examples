from dolfin import *
from rbnics import *
import matplotlib.pyplot as plt

@DEIM()
# @SCM()
@PullBackFormsToReferenceDomain()
@ShapeParametrization(
    ("x[0]", "x[1]"),  # subdomain 1
    ("(0.5*(0.5-sqrt(x[0]**2 + x[1]**2)) + mu[2]*(sqrt(x[0]**2 + x[1]**2)-0.25))/(0.5*(0.5-0.25))*x[0]",
     "(0.5*(0.5-sqrt(x[0]**2 + x[1]**2)) + mu[2]*(sqrt(x[0]**2 + x[1]**2)-0.25))/(0.5*(0.5-0.25))*x[1]"),  # subdomain 2
    ("(0.5*(0.5-sqrt(x[0]**2 + x[1]**2)) + mu[2]*(sqrt(x[0]**2 + x[1]**2)-0.75))/(0.5*(0.5-0.75))*x[0]",
     "(0.5*(0.5-sqrt(x[0]**2 + x[1]**2)) + mu[2]*(sqrt(x[0]**2 + x[1]**2)-0.75))/(0.5*(0.5-0.75))*x[1]"),  # subdomain 3
    ("x[0]", "x[1]"),  # subdomain 4
)
class UnitSquare(EllipticCoerciveProblem):

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
        # Customize eigen solver parameters (remove when not using SCM)
        # self._eigen_solver_parameters.update({
        #     "bounding_box_minimum": {
        #         "problem_type": "gen_hermitian", "spectral_transform": "shift-and-invert",
        #         "spectral_shift": 1.e-5, "linear_solver": "mumps"
        #     },
        #     "bounding_box_maximum": {
        #         "problem_type": "gen_hermitian", "spectral_transform": "shift-and-invert",
        #         "spectral_shift": 1.e5, "linear_solver": "mumps"
        #     },
        #     "stability_factor": {
        #         "problem_type": "gen_hermitian", "spectral_transform": "shift-and-invert",
        #         "spectral_shift": 1.e-5, "linear_solver": "mumps"
        #     }
        # })

    # Return custom problem name
    def name(self):
        return "UnitSquare"

    # Return the alpha_lower bound.
    def get_stability_factor_lower_bound(self):
        return 1.

    # Return theta multiplicative terms of the affine expansion of the problem.
    # @compute_theta_for_stability_factor
    def compute_theta(self, term):
        mu = self.mu
        if term == "a":
            theta_a0 = mu[0]
            theta_a1 = mu[0]
            theta_a2 = 1.
            theta_a3 = 1.
            return (theta_a0, theta_a1, theta_a2, theta_a3)
        elif term == "f":
            theta_f0 = mu[1]
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
            a3 = inner(grad(u), grad(v)) * dx(4)
            return (a0, a1, a2, a3)
        elif term == "f":
            ds = self.ds
            f0 = v * ds(1)
            return (f0,)
        elif term == "dirichlet_bc":
            bc0 = [DirichletBC(self.V, Constant(0.0), self.boundaries, 3)]
            return (bc0,)
        elif term == "inner_product":
            u = self.u
            x0 = inner(grad(u), grad(v)) * dx
            return (x0,)
        else:
            raise ValueError("Invalid term for assemble_operator().")


mesh = Mesh("data/unit_square.xml")
subdomains = MeshFunction("size_t", mesh, "data/unit_square_physical_region.xml")
boundaries = MeshFunction("size_t", mesh, "data/unit_square_facet_region.xml")

V = FunctionSpace(mesh, "Lagrange", 1)

problem = UnitSquare(V, subdomains=subdomains, boundaries=boundaries)
mu_range = [(0.1, 10.0), (-1.0, 1.0), (0.4, 0.6)]
problem.set_mu_range(mu_range)

reduction_method = ReducedBasis(problem)
reduction_method.set_Nmax(30, DEIM=21, SCM=25)
reduction_method.set_tolerance(1e-4, DEIM=1e-3, SCM=0.75)

reduction_method.initialize_training_set(50, DEIM=100, SCM=100)
reduced_problem = reduction_method.offline()

online_mu = (1.0, -1.0, 0.6)
reduced_problem.set_mu(online_mu)
reduced_solution = reduced_problem.solve()
P = plot(reduced_solution, reduced_problem=reduced_problem)

P.set_cmap('jet')
plt.ylabel(r'$x_2$', fontsize=14)
plt.xlabel(r'$x_1$', fontsize=14)
plt.colorbar(P)
# plt.savefig('rb_solution.pdf')
plt.show()

# reduction_method.initialize_testing_set(100)
# reduction_method.error_analysis()

# reduction_method.initialize_testing_set(100)
# reduction_method.speedup_analysis()
