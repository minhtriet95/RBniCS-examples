from dolfin import *
from rbnics import *
import matplotlib.pyplot as plt


# @EIM()
# @PullBackFormsToReferenceDomain()
# @ShapeParametrization(
#     ("mu[0]*x[0]", "mu[0]*x[1]"),  # subdomain 1
# )
class ParameterizedCurve(EllipticCoerciveProblem):

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
        self.dx = Measure("dx")(subdomain_data=self.subdomains)
        self.ds = Measure("ds")(subdomain_data=self.boundaries)
        # ...
        self.n = FacetNormal(mesh)
        self.E = 1.0
        self.nu = 0.3
        self.lambda_1 = self.E * self.nu / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu))
        self.lambda_2 = self.E / (2.0 * (1.0 + self.nu))

    # Return custom problem name
    def name(self):
        return "ParameterizedCurve"

    # Return the alpha lower bound.
    def get_stability_factor_lower_bound(self):
        return 1.0

    # Return theta multiplicative terms of the affine expansion of the problem.
    def compute_theta(self, term):
        mu = self.mu
        if term == "a":
            theta_a0 = 1.0
            return (theta_a0, )
        elif term == "f":
            theta_f0 = mu[1]
            return (theta_f0, )
        else:
            raise ValueError("Invalid term for compute_theta().")

    # Return forms resulting from the discretization of the affine expansion of the problem operators.
    def assemble_operator(self, term):
        v = self.v
        dx = self.dx
        if term == "a":
            u = self.u
            a0 = self.elasticity(u, v) * dx
            return (a0, )
        elif term == "f":
            ds = self.ds
            n = self.n
            f0 = inner(Constant(-1.0)*n, v) * ds(4)
            return (f0, )
        elif term == "dirichlet_bc":
            bc0 = [DirichletBC(self.V.sub(1), Constant(0.0), self.boundaries, 1),
                   DirichletBC(self.V.sub(0), Constant(0.0), self.boundaries, 3)]
            return (bc0,)
        elif term == "inner_product":
            u = self.u
            x0 = inner(u, v) * dx + inner(grad(u), grad(v)) * dx
            return (x0,)
        else:
            raise ValueError("Invalid term for assemble_operator().")

    # Auxiliary function to compute the elasticity bilinear form
    def elasticity(self, u, v):
        lambda_1 = self.lambda_1
        lambda_2 = self.lambda_2
        return 2.0 * lambda_2 * inner(sym(grad(u)), sym(grad(v))) + lambda_1 * tr(sym(grad(u))) * tr(sym(grad(v)))

mesh = Mesh("data/arch_test.xml")
subdomains = MeshFunction("size_t", mesh, "data/arch_test_physical_region.xml")
boundaries = MeshFunction("size_t", mesh, "data/arch_test_facet_region.xml")

V = VectorFunctionSpace(mesh, "Lagrange", 1)

problem = ParameterizedCurve(V, subdomains=subdomains, boundaries=boundaries)
mu_range = [(0.5, 2.0), (-1.0, 1.0)]
problem.set_mu_range(mu_range)

reduction_method = ReducedBasis(problem)
reduction_method.set_Nmax(20, EIM=21)
reduction_method.set_tolerance(1e-5, EIM=1e-4)

reduction_method.initialize_training_set(100, EIM=100)
reduced_problem = reduction_method.offline()

online_mu = (0.5, 1.0)
reduced_problem.set_mu(online_mu)
reduced_solution = reduced_problem.solve()  # Obtain coefficient "weights"
reduced_problem.export_solution(filename="online_solution")
P = plot(reduced_solution, reduced_problem=reduced_problem)

P.set_cmap('jet')
plt.ylabel(r'$x_2$', fontsize=14)
plt.xlabel(r'$x_1$', fontsize=14)
plt.colorbar(P)
plt.show()
