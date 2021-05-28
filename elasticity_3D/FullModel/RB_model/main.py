from dolfin import *
from rbnics import *
import timeit
from subfiles.clear_offline_data import Cleaner

start_time = timeit.default_timer()

# Units
GPa = 10**9
MPa = 10**6

# Define parameters for mapping functions
# offset = "1.5131"
#* thickness = mu[0], radius = mu[1]

# # For cylindrical body (subdomain 1)
# x_norm_sub1 = "sqrt(x[0]**2 + x[1]**2)"
# mapping_x_sub1 = f"x[0] + (0.05-mu[0])*x[0]/0.2294*(-20*{x_norm_sub1}+1397./250) + (mu[1]-0.2794)*x[0]/{x_norm_sub1}"
# mapping_y_sub1 = f"x[1] + (0.05-mu[0])*x[1]/0.2294*(-20*{x_norm_sub1}+1397./250) + (mu[1]-0.2794)*x[1]/{x_norm_sub1}"
# mapping_z_sub1 = "x[2]"

# # For left spherical head (subdomain 2)
# z_sub2 = "(x[2] + 1.5131)"
# x_norm_sub2 = "sqrt(x[0]**2 + x[1]**2 + (x[2] + 1.5131)**2)"
# # x_norm_sub2 = "sqrt(x[0]**2 + x[1]**2 + x[2]**2)"
# mapping_x_sub2 = f"x[0] + (0.05-mu[0])*x[0]/0.2294*(-20*{x_norm_sub2}+1397./250) + (mu[1]-0.2794)*x[0]/{x_norm_sub2}"
# mapping_y_sub2 = f"x[1] + (0.05-mu[0])*x[1]/0.2294*(-20*{x_norm_sub2}+1397./250) + (mu[1]-0.2794)*x[1]/{x_norm_sub2}"
# mapping_z_sub2 = f"{z_sub2} + (0.05-mu[0])*{z_sub2}/0.2294*(-20*{x_norm_sub2}+1397./250) + (mu[1]-0.2794)*{z_sub2}/{x_norm_sub2} - 1.5131"

# # For right spherical head (subdomain 3)
# z_sub3 = "(x[2] - 1.5131)"
# x_norm_sub3 = "sqrt(x[0]**2 + x[1]**2 + (x[2] - 1.5131)**2)"
# mapping_x_sub3 = f"x[0] + (0.05-mu[0])*x[0]/0.2294*(-20*{x_norm_sub3}+1397./250) + (mu[1]-0.2794)*x[0]/{x_norm_sub3}"
# mapping_y_sub3 = f"x[1] + (0.05-mu[0])*x[1]/0.2294*(-20*{x_norm_sub3}+1397./250) + (mu[1]-0.2794)*x[1]/{x_norm_sub3}"
# mapping_z_sub3 = f"{z_sub3} + (0.05-mu[0])*{z_sub3}/0.2294*(-20*{x_norm_sub3}+1397./250) + (mu[1]-0.2794)*{z_sub3}/{x_norm_sub3} + 1.5131"

# # For nozzle (subdomain 4)
# z_sub4 = "(x[2] - 1.5131)"
# x_norm_sub4 = "sqrt(x[0]**2 + x[1]**2 + (x[2] - 1.5131)**2)"
# mapping_x_sub4 = f"x[0] + (mu[1] - 0.2794)*x[0]/{x_norm_sub4}"
# mapping_y_sub4 = f"x[1] + (mu[1] - 0.2794)*x[1]/{x_norm_sub4}"
# mapping_z_sub4 = f"{z_sub4} + (mu[1] - 0.2794)*{z_sub4}/{x_norm_sub4} + 1.5131"

# @DEIM()
# @PullBackFormsToReferenceDomain()
# @ShapeParametrization(
#     (mapping_x_sub1, mapping_y_sub1, mapping_z_sub1),
#     (mapping_x_sub2, mapping_y_sub2, mapping_z_sub2),
#     (mapping_x_sub3, mapping_y_sub3, mapping_z_sub3),
#     (mapping_x_sub4, mapping_y_sub4, mapping_z_sub4),
# )
class ElasticBlock(EllipticCoerciveProblem):

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
        self.n = FacetNormal(mesh)
        self.dx = Measure("dx")(subdomain_data=self.subdomains)
        self.ds = Measure("ds")(subdomain_data=self.boundaries)
        # ...
        self.E = 201.0*GPa
        self.nu = 0.306
        self.lambda_1 = self.E * self.nu / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu))
        self.lambda_2 = self.E / (2.0 * (1.0 + self.nu))
        # self.lambda_1 = ParametrizedExpression(self, "mu[0]*pow(10,9) * mu[1] / ((1.0 + mu[1]) * (1.0 - 2.0 * mu[1]))", mu=(0.,0.), element=V.ufl_element())
        # self.lambda_2 = ParametrizedExpression(self, "mu[0]*pow(10,9) / (2.0 * (1.0 + mu[1]))", mu=(0.,0.), element=V.ufl_element())

    # Return custom problem name
    def name(self):
        return "HydrogenTank"

    # Return the alpha lower bound.
    def get_stability_factor_lower_bound(self):
        return 1.0

    # Return theta multiplicative terms of the affine expansion of the problem.
    def compute_theta(self, term):
        mu = self.mu
        if term == "a":
            theta_a0 = mu[0]
            theta_a1 = mu[0]
            theta_a2 = mu[0]
            theta_a3 = mu[0]
            return (theta_a0, theta_a1, theta_a2, theta_a3)
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
            a0 = self.elasticity(u, v) * dx(1)
            a1 = self.elasticity(u, v) * dx(2)
            a2 = self.elasticity(u, v) * dx(3)
            a3 = self.elasticity(u, v) * dx(4)
            return (a0, a1, a2, a3)
        elif term == "f":
            ds = self.ds
            f = -40.0*MPa
            n = self.n
            f0 = inner(f*n, v) * ds(1)
            return (f0, )
        elif term == "dirichlet_bc":
            bc0 = [DirichletBC(self.V, Constant((0.0, 0.0, 0.0)), self.boundaries, 2)]
            return (bc0,)
        elif term == "inner_product":
            u = self.u
            x0 = inner(u, v) * dx + inner(grad(u), grad(v)) * dx
            return (x0, )
        else:
            raise ValueError("Invalid term for assemble_operator().")

    # Auxiliary function to compute the elasticity bilinear form
    def elasticity(self, u, v):
        lambda_1 = self.lambda_1
        lambda_2 = self.lambda_2
        return 2.0 * lambda_2 * inner(sym(grad(u)), sym(grad(v))) + lambda_1 * tr(sym(grad(u))) * tr(sym(grad(v)))

mesh = Mesh("data/full_model.xml")
subdomains = MeshFunction("size_t", mesh, "data/full_model_physical_region.xml")
boundaries = MeshFunction("size_t", mesh, "data/full_model_facet_region.xml")

V = VectorFunctionSpace(mesh, "Lagrange", 2)

problem = ElasticBlock(V, subdomains=subdomains, boundaries=boundaries)
mu_range = [(0.5, 1.5), (1.0, 1.0)]
problem.set_mu_range(mu_range)

reduction_method = ReducedBasis(problem)
reduction_method.set_Nmax(30)
reduction_method.set_tolerance(1e-5)

solve_stage = "offline"
reduction_method.initialize_training_set(100)

# Delete old offline data
clean = Cleaner(problem.name())
names = clean.list_names(clean.path)
clean.delete_offline_folder() if solve_stage == "offline" and problem.name() in names else None

reduced_problem = reduction_method.offline()

online_mu = (1.0, 1.0)
reduced_problem.set_mu(online_mu)
reduced_solution = reduced_problem.solve()  # Obtain coefficient "weights"
reduced_problem.export_solution(filename="online_solution")
# P = plot(reduced_solution, reduced_problem=reduced_problem)

stop_time = timeit.default_timer()
run_time = round(stop_time - start_time, 2)

print(f'FINISH SOLVING IN {run_time}s')