from dolfin import *
from rbnics import *
from subfiles.clear_offline_data import Cleaner


# Define parameters for mapping functions
offset = "1.5131"
#* thickness = mu[0], radius = mu[1]

# For cylindrical body (subdomain 1)
x_norm_sub1 = "sqrt(x[0]**2 + x[1]**2)"
mapping_x_sub1 = f"x[0] + (0.05-mu[0])*x[0]/0.2294*(-20*{x_norm_sub1}+1397./250) + (mu[1]-0.2794)*x[0]/{x_norm_sub1}"
mapping_y_sub1 = f"x[1] + (0.05-mu[0])*x[1]/0.2294*(-20*{x_norm_sub1}+1397./250) + (mu[1]-0.2794)*x[1]/{x_norm_sub1}"
mapping_z_sub1 = "x[2]"

# For left spherical head (subdomain 2)
z_sub2 = f"x[2] + {offset}"
x_norm_sub2 = f"sqrt(x[0]**2 + x[1]**2 + {z_sub2}**2)"
mapping_x_sub2 = f"x[0] + (0.05-mu[0])*x[0]/0.2294*(-20*{x_norm_sub2}+1397./250) + (mu[1]-0.2794)*x[0]/{x_norm_sub2}"
mapping_y_sub2 = f"x[1] + (0.05-mu[0])*x[1]/0.2294*(-20*{x_norm_sub2}+1397./250) + (mu[1]-0.2794)*x[1]/{x_norm_sub2}"
mapping_z_sub2 = f"{z_sub2} + (0.05-mu[0])*{z_sub2}/0.2294*(-20*{x_norm_sub2}+1397./250) + (mu[1]-0.2794)*{z_sub2}/{x_norm_sub2} - {offset}"

# For right spherical head (subdomain 3)
z_sub3 = f"x[2] - {offset}"
x_norm_sub3 = f"sqrt(x[0]**2 + x[1]**2 + {z_sub3}**2)"
mapping_x_sub3 = f"x[0] + (0.05-mu[0])*x[0]/0.2294*(-20*{x_norm_sub3}+1397./250) + (mu[1]-0.2794)*x[0]/{x_norm_sub3}"
mapping_y_sub3 = f"x[1] + (0.05-mu[0])*x[1]/0.2294*(-20*{x_norm_sub3}+1397./250) + (mu[1]-0.2794)*x[1]/{x_norm_sub3}"
mapping_z_sub3 = f"{z_sub3} + (0.05-mu[0])*{z_sub3}/0.2294*(-20*{x_norm_sub3}+1397./250) + (mu[1]-0.2794)*{z_sub3}/{x_norm_sub3} + {offset}"

# For nozzle (subdomain 4)
z_sub4 = f"x[2] - {offset}"
x_norm_sub4 = f"sqrt(x[0]**2 + x[1]**2 + {z_sub4}**2)"
mapping_x_sub4 = f"x[0] + (mu[1] - 0.2794)*x[0]/{x_norm_sub4}"
mapping_y_sub4 = f"x[1] + (mu[1] - 0.2794)*x[1]/{x_norm_sub4}"
mapping_z_sub4 = f"{z_sub4} + (mu[1] - 0.2794)*{z_sub4}/{x_norm_sub4} + {offset}"

print("Mapping functions: done")

@DEIM()
@PullBackFormsToReferenceDomain()
@ShapeParametrization(
    (mapping_x_sub1, mapping_y_sub1, mapping_z_sub1),
    (mapping_x_sub2, mapping_y_sub2, mapping_z_sub2),
    (mapping_x_sub3, mapping_y_sub3, mapping_z_sub3),
    (mapping_x_sub4, mapping_y_sub4, mapping_z_sub4),
)
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
        self.dx = Measure("dx")(subdomain_data=self.subdomains)
        self.ds = Measure("ds")(subdomain_data=self.boundaries)
        # ...
        self.f = Constant((0.0, 1.0, 0.0))
        self.E = 1.0
        self.nu = 0.3
        self.lambda_1 = self.E * self.nu / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu))
        self.lambda_2 = self.E / (2.0 * (1.0 + self.nu))

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
            theta_a0 = 1.0
            theta_a1 = 1.0
            theta_a2 = 1.0
            theta_a3 = 1.0
            return (theta_a0, theta_a1, theta_a2, theta_a3)
        elif term == "f":
            theta_f0 = 10.0
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
            f = self.f
            f0 = inner(f, v) * ds(2)
            return (f0, )
        elif term == "dirichlet_bc":
            bc0 = [DirichletBC(self.V, Constant((0.0, 0.0, 0.0)), self.boundaries, 3)]
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

mesh = Mesh("data/vessel.xml")
subdomains = MeshFunction("size_t", mesh, "data/vessel_physical_region.xml")
boundaries = MeshFunction("size_t", mesh, "data/vessel_facet_region.xml")

V = VectorFunctionSpace(mesh, "Lagrange", 1)

problem = ElasticBlock(V, subdomains=subdomains, boundaries=boundaries)

print("Create instance: done")

mu_range = [(0.04, 0.06), (0.27, 0.28)]
problem.set_mu_range(mu_range)

reduction_method = ReducedBasis(problem)

print("Set training parameters")

reduction_method.set_Nmax(30, DEIM=20)
reduction_method.set_tolerance(1e-5, DEIM=1e-3)

print("Create training set")
solve_stage = "offline"
reduction_method.initialize_training_set(100, DEIM=25)
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

print("Begin offline training")
reduced_problem = reduction_method.offline()

print("Offline training: done")

online_mu = (0.04, 0.2794)
reduced_problem.set_mu(online_mu)
reduced_solution = reduced_problem.solve()  # Obtain coefficient "weights"
reduced_problem.export_solution(filename="online_solution")
# P = plot(reduced_solution, reduced_problem=reduced_problem)
