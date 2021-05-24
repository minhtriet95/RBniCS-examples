from dolfin import *
from rbnics import *
from subfiles.clear_offline_data import Cleaner
import matplotlib.pyplot as plt

from matplotlib import rc
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)


@SCM()
@PullBackFormsToReferenceDomain()
@ShapeParametrization(
    ("2*mu[0]*x[0]", "2*mu[0]*x[1]"),  # subdomain 1
    ("2*mu[0]*x[0]", "2*mu[0]*x[1]"),  # subdomain 2
    ("2*mu[0]*x[0]", "2*mu[0]*x[1]"),  # subdomain 3
    ("2*mu[0]*x[0]", "2*mu[0]*x[1]"),  # subdomain 4
    ("2*mu[0]*x[0]", "2*mu[0]*x[1]"),  # subdomain 5
    ("2*mu[0]*x[0]", "2*mu[0]*x[1]"),  # subdomain 6
    ("2*mu[0]*x[0]", "2*mu[0]*x[1]"),  # subdomain 7
    ("2*mu[0]*x[0]", "2*mu[0]*x[1]"),  # subdomain 8
    ("2*mu[0]*x[0]", "2*mu[0]*x[1]"),  # subdomain 9
    ("(3*mu[0] - 1/2)*x[0] + (3 - 6*mu[0])*x[1] + (3/2 - 3*mu[0])", 
     "(mu[0] - 1/2)*x[0] + (3 - 4*mu[0])*x[1] + (3/2 - 3*mu[0])"),  # subdomain 10
    ("(3 - 4*mu[0])*x[0] + (mu[0] - 1/2)*x[1] + (3/2 - 3*mu[0])",
     "(3 - 6*mu[0])*x[0] + (3*mu[0] - 1/2)*x[1] + (3/2 - 3*mu[0])"),  # subdomain 11
    ("(3 - 4*mu[0])*x[0] + (1/2 - mu[0])*x[1] + (3/2 - 3*mu[0])",
     "(6*mu[0] - 3)*x[0] + (3*mu[0] - 1/2)*x[1] + (3*mu[0] - 3/2)"),  # subdomain 12
    ("(3*mu[0] - 1/2)*x[0] + (6*mu[0] - 3)*x[1] + (3/2 - 3*mu[0])",
     "(1/2 - mu[0])*x[0] + (3 - 4*mu[0])*x[1] + (3*mu[0] - 3/2)"),  # subdomain 13
    ("(3*mu[0] - 1/2)*x[0] + (3 - 6*mu[0])*x[1] + (3*mu[0] - 3/2)",
     "(mu[0] - 1/2)*x[0] + (3 - 4*mu[0])*x[1] + (3*mu[0] - 3/2)"),  # subdomain 14
    ("(3 - 4*mu[0])*x[0] + (mu[0] - 1/2)*x[1] + (3*mu[0] - 3/2)",
     "(3 - 6*mu[0])*x[0] + (3*mu[0] - 1/2)*x[1] + (3*mu[0] - 3/2)"),  # subdomain 15
    ("(3 - 4*mu[0])*x[0] + (1/2 - mu[0])*x[1] + (3*mu[0] - 3/2)",
     "(6*mu[0] - 3)*x[0] + (3*mu[0] - 1/2)*x[1] + (3/2 - 3*mu[0])"),  # subdomain 16
    ("(3*mu[0] - 1/2)*x[0] + (6*mu[0] - 3)*x[1] + (3*mu[0] - 3/2)",
     "(1/2 - mu[0])*x[0] + (3 - 4*mu[0])*x[1] + (3/2 - 3*mu[0])"),  # subdomain 17
    ("x[0]", "(2 - 2*mu[0])*x[1] + (1 - 2*mu[0])"),  # subdomain 18
    ("(2 - 2*mu[0])*x[0] + (1 - 2*mu[0])", "x[1]"),  # subdomain 19
    ("x[0]", "(2 - 2*mu[0])*x[1] + (2*mu[0] - 1)"),  # subdomain 20
    ("(2 - 2*mu[0])*x[0] + (2*mu[0] - 1)", "x[1]"),  # subdomain 21
)
class UnitCell(EllipticCoerciveProblem):

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
        # Customize eigen solver parameters to enhance the convergence (remove when not using SCM)
        # File location: /usr/include/dolfin/la/SLEPcEigenSolver.h
        self._eigen_solver_parameters.update({
            "bounding_box_minimum": {
                "solver": "krylov-schur", "problem_type": "gen_hermitian",
                "spectral_transform": "shift-and-invert", "spectral_shift": 1.e-5,
                "linear_solver": "mumps"
            },
            "bounding_box_maximum": {
                "solver": "krylov-schur", "problem_type": "gen_hermitian",
                "spectral_transform": "shift-and-invert", "spectral_shift": 1.e5,
                "linear_solver": "mumps"
            },
            "stability_factor": {
                "solver": "krylov-schur", "problem_type": "gen_hermitian",
                "spectral_transform": "shift-and-invert", "spectral_shift": 1.e-5,
                "linear_solver": "mumps"
            }
        })

    # Return custom problem name
    def name(self):
        return "CellTest"

    # Return theta multiplicative terms of the affine expansion of the problem.
    @compute_theta_for_stability_factor
    def compute_theta(self, term):
        mu = self.mu
        if term == "a":
            theta_a0 = 10.0
            theta_a1 = 1.0
            theta_a2 = 1.0
            theta_a3 = 1.0
            theta_a4 = 1.0
            theta_a5 = 1.0
            theta_a6 = 1.0
            theta_a7 = 1.0
            theta_a8 = 1.0
            theta_a9 = 1.0
            theta_a10 = 1.0
            theta_a11 = 1.0
            theta_a12 = 1.0
            theta_a13 = 1.0
            theta_a14 = 1.0
            theta_a15 = 1.0
            theta_a16 = 1.0
            theta_a17 = 1.0
            theta_a18 = 1.0
            theta_a19 = 1.0
            theta_a20 = 1.0
            return (theta_a0, theta_a1, theta_a2,
                    theta_a3, theta_a4, theta_a5,
                    theta_a6, theta_a7, theta_a8,
                    theta_a9, theta_a10, theta_a11,
                    theta_a12, theta_a13, theta_a14,
                    theta_a15, theta_a16, theta_a17,
                    theta_a18, theta_a19, theta_a20)
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
            a9 = inner(grad(u), grad(v)) * dx(10)
            a10 = inner(grad(u), grad(v)) * dx(11)
            a11 = inner(grad(u), grad(v)) * dx(12)
            a12 = inner(grad(u), grad(v)) * dx(13)
            a13 = inner(grad(u), grad(v)) * dx(14)
            a14 = inner(grad(u), grad(v)) * dx(15)
            a15 = inner(grad(u), grad(v)) * dx(16)
            a16 = inner(grad(u), grad(v)) * dx(17)
            a17 = inner(grad(u), grad(v)) * dx(18)
            a18 = inner(grad(u), grad(v)) * dx(19)
            a19 = inner(grad(u), grad(v)) * dx(20)
            a20 = inner(grad(u), grad(v)) * dx(21)
            return (a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20)
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
mesh = Mesh("data/cellAffine.xml")
subdomains = MeshFunction("size_t", mesh, "data/cellAffine_physical_region.xml")
boundaries = MeshFunction("size_t", mesh, "data/cellAffine_facet_region.xml")

# 2. Create Finite Element space for this problem
V = FunctionSpace(mesh, "Lagrange", 1)

# 3. Allocate an object for UnitCell class
problem = UnitCell(V, subdomains=subdomains, boundaries=boundaries)
mu_range = [(0.3, 0.6)]
problem.set_mu_range(mu_range)

# 4. Prepare reduction with a ReducedBasis method
reduction_method = ReducedBasis(problem)
reduction_method.set_Nmax(50, SCM=20)
reduction_method.set_tolerance(1e-4, SCM=0.75)

# 5. Perform the offline phase
solve_stage = "offline"
reduction_method.initialize_training_set(50, SCM=200)
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
plt.savefig('rb_solution_max.pdf')
plt.show()


# 8. Perform an error analysis
# reduction_method.initialize_testing_set(100, EIM=100)
# reduction_method.error_analysis()

# reduction_method.initialize_testing_set(100, EIM=100)
# reduction_method.speedup_analysis()
