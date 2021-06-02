from dolfin import *
from rbnics import *
import timeit
from tools.cleaner import Cleaner


start_time = timeit.default_timer()

# Define units
GPa = 10**9
MPa = 10**6

# Define parameters
#* Young's modulus = mu[0], Possion ratio = mu[1], Pressure = mu[2]


class PressureVessel(EllipticCoerciveProblem):

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

    # Return custom problem name
    def name(self):
        return "PressureVessel"

    # Return the alpha lower bound.
    def get_stability_factor_lower_bound(self):
        return 1.0

    # Return theta multiplicative terms of the affine expansion of the problem.
    def compute_theta(self, term):
        mu = self.mu
        if term == "a":
            theta_a0 = (mu[0]*GPa)*mu[1] / ((1 + mu[1])*(1 - 2*mu[1]))
            theta_a1 = (mu[0]*GPa) / (2*(1 + mu[1]))
            return (theta_a0, theta_a1)
        elif term == "f":
            theta_f0 = mu[2]
            return (theta_f0, )
        else:
            raise ValueError("Invalid term for compute_theta().")

    # Return forms resulting from the discretization of the affine expansion of the problem operators.
    def assemble_operator(self, term):
        v = self.v
        dx = self.dx
        if term == "a":
            u = self.u
            a0 = tr(sym(grad(u))) * tr(sym(grad(v))) * dx
            a1 = 2 * inner(sym(grad(u)), sym(grad(v))) * dx
            return (a0, a1)
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


#* 4. Main program
###* 4.1 Read the mesh for this problem
mesh = Mesh("data/full_model.xml")
subdomains = MeshFunction("size_t", mesh, "data/full_model_physical_region.xml")
boundaries = MeshFunction("size_t", mesh, "data/full_model_facet_region.xml")

###* 4.2 Create Finite Element space (Lagrange P2)
V = VectorFunctionSpace(mesh, "Lagrange", 2)

###* 4.3 Allocate an object of the PressureVessel class
problem = PressureVessel(V, subdomains=subdomains, boundaries=boundaries)
mu_range = [(200., 210.), (0.15, 0.35), (0.5, 1.5)]
problem.set_mu_range(mu_range)

###* 4.4 Prepare reduction with a reduced basis method
reduction_method = ReducedBasis(problem)
reduction_method.set_Nmax(30)
reduction_method.set_tolerance(1e-5)

###* 4.5 Perform the offline phase
def retrain_offline_phase(problem, Cleaner, is_offline=bool):
    assert isinstance(is_offline, bool), "is_offline must be True or False"
    foldername = problem.name()
    clean = Cleaner(foldername)
    names = clean.get_names(clean.path)
    if is_offline and foldername in names:
        clean.delete_offline_folder()

retrain_offline_phase(problem, Cleaner, True)
reduction_method.initialize_training_set(100)
reduced_problem = reduction_method.offline()

stop_time_offline = timeit.default_timer()

###* 4.6 Perform an online solve
online_mu = (201.0, 0.306, 1.0)
reduced_problem.set_mu(online_mu)
reduced_solution = reduced_problem.solve()
reduced_problem.export_solution(filename="online_solution")

stop_time_online = timeit.default_timer()

###* 4.7 Post-processing
# Get u (\curly{N}) solution from the reduced solution (N)
N_max = reduced_solution.N
basis_functions = reduced_problem.basis_functions[:N_max]
u = basis_functions * reduced_solution

# Get lambda_1 and lambda_2 from the theta terms
lambda_ = reduced_problem.compute_theta("a")

def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)

def sigma(u, lambda_):
    dim = u.geometric_dimension()  # domain dimension
    return lambda_[0]*tr(epsilon(u))*Identity(dim) + 2.0*lambda_[1]*epsilon(u)

run_time_offline = round(stop_time_offline - start_time, 2)
run_time_online = round(stop_time_online - stop_time_offline, 2)
print(f'Finish solving offline in {run_time_offline}s')
print(f'Finish solving online in {run_time_online}s')

# Compute stresses
dim = u.geometric_dimension()
s = sigma(u, lambda_) - (1./3)*tr(sigma(u, lambda_))*Identity(dim)
von_Mises = sqrt(3./2*inner(s, s))
V = FunctionSpace(mesh, 'Lagrange', 2)
von_Mises = project(von_Mises, V)
von_Mises.rename("von_Mises", "Von Mises")

# Export stress results
file_results = XDMFFile(problem.name() + "/stress_results.xdmf")
file_results.parameters["flush_output"] = True
file_results.parameters["functions_share_mesh"] = True
file_results.write(von_Mises, 0.)
