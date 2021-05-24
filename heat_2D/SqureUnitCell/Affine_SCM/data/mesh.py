from dolfin import *
from mshr import *
import numpy
import meshio
import matplotlib.pyplot as plt
from rbnics.shape_parametrization.utils.symbolic import VerticesMappingIO


# Import mesh
mesh = meshio.read("cellAffine.msh")
# points = mesh.points
points_3d = mesh.points
points = numpy.delete(points_3d, 2, axis=1)  # Remove z column

# Convert to xdmf files
meshio.write("xdmf/mesh.xdmf", meshio.Mesh(
    points, cells={"triangle": mesh.cells_dict["triangle"]}))
meshio.write("xdmf/subdomains.xdmf", meshio.Mesh(
    points, cells={"triangle": mesh.cells_dict["triangle"]},
    cell_data={"name_to_read": [mesh.cell_data["gmsh:physical"][1]]}))

mesh = Mesh()
with XDMFFile("xdmf/mesh.xdmf") as infile:
    infile.read(mesh)
mvc = MeshValueCollection("size_t", mesh, 2)
with XDMFFile("xdmf/subdomains.xdmf") as infile:
    infile.read(mvc, "name_to_read")
subdomains = MeshFunction("size_t", mesh, mvc)


# Define vertices mappings of affine shape parametrization.
sqrt_0_125 = sqrt(0.125)
vertices_mappings = [
    {
        ("0", "0"): ("0", "0"),
        ("0", "0.5"): ("0", "mu[0]"),
        ("0.5", "0"): ("mu[0]", "0")
    },  # subdomain 1
    {
        (f"-{sqrt_0_125}", f"-{sqrt_0_125}"): ("-mu[0]/sqrt(2)", "-mu[0]/sqrt(2)"),
        ("-0.6", "-0.6"): ("(-6*mu[0])/(5*sqrt(2))", "(-6*mu[0])/(5*sqrt(2))"),
        ("0", "-0.5"): ("0", "-mu[0]")
    },  # subdomain 2
    {
        ("-0.5", "0"): ("-mu[0]", "0"),
        ("-0.6", "-0.6"): ("(-6*mu[0])/(5*sqrt(2))", "(-6*mu[0])/(5*sqrt(2))"),
        (f"-{sqrt_0_125}", f"-{sqrt_0_125}"): ("-mu[0]/sqrt(2)", "-mu[0]/sqrt(2)")
    },  # subdomain 3
    {
        ("-0.6", "0.6"): ("(-6*mu[0])/(5*sqrt(2))", "(6*mu[0])/(5*sqrt(2))"),
        ("-0.5", "0"): ("-mu[0]", "0"),
        (f"-{sqrt_0_125}", f"{sqrt_0_125}"): ("-mu[0]/sqrt(2)", "mu[0]/sqrt(2)")
    },  # subdomain 4
    {
        ("-0.6", "0.6"): ("(-6*mu[0])/(5*sqrt(2))", "(6*mu[0])/(5*sqrt(2))"),
        (f"-{sqrt_0_125}", f"{sqrt_0_125}"): ("-mu[0]/sqrt(2)", "mu[0]/sqrt(2)"),
        ("0", "0.5"): ("0", "mu[0]")
    },  # subdomain 5
    {
        ("0", "0.5"): ("0", "mu[0]"),
        (f"{sqrt_0_125}", f"{sqrt_0_125}"): ("mu[0]/sqrt(2)", "mu[0]/sqrt(2)"),
        ("0.6", "0.6"): ("(6*mu[0])/(5*sqrt(2))", "(6*mu[0])/(5*sqrt(2))")
    },  # subdomain 6
    {
        (f"{sqrt_0_125}", f"{sqrt_0_125}"): ("mu[0]/sqrt(2)", "mu[0]/sqrt(2)"),
        ("0.5", "0"): ("mu[0]", "0"),
        ("0.6", "0.6"): ("(6*mu[0])/(5*sqrt(2))", "(6*mu[0])/(5*sqrt(2))")
    },  # subdomain 7
    {
        ("0.5", "0"): ("mu[0]", "0"),
        (f"{sqrt_0_125}", f"-{sqrt_0_125}"): ("mu[0]/sqrt(2)", "-mu[0]/sqrt(2)"),
        ("0.6", "-0.6"): ("(6*mu[0])/(5*sqrt(2))", "(-6*mu[0])/(5*sqrt(2))")
    },  # subdomain 8
    {
        (f"{sqrt_0_125}", f"-{sqrt_0_125}"): ("mu[0]/sqrt(2)", "-mu[0]/sqrt(2)"),
        ("0", "-0.5"): ("0", "-mu[0]"),
        ("0.6", "-0.6"): ("(6*mu[0])/(5*sqrt(2))", "(-6*mu[0])/(5*sqrt(2))")
    },  # subdomain 9
    {
        ("-0.6", "-0.6"): ("(-6*mu[0])/(5*sqrt(2))", "(-6*mu[0])/(5*sqrt(2))"),
        ("-1", "-1"): ("-1", "-1"),
        ("0", "-0.5"): ("0", "-mu[0]")
    },  # subdomain 10
    {
        ("-0.5", "0"): ("-mu[0]", "0"),
        ("-1", "-1"): ("-1", "-1"),
        ("-0.6", "-0.6"): ("(-6*mu[0])/(5*sqrt(2))", "(-6*mu[0])/(5*sqrt(2))")
    },  # subdomain 11
    {
        ("-1", "1"): ("-1", "1"),
        ("-0.5", "0"): ("-mu[0]", "0"),
        ("-0.6", "0.6"): ("(-6*mu[0])/(5*sqrt(2))", "(6*mu[0])/(5*sqrt(2))")
    },  # subdomain 12
    {
        ("-1", "1"): ("-1", "1"),
        ("-0.6", "0.6"): ("(-6*mu[0])/(5*sqrt(2))", "(6*mu[0])/(5*sqrt(2))"),
        ("0", "0.5"): ("0", "mu[0]")
    },  # subdomain 13
    {
        ("0", "0.5"): ("0", "mu[0]"),
        ("0.6", "0.6"): ("(6*mu[0])/(5*sqrt(2))", "(6*mu[0])/(5*sqrt(2))"),
        ("1", "1"): ("1", "1")
    },  # subdomain 14
    {
        ("0.6", "0.6"): ("(6*mu[0])/(5*sqrt(2))", "(6*mu[0])/(5*sqrt(2))"),
        ("0.5", "0"): ("mu[0]", "0"),
        ("1", "1"): ("1", "1")
    },  # subdomain 15
    {
        ("0.5", "0"): ("mu[0]", "0"),
        ("1", "-1"): ("1", "-1"),
        ("0.6", "-0.6"): ("(6*mu[0])/(5*sqrt(2))", "(-6*mu[0])/(5*sqrt(2))")
    },  # subdomain 16
    {
        ("0", "-0.5"): ("0", "-mu[0]"),
        ("1", "-1"): ("1", "-1"),
        ("0.6", "-0.6"): ("(6*mu[0])/(5*sqrt(2))", "(-6*mu[0])/(5*sqrt(2))")
    },  # subdomain 17
    {
        ("0", "-0.5"): ("0", "-mu[0]"),
        ("-1", "-1"): ("-1", "-1"),
        ("1", "-1"): ("1", "-1")
    },  # subdomain 18
    {
        ("-1", "1"): ("-1", "1"),
        ("-1", "-1"): ("-1", "-1"),
        ("-0.5", "0"): ("-mu[0]", "0")
    },  # subdomain 19
    {
        ("-1", "1"): ("-1", "1"),
        ("0", "0.5"): ("0", "mu[0]"),
        ("1", "1"): ("1", "1")
    },  # subdomain 20
    {
        ("0.5", "0"): ("mu[0]", "0"),
        ("1", "-1"): ("1", "-1"),
        ("1", "1"): ("1", "1")
    }  # subdomain 21
]


# Create boundaries
class Left(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and abs(x[0] + 1.0) < DOLFIN_EPS


class Right(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and abs(x[0] - 1.0) < DOLFIN_EPS


class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and abs(x[1] + 1.0) < DOLFIN_EPS


class Top(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and abs(x[1] - 1.0) < DOLFIN_EPS


boundaries = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
boundaries.set_all(0)
bottom = Bottom()
bottom.mark(boundaries, 1)
left = Left()
left.mark(boundaries, 4)  ##########
right = Right()
right.mark(boundaries, 2)
top = Top()
top.mark(boundaries, 3)  ##########


VerticesMappingIO.save_file(vertices_mappings, ".", "cellAffine_vertices_mapping.vmp")
File("cellAffine.xml") << mesh
File("cellAffine_physical_region.xml") << subdomains
File("cellAffine_facet_region.xml") << boundaries

# plot(mesh)
# plt.show()
