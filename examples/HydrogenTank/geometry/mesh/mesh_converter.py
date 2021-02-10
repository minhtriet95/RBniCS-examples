from dolfin import *
from mshr import *
import os
import meshio
import numpy as np


# Extract file extension
filePath = 'vessel.msh'
fileName, fileExtension = os.path.splitext(filePath)

# Import mesh
msh = meshio.read(filePath)

dim = None
elements = dict()
keys = list()  # Element types
for item in msh.cells:
    elements[str(item.type)] = len(item.data)  # Store element info into dict()
    
    #! Only support 'line', 'triangle', and 'tetra'
    keys.append(item.type)
    if item.type == 'tetra':
        dim = 3  # 3-Dimentional case
    else:
        dim = 2  # 2-Dimentional case

points = msh.points
if dim == 2: points = np.delete(points, 2, axis=1)  # Remove z column

dash = '-' * 40
print(dash)
print('{:^40}'.format('MESH INFO'))
print(dash)

print('{:<25}{:>15}'.format('Case:', str(dim) + '-D'))
print('{:<25}{:>15}'.format('Number of nodes:', len(points)))
print('{:<25}'.format('Number of elements:'))
for key, value in elements.items():
    print('{:<15}{:<15}{:>10}'.format('', key, value))
    
print(msh.cell_data)
    
print('{:<25}{:>15}'.format('Number of subdomains:', len(msh.cells[1][1])))
print(dash)


print(type(msh.cell_data["gmsh:physical"][1]))

# cell_in_block_1 = msh.cells[0].data  # : list of subdomains with corresponding data (type, number_of_cell)
# cell_in_block_2 = msh.cells[1].data
# cell_in_block_3 = msh.cells[2].data
# cell_in_block_4 = msh.cells[3].data

# blocks = [cell_in_block_1, cell_in_block_2, cell_in_block_3, cell_in_block_4]

# count = 1
# block_data = list()  # Use list() for fast
# for block in blocks:  # Loop through all subdomains
#     for elem in block:  # Loop through all elements in current subdomain
#         block_data.append(count)
#     count += 1

# cell_data = np.array(block_data)

# Convert to xdmf files
meshio.write("xdmf/mesh.xdmf", meshio.Mesh(
    points, cells={keys[1]: msh.cells_dict[keys[1]]}))
meshio.write("xdmf/subdomains.xdmf", meshio.Mesh(
    points, cells={keys[1]: msh.cells_dict[keys[1]]},
    cell_data={"name_to_read": [msh.cell_data["gmsh:physical"][1]]}))
meshio.write("xdmf/boundaries.xdmf", meshio.Mesh(
    points, cells={keys[0]: msh.cells_dict[keys[0]]},
    cell_data={"name_to_read": [msh.cell_data["gmsh:physical"][0]]}))

# print(mesh.cell_data["gmsh:physical"])

# Reload xdmf files
mesh = Mesh()
with XDMFFile("xdmf/mesh.xdmf") as infile:
    infile.read(mesh)

# Physical (subdomain) regions
mvc_block = MeshValueCollection("size_t", mesh, dim)
with XDMFFile("xdmf/subdomains.xdmf") as infile:
    infile.read(mvc_block, "name_to_read")
subdomains = MeshFunction("size_t", mesh, mvc_block)

# Facet (boundary) regions
mvc_facet = MeshValueCollection("size_t", mesh, dim - 1)
with XDMFFile("xdmf/boundaries.xdmf") as infile:
    infile.read(mvc_facet, "name_to_read")
boundaries = MeshFunction("size_t", mesh, mvc_facet)

cell_data = msh.cell_data["gmsh:physical"][1]
for index, value in enumerate(mesh.cells()):
    mesh.domains().set_marker([index, cell_data[index]], dim)

sideset_1 = boundaries.where_equal(1)
sideset_2 = boundaries.where_equal(2)
sideset_3 = boundaries.where_equal(3)

boundaries.set_all(0)

for cell_index in sideset_1:
    boundaries.set_value(cell_index, 1)
for cell_index in sideset_2:
    boundaries.set_value(cell_index, 2)
for cell_index in sideset_3:
    boundaries.set_value(cell_index, 3)

File(fileName + ".xml") << mesh
File(fileName + "_physical_region.xml") << subdomains
File(fileName + "_facet_region.xml") << boundaries
