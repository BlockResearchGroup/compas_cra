

import compas
import compas_cra
import os

from compas_assembly_old.datastructures import Assembly as AssemblyOld
from compas_cra.datastructures import CRA_Assembly

from compas_assembly.datastructures import Assembly

FILE_I = './armadillo_cra.json'
assembly_old = compas.json_load(os.path.join(compas_cra.DATA, './old/', FILE_I))
assembly_old = assembly_old.copy(cls=AssemblyOld)

print(assembly_old)

# assembly = Assembly()
assembly = CRA_Assembly()

for node in assembly_old.nodes():
    block = assembly_old.node_attribute(node, "block")
    is_support = assembly_old.node_attribute(node, "is_support")
    assembly.add_block(block, node=node, is_support=is_support)

for edge in assembly_old.edges():
    interface = assembly_old.edge_attribute(edge, "interface")
    assembly.add_interface_to_interfaces(edge[0], edge[1], interface)

compas.json_dump(assembly, os.path.join(compas_cra.DATA, FILE_I))
print(assembly)




