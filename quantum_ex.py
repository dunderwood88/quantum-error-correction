from collections import OrderedDict, defaultdict
from typing import List, Union

from src.quantum.error_correction.helpers import convert_qubit_list_to_binary, convert_binary_to_qubit_list
from src.quantum.error_correction.codes.rotated_planar import RPlanarCode

from src.quantum.error_correction.decoders.union_find.uf_functions import generate_spanning_trees, syndrome_validation_naive, tree_peeler

### --- D = 3  --- ###
code = RPlanarCode(3)
x_error = [4]
### -------------- ###

### --- D = 5  --- ###
# code = RPlanarCode(5)
# x_error = [12, 16]
### -------------- ###

## --- D = 7  --- ###
# code = RPlanarCode(7)
# x_error = [18, 24, 30]
# x_error = [23, 24, 18]
# x_error = [18, 24]
## -------------- ###

### --- D = 9  --- ###
# code = RPlanarCode(9)
# x_error = [18, 24, 30]
# x_error = [23, 24, 18]
# x_error = [30, 38, 39, 33]
### -------------- ###

### --- D = 11 --- ###
# code = RPlanarCode(11)
# x_error = [40, 50, 60]
### -------------- ###

### --- D = 15 --- ###
code = RPlanarCode(15)
# x_error = [73, 85, 97, 109, 108, 94]
# x_error = [123, 124]
x_error = [33, 49]
# x_error = [32, 48]
# x_error = [65, 81, 97, 113]
# x_error = [96, 97]
# x_error = [111, 112]
# x_error = [65, 81, 97]
### -------------- ###


x_error_bin = convert_qubit_list_to_binary(x_error)
error_syn = code.generate_syndrome(x_error_bin, error_type="x")
print(error_syn)


print("BEFORE UNION FIND")
code.draw(
    x_data_string=x_error_bin,
    restrict_graph="z"
)
# exit()

clusters, count = syndrome_validation_naive(error_syn, code)

data = 0
syn = 0

for i in clusters.values():
    data |= i[0]
    syn |= i[1]

# data = clusters[57][0]
# syn = clusters[57][1]

print("AFTER UNION FIND")
print(count)
print(clusters)
code.draw(
    x_data_string=data,
    z_syndrome_string=syn,
    restrict_graph="z"
)

print(f"Error syndrome: {convert_binary_to_qubit_list(error_syn)}")
print()


# def generate_spanning_trees(code, clusters):

#     spanning_trees = {}

#     syn = error_syn

#     for root, (data_qubits, syndrome_qubits) in clusters.items():

#         syn_list = convert_binary_to_qubit_list(syndrome_qubits)

#         stack = [(syn_list[0],)]
#         visited = OrderedDict()

#         while stack:

#             # get top item in the stack
#             current_node = stack.pop()

#             if current_node[0] in visited:
#                 continue

#             # get the data qubits of the stabilizer
#             stab = code.get_stabilizers(current_node[0], "z")

#             # find adjacent nodes that are part of the syndrome
#             neighbours = [
#                 (
#                     n,
#                     current_node[0],
#                     convert_binary_to_qubit_list(
#                         stab & code.get_stabilizers(n, "z")
#                     )[0]
#                 ) for n in convert_binary_to_qubit_list(
#                     code.generate_syndrome(stab)
#                 ) if n in syn_list
#             ]

#             # add unvisted syndrome nodes to the top of the stack
#             stack.extend(neighbours)

#             visited[current_node[0]] = [
#                 current_node[2],  # edge
#                 1 if (1 << current_node[1]) & syn else 0,
#                 1 if (1 << current_node[0]) & syn else 0

#              ] if len(current_node) > 1 else None

#         spanning_trees[root] = visited
#         syn ^= syndrome_qubits

#     return spanning_trees


spanning_trees = generate_spanning_trees(clusters, code, error_syn)





for root, tree in spanning_trees.items():
    print(tree_peeler(list(tree.values())))
    print()
