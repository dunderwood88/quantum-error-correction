from src.quantum.error_correction.helpers import convert_qubit_list_to_binary
from src.quantum.error_correction.codes.toric_code import ToricCode

from src.quantum.error_correction.decoders.union_find.uf_functions import (
    generate_spanning_trees, syndrome_validation_naive, tree_peeler)

### --- D = 3  --- ###
# code = ToricCode(3)
# # x_error = [1, 14, 17]
# x_error = [1]
### -------------- ###

### --- D = 5  --- ###
code = ToricCode(5)
x_error = [1, 14, 17]
x_error = [1, 15]
x_error = [1, 7, 12]
x_error = [3, 13, 43]
### -------------- ###

## --- D = 7  --- ###
# code = ToricCode(7)
# x_error = [2, 10, 17]
# x_error = [4, 18, 32, 46, 60]
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
# code = RPlanarCode(15)
# x_error = [73, 85, 97, 109, 108, 94]
# x_error = [123, 124]
# x_error = [33, 49]
# x_error = [32, 48]
# x_error = [65, 81, 97, 113]
# x_error = [96, 97]
# x_error = [111, 112]
# x_error = [65, 81, 97]
### -------------- ###


# 1) get the initial syndrome
x_error_bin = convert_qubit_list_to_binary(x_error)
error_syn = code.generate_syndrome(x_error_bin, error_type="x")
print(f"Error syndrome is {error_syn}")
print()

# print the code with syndrome
print("BEFORE UNION FIND")
code.draw(
    x_data_string=x_error_bin,
    restrict_graph="z"
)
# exit()

# 2) grow clusters
clusters, count = syndrome_validation_naive(error_syn, code)

data = 0
syn = 0
for i in clusters.values():
    data |= i[0]
    syn |= i[1]


# print the clusters
print("AFTER UNION FIND")
print(count)
print(clusters)
code.draw(
    x_data_string=data,
    z_syndrome_string=syn,
    restrict_graph="z"
)


# 3) create spanning tree for each cluster
spanning_trees = generate_spanning_trees(clusters, code, error_syn)

print(spanning_trees)
print()

# 4) peeling stage to find corrections
for root, tree in spanning_trees.items():
    print(tree_peeler(list(tree.values())))
    print()
