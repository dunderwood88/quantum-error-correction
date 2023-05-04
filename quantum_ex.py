from collections import defaultdict
from typing import List, Union

from src.quantum.error_correction.helpers import convert_qubit_list_to_binary, convert_binary_to_qubit_list
from src.quantum.error_correction.codes.rotated_planar import RPlanarCode


def get_set_count(num: int) -> int:

    count = 0
    while num:
        if num & 1:
            count += 1
        num >>= 1
    return count


def union_find_naive(syndrome_string: Union[int, List[int]], syndrome_type: str="z"):
    """
    growth
    1: half-step
    for each stabilizer:
      data_string |= stabilizer
      store neighbouring stabilizers for step 2

    # 2: full-step
    if cluster still odd
    for each neighbour stabilizer_index:
      syndrome_string |= stabilizer_index
    """

    # clusters
    # key: root syndrome qubit 
    # value: (data tree, syndrome tree)
    even_clusters = {}
    odd_clusters = defaultdict(lambda: (0, 0))
    
    syn = syndrome_string
    full_step = 0
    first_pass = 1
    stabilizer_index = 0
    count = 0
    while True:
        if stabilizer_index not in even_clusters.keys() and (syn & 1):
            if not full_step:
                data_tree, syndrome_tree = odd_clusters[stabilizer_index]
                if first_pass:
                    data_tree |= code._stabilizers[syndrome_type][stabilizer_index]
                else:
                    for stabilizer in convert_binary_to_qubit_list(syndrome_tree):
                        data_tree |= code._stabilizers[syndrome_type][stabilizer]
                odd_clusters[stabilizer_index] = (data_tree, syndrome_tree)
            else:
                data_tree, syndrome_tree = odd_clusters[stabilizer_index]
                update_syndrome = code.generate_syndrome(
                    data_tree, show_all_adjacent=True
                )
                odd_clusters[stabilizer_index] = (
                    data_tree, syndrome_tree | update_syndrome
                )
                first_pass = 0

            root_merge = 0
            for root, trees in {
                k: v for k, v in odd_clusters.items() if k != stabilizer_index
            }.items():
                if get_set_count(
                    odd_clusters[stabilizer_index][full_step] & trees[full_step]
                ) >= 1:
                    root_merge = root
                    break

            if root_merge:
                data_tree, syndrome_tree = odd_clusters[root_merge]
                even_clusters[root_merge] = (
                    odd_clusters[stabilizer_index][0] | data_tree,
                    syndrome_tree
                    | odd_clusters[stabilizer_index][1]
                    | (1 << stabilizer_index)
                    | (1 << root_merge)
                )
                del odd_clusters[stabilizer_index]
                del odd_clusters[root_merge]
                syndrome_string &= ~((1 << stabilizer_index) | (1 << root_merge))

        syn >>= 1

        if syn == 0:
            count += 1
            if odd_clusters:
                syn = syndrome_string
                full_step = not full_step
                stabilizer_index = 0
            else:
                break
        else:
            stabilizer_index += 1

    return even_clusters, count

### --- D = 7  --- ###
# code = RPlanarCode(7)
# x_error = [18, 24, 30]
# x_error = [23, 24, 18]
### -------------- ###

### --- D = 9  --- ###
code = RPlanarCode(9)
# x_error = [18, 24, 30]
# x_error = [23, 24, 18]
x_error = [30, 38, 39, 33]
### -------------- ###

### --- D = 11 --- ###
# code = RPlanarCode(11)
# x_error = [40, 50, 60]
### -------------- ###

### --- D = 15 --- ###
code = RPlanarCode(15)
x_error = [73, 85, 97, 109, 108, 94]
### -------------- ###


x_error_bin = convert_qubit_list_to_binary(x_error)
syn = code.generate_syndrome(x_error_bin, error_type="x")


print("BEFORE UNION FIND")
code.draw(
    x_error_string=x_error_bin,
    restrict_graph="z"
)
# exit()

clusters, count = union_find_naive(syn)

data = 0
syn = 0

for i in clusters.values():
    data |= i[0]
    syn |= i[1]

print("AFTER UNION FIND")
print(count)
print(clusters)
code.draw(
    x_error_string=data,
    z_syndrome_string=syn,
    restrict_graph="z"
)