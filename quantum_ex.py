from collections import defaultdict
from typing import List, Union

from quantum.error_correction.helpers import convert_qubit_list_to_binary
from quantum.error_correction.codes.rotated_planar import RPlanarCode


# code = RPlanarCode(7)
# # code.draw(x_error_string=69793742848, z_error_string=4294967552)
# # code.draw(x_error_string=[19, 30, 36], z_error_string=[32, 8])
# code.draw(z_error_string=[1, 47, 43], x_error_string=[7, 21, 28, 20, 34, 41])

code = RPlanarCode(7)
# code.draw(
#     x_error_string=[20, 42, 50, 33],
#     z_syndrome_string=[11],
#     restrict_graph="z"
# )

# x_error = [20, 33, 8]
# x_error_bin = convert_qubit_list_to_binary(x_error)
# syn = code.generate_syndrome(x_error_bin, error_type="x")

# print(syn)

# code.draw(
#     x_error_string=x_error_bin,
#     restrict_graph="z"
# )

# i = 0
# while syn:
#     if syn & 1:
#         x_error_bin |= code._stabilizers["z"][i]
#     syn >>= 1
#     i += 1

# code.draw(
#     x_error_string=x_error_bin,
#     restrict_graph="z"
# )


# growth
# 1: half-step
# for each stabilizer:
#   data_string |= stabilizer
#   store neighbouring stabilizers for step 2
#
# 2: full-step
# if cluster still odd
# for each neighbour stabilizer_index:
#   syndrome_string |= stabilizer_index


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
    
    # half step
    syn = syndrome_string
    flag = 0
    i = 0
    while True:
        if syn & 1:
            if not flag:
                stabilizer = code._stabilizers[syndrome_type][i]
                data_tree, syndrome_tree = odd_clusters[i]
                odd_clusters[i] = (data_tree | stabilizer, syndrome_tree)
            else:
                data_tree, syndrome_tree = odd_clusters[i]
                update_syndrome = code.generate_syndrome(data_tree)
                odd_clusters[i] = (data_tree, syndrome_tree | update_syndrome)

            root_merge = 0
            for root, trees in odd_clusters.items():
                if get_set_count(odd_clusters[i][flag] & trees[flag]) == 1:
                    root_merge = root
                    break
            
            if root_merge:
                data_tree, syndrome_tree = odd_clusters[root_merge]
                even_clusters[root_merge] = (
                    odd_clusters[i][0] | data_tree,
                    syndrome_tree | odd_clusters[i][1] | (1 << i) | (1 << root_merge)
                )
                del odd_clusters[i]
                del odd_clusters[root_merge]

        syn >>= 1

        if syn == 0:
            if odd_clusters:
                syn = syndrome_string
                flag = not flag
                i = 0
            else:
                break
        else:
            i += 1

    return even_clusters


# x_error = [20, 33, 8]
# x_error = [8, 16]
x_error = [26, 32]
x_error_bin = convert_qubit_list_to_binary(x_error)
syn = code.generate_syndrome(x_error_bin, error_type="x")


print("BEFORE UNION FIND")
code.draw(
    x_error_string=x_error_bin,
    restrict_graph="z"
)


clusters = union_find_naive(syn)

data = 0
syn = 0

for i in clusters.values():
    data |= i[0]
    syn |= i[1]

print("AFTER UNION FIND")
print(clusters)
code.draw(
    x_error_string=data,
    z_syndrome_string=syn,
    restrict_graph="z"
)