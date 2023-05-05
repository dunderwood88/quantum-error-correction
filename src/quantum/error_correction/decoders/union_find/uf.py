from typing import Dict, List, Union

from src.quantum.error_correction.helpers import convert_binary_to_qubit_list, convert_qubit_list_to_binary

def syndrome_validation_naive(
    syndrome_string: Union[int, List[int]],
    code,
    syndrome_type: str="z"
):
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
    odd_clusters = {}

    if isinstance(syndrome_string, List):
        syndrome_string = convert_qubit_list_to_binary(syndrome_string)

    syn = syndrome_string
    full_step = 0
    first_pass = 1
    count = 0

    while syndrome_string > 0:
        for stabilizer_index, stabilizer in enumerate(
            code._stabilizers[syndrome_type]
        ):
            if (syn & 1) and stabilizer_index not in even_clusters:
                if not full_step:
                    data_tree, syndrome_tree = odd_clusters.get(
                        stabilizer_index, (0, 0)
                    )
                    if first_pass:
                        data_tree |= stabilizer
                    else:
                        for stab in convert_binary_to_qubit_list(syndrome_tree):
                            data_tree |= code._stabilizers[syndrome_type][stab]
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
                for root, trees in odd_clusters.items():
                    if root != stabilizer_index:
                        if bin(
                            odd_clusters[stabilizer_index][full_step] &
                            trees[full_step]
                        ).count("1") >= 1:
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
                    syndrome_string &= ~(
                        (1 << stabilizer_index) | (1 << root_merge)
                    )

            syn >>= 1

        if syn == 0:
            count += 1
            if odd_clusters:
                syn = syndrome_string
                full_step = not full_step

    return even_clusters, count



        

    return even_clusters
