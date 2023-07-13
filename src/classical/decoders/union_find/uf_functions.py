from collections import OrderedDict
from typing import Dict, List, Tuple, Union

from src.classical.helpers import (convert_binary_to_qubit_list,
                                   convert_qubit_list_to_binary)
from src.classical.periodic_grid_graph import PeriodicGridGraph
from src.quantum.codes.abstract_surface_code import AbstractSurfaceCode


def grow_clusters(
    marked_vertices: List[int],
    graph: PeriodicGridGraph,
) -> Tuple[Dict[int, Tuple[List[int], List[int]]], int]:
    """Grows and combines clusters from marked vertex roots until they are even.
    Provides a simplified interface to the syndrome_validation function which
    helps describe the problem structure without QEC language.


    Parameters
    ----------
    marked_vertices : List[int]
        list of indexes of marked vertices from which to grow clusters
    graph : PeriodicGridGraph
        the code over which the clusters are grown

    Returns
    -------
    Tuple[Dict[int, Tuple[List[int], List[int]]], int]
        the resulting clusters dictionary, with keys as cluster roots and
        values as edges/vertices within the cluster, and the total number of
        growth steps
    """

    return syndrome_validation_naive(marked_vertices, graph)


def syndrome_validation_naive(
    syndrome_string: Union[int, List[int]],
    code: AbstractSurfaceCode,
    syndrome_type: str = "z"
) -> Tuple[Dict[int, Tuple[int, int]], int]:
    """Naive implementation of the syndrome validation step as outlined in
    https://arxiv.org/pdf/1709.06218.pdf.

    Takes an error syndrome of a particular type (x or z) and performs
    cluster growth from each active syndrome root.
    Starting with all clusters odd (single active syndrome), performs the
    following:

    1) create a list of odd clusters
    2) while an odd cluster exists:
        3) for all odd cluster C:
            a) grow C by a half-edge
            b) if C meets another cluster, fuse and update parity (odd or even)
            c) if C is even, remove from odd cluster list
    4) returns the resulting list of even clusters

    Parameters
    ----------
    syndrome_string : Union[int, List[int]]
        the binary string representing the syndrome of active syndrome qubit
        indexes
    code : AbstractSurfaceCode
        the code over which the clusters are grown, specifies the graph layout
    syndrome_type : str, optional
        the type of syndrome qubits involved (x or z), by default "z"

    Returns
    -------
    Tuple[Dict[int, Tuple[int, int]], int]
        the resulting even clusters dictionary, with keys as cluster roots and
        values as edges/vertices within the cluster (syndrome and data qubits),
        and the total number of growth steps
    """

    # set up the odd and even clusters
    # key: cluster root syndrome qubit
    # value: (data edges, syndrome vertices)
    even_clusters = {}
    odd_clusters = {}

    if isinstance(syndrome_string, List):
        syndrome_string = convert_qubit_list_to_binary(syndrome_string)

    syn = syndrome_string
    full_step = 0  # start with a half step
    first_step = 1  # flag to indicate the initial full step growth
    count = 0  # counter for total number of growh steps

    while syndrome_string > 0:
        for stabilizer_index, stabilizer in enumerate(
            code.get_stabilizers(stabilizer_type=syndrome_type)
        ):  # while there are odd clusters
            if (syn & 1) and stabilizer_index not in even_clusters:
                if not full_step:  # grow by a half step
                    data_tree, syndrome_tree = odd_clusters.get(
                        stabilizer_index, (0, 0)
                    )
                    if first_step:
                        data_tree |= stabilizer
                    else:
                        for stab in convert_binary_to_qubit_list(syndrome_tree):
                            data_tree |= code.get_stabilizers(
                                stab, syndrome_type)
                    odd_clusters[stabilizer_index] = (data_tree, syndrome_tree)
                else:
                    data_tree, syndrome_tree = odd_clusters[stabilizer_index]
                    update_syndrome = code.generate_syndrome(
                        data_tree, show_all_adjacent=True
                    )
                    odd_clusters[stabilizer_index] = (
                        data_tree, syndrome_tree | update_syndrome
                    )
                    first_step = 0

                root_merge = None  # if cluster meets another...
                for root, trees in odd_clusters.items():
                    if root != stabilizer_index:
                        if bin(
                            odd_clusters[stabilizer_index][full_step] &
                            trees[full_step]
                        ).count("1") >= 1:
                            root_merge = root
                            break

                if root_merge is not None:  # ...fuse and remove from odd cluster list
                    root_data_tree, root_syndrome_tree = odd_clusters[root_merge]
                    other_data_tree, other_syndrome_tree = odd_clusters[stabilizer_index]

                    even_clusters[root_merge] = (
                        other_data_tree | root_data_tree,
                        other_syndrome_tree | root_syndrome_tree
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


def generate_spanning_trees(
    clusters: Dict[int, Tuple[int, int]],
    code: AbstractSurfaceCode,
    original_syndrome: int,
    original_syndrome_type: str = "z"
) -> Dict[int, List[int]]:
    """Generates the spanning trees of the provided even clusters, defined over
    a specific code and arising from a given original syndrome

    Parameters
    ----------
    code : AbstractSurfaceCode
        the code over which the clusters are defined, specifies the graph layout
    clusters : Dict[int, Tuple[int, int]]
        a dictionary of even clusters, with keys as cluster roots and values
        as edges/vertices within the cluster (syndrome and data qubits)
    original_syndrome : int
        the original syndrome that gave rise to the clusters
    original_syndrome_type : int
        the type (x or z) of the original syndrome

    Returns
    -------
    Dict[int, Dict[int, List[int]]]
        an dictionary of original cluster root syndrome qubits mapped to 
        directed spanning tree data which itself is a dictionary mapping
        visited syndrome qubit indexes to spanning tree edge data
    """

    spanning_trees = {}

    if isinstance(original_syndrome, List):
        original_syndrome = convert_qubit_list_to_binary(original_syndrome)

    syn = original_syndrome

    # for each cluster, find the spanning tree
    for root, (data_qubits, syndrome_qubits) in clusters.items():

        syn_list = convert_binary_to_qubit_list(syndrome_qubits)

        stack = [(syn_list[0],)]  # start at lowest-indexed syndrome qubit
        visited = OrderedDict()

        # depth first search across cluster edges
        while stack:

            # get top item in the stack
            current_node = stack.pop()

            if current_node[0] in visited:
                continue

            # get the data qubits of the stabilizer
            stab = code.get_stabilizers(
                current_node[0], original_syndrome_type)

            # find adjacent nodes that are part of the syndrome
            neighbours = [
                (
                    n,
                    current_node[0],
                    convert_binary_to_qubit_list(
                        stab & code.get_stabilizers(n, original_syndrome_type)
                    )[0]
                ) for n in convert_binary_to_qubit_list(
                    code.generate_syndrome(stab & data_qubits)
                ) if n in syn_list
            ]

            # add unvisted syndrome nodes to the top of the stack
            stack.extend(neighbours)

            visited[current_node[0]] = [
                current_node[1],
                current_node[0],
                current_node[2]  # edge

            ] if len(current_node) > 1 else None

        spanning_trees[root] = list(visited.values())
        syn ^= syndrome_qubits

    return spanning_trees


def tree_peeler(spanning_trees, original_syndrome):

    if isinstance(original_syndrome, List):
        original_syndrome = convert_qubit_list_to_binary(original_syndrome)

    corrections = {}
    syn = original_syndrome

    for root, spanning_tree in spanning_trees.items():

        corrections[root] = []

        while spanning_tree:

            edge = spanning_tree.pop()

            if not edge:
                break

            if (1 << edge[1]) & syn:
                corrections[root].append(edge[2])
                syn ^= (1 << edge[1]) | (1 << edge[0])

    return corrections
