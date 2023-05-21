from src.quantum.error_correction.codes.simple_toric_code import RepeatingGraph
from src.quantum.error_correction.decoders.union_find.uf_functions import generate_spanning_trees, syndrome_validation_naive, tree_peeler

graph = RepeatingGraph(6)

# print graph without marked edges or vertices
graph.draw_graph()

edges = [5, 32, 33]
marked_vertices = graph.mark_vertices(edges)

# print graph without marked edges or vertices
graph.draw_graph(edges, marked_vertices)

clusters, count = syndrome_validation_naive(marked_vertices, graph)
spanning_trees = generate_spanning_trees(clusters, graph, marked_vertices)
for root, tree in spanning_trees.items():
    print(tree_peeler(list(tree.values())))
    print()
