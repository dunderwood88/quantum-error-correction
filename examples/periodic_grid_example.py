from src.classical.periodic_grid_graph import PeriodicGridGraph
from src.classical.decoders.union_find.uf_functions import (
    generate_spanning_trees, grow_clusters, tree_peeler)

graph = PeriodicGridGraph(6)

# print graph without marked edges or vertices
graph.draw_graph()

edges = [5, 32, 33]
marked_vertices = graph.mark_vertices(edges)

# print graph without marked edges or vertices
graph.draw_graph(edges, marked_vertices)

clusters, count = grow_clusters(marked_vertices, graph)
spanning_trees = generate_spanning_trees(clusters, graph, marked_vertices)
for root, tree in spanning_trees.items():
    print(tree_peeler(list(tree.values())))
    print()
