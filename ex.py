from graphs.graph import Graph
from graphs.basic_search.bfs import breadth_first_search

g = Graph()
g.add_edge(0, 1)
g.add_edge(1, 0)
g.add_edge(1, 10)
g.add_edge(1, 10)
g.add_edge(1, 11)
g.add_edge(1, 12)
g.add_edge(0, 2)
g.add_edge(2, 4)
g.add_edge(2, 5)
g.add_edge(2, 6)
g.add_edge(0, 3)
g.add_edge(3, 7)
g.add_edge(3, 8)
g.add_edge(3, 9)

print(g.get_definition())

path = breadth_first_search(g, 1)

print(path)