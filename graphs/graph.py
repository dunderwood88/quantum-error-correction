from collections import defaultdict


class Graph():
    """Class to represent the structure of a graph
    """

    def __init__(self) -> None:
        self._vertex_map = defaultdict(list)
        self._num_vertices = 0

    def add_edge(self, i_vertex: int, f_vertex: int):
        """_summary_

        Parameters
        ----------
        i_vertex : int
            _description_
        f_vertex : int
            _description_
        """

        self._vertex_map[i_vertex].append(f_vertex)
        n = max(i_vertex, f_vertex)
        if n > self._num_vertices:
            self._num_vertices = n + 1

    def get_definition(self):
        return self._vertex_map

    def get_num_vertices(self):
        return self._num_vertices