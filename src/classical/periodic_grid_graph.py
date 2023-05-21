from typing import List

from src.classical.helpers import convert_binary_to_qubit_list
from src.quantum.codes.toric_code import ToricCode


class PeriodicGridGraph(ToricCode):
    """Implements a graph of vertices and edges of a certain dimension which
    has periodic boundary conditions - left/right and top/bottom edges are
    identical.
    Provides a simplified interface to the Toric Code which helps describe the
    problem structure without QEC language.

    Defined as coordinates over:
        e: edges
        V: vertices
    Indices run from left-to-right, top-to-bottom, for both edges/vertices.

    Example for dimension = 5

        e0      e1      e2      e3      e4

    e5  V0  e6  V1  e7  V2  e8  V3  e9  V4  e5

        e10     e11     e12     e13     e14

    e15 V5  e16 V6  e17 V7  e18 V8  e19 V9  e15

        e20     e21     e22     e23     e24

    e25 V10 e26 V11 e27 V12 e28 V13 e29 V14 e25

        e30     e31     e32     e33     e34

    e35 V15 e36 V16 e37 V17 e38 V18 e39 V19 e35

        e40     e41     e42     e43     e44

    e45 V20 e46 V21 e47 V22 e48 V23 e49 V24 e45

        e0      e1      e2      e3      e4

    """

    def __init__(self, dimension: int) -> None:
        super().__init__(dimension)
        self._name = f"D = {dimension} periodic grid graph"

    def get_neighbour_edges(self, index: int) -> List[int]:
        """Gets the neighbouring edges for a given vertex index.

        Parameters
        ----------
        index : int
            the index of the vertex to get the neighbours of

        Returns
        -------
        List[int]
            list of the edge indexes
        """
        return convert_binary_to_qubit_list(self.get_stabilizers(index))

    def mark_vertices(self, edges: List[int]) -> List[int]:
        """Given a list of edges, marks adjecent vertices. Marking a vertex an
        even number of times reverts that vertex to an unmarked state.

        Parameters
        ----------
        edges : List[int]
            a list of edge indexes

        Returns
        -------
        List[int]
            a list of resulting marked vertices
        """
        return convert_binary_to_qubit_list(self.generate_syndrome(edges))

    def draw_graph(
        self,
        edges: List[int] = 0,
        vertices: List[int] = 0,
    ) -> None:
        """

        Parameters
        ----------
        edges : List[int], optional
            a list of edges to highlight, by default 0
        vertices : List[int], optional
            a list of vertices to highlight, by default 0
        """

        self.draw(
            x_data_string=edges,
            z_syndrome_string=vertices,
            restrict_graph="z",
            simplify=True
        )
