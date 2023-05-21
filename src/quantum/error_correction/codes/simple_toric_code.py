from typing import List
from src.quantum.error_correction.codes.toric_code import ToricCode
from src.quantum.error_correction.helpers import convert_binary_to_qubit_list


class RepeatingGraph(ToricCode):

    def __init__(self, dimension: int) -> None:
        super().__init__(dimension)

    def get_neighbour_edges(self, index: int):
        return self.get_stabilizers(index)

    def mark_vertices(self, edges: List[int]) -> List[int]:
        return convert_binary_to_qubit_list(self.generate_syndrome(edges))

    def draw_graph(
        self,
        edges: List[int] = 0,
        vertices: List[int] = 0,
    ) -> None:

        self.draw(
            x_data_string=edges,
            z_syndrome_string=vertices,
            restrict_graph="z",
            simplify=True
        )
