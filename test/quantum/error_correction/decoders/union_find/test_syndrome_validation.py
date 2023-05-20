import json
import unittest

from src.quantum.error_correction.codes.rotated_planar_code import RPlanarCode
from src.quantum.error_correction.decoders.union_find.uf_functions import syndrome_validation_naive
from src.quantum.error_correction.helpers import convert_qubit_list_to_binary, convert_binary_to_qubit_list

class TestUnionFind(unittest.TestCase):

    def setUp(self) -> None:

        self._codes = {
            "rotated planar": RPlanarCode
        }

        return super().setUp()

    def test_syndrome_validation(self):

        with open(
            "test/quantum/error_correction/decoders/union_find/test_data/syndrome_validation.json",
            "r"
        ) as syndrome_cluster_file:

            syndrome_cluster_data = json.load(syndrome_cluster_file)

            for code in syndrome_cluster_data["codes"]:
                code_class = self._codes[code["code_type"]]

                for dimension in code["dimensions"]:
                    code_instance = code_class(dimension["size"])

                    for syndrome_cluster in dimension["syndrome_clusters"]:

                        uf_clusters, count = syndrome_validation_naive(
                            syndrome_cluster["syndrome"],
                            code_instance
                        )

                        for expected_cluster in syndrome_cluster["clusters"]:
                            self.assertIn(
                                expected_cluster["root"], uf_clusters.keys()
                            )
                            uf_cluster = uf_clusters[expected_cluster["root"]]
                            self.assertEqual(
                                uf_cluster[0], 
                                convert_qubit_list_to_binary(
                                    expected_cluster["data_qubits"]
                                )
                            )
                            self.assertEqual(
                                uf_cluster[1], 
                                convert_qubit_list_to_binary(
                                    expected_cluster["syndrome_qubits"]
                                )
                            )


if __name__ == '__main__':
    unittest.main()