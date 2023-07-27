import math
from typing import List, Union

from src.classical.helpers import convert_qubit_list_to_binary
from src.quantum.codes.abstract_surface_code import AbstractSurfaceCode


class RPlanarCode(AbstractSurfaceCode):
    """Rotated Planar surface code, defined as coordinates over:
        D: data qubits
        X: X-type parity check qubits
        Z: Z-type parity check qubits
    Indices run from left-to-right, top-to-bottom, for all qubit types.

    Example for dimension = 3
    (a.k.a. Surface-17 https://arxiv.org/pdf/1612.08208.pdf)

    combined X and Z parity checks:
                    X0
        D0      D1      D2
    Z0      X1      Z1
        D3      D4      D5
            Z2      X2      Z3
        D6      D7      D8
            X3

    Example for dimension = 5

    combined X and Z parity checks:
                    X0              X1
        D0      D1      D2      D3      D4
    Z0      X2      Z1      X3      Z2
        D5      D6      D7      D8      D9
            Z3      X4      Z4      X5      Z5
        D10     D11     D12     D13     D14
    Z6      X6      Z7      X7      Z8
        D15     D16     D17     D18     D19
            Z9      X8      Z10     X9      Z11
        D20     D21     D22     D23     D24
            X10             X11

    """

    def __init__(self, dimension: int) -> None:
        super().__init__(dimension)

        self._num_parity_check_qubits = int((dimension**2 - 1) / 2)
        self._name = f"D = {dimension} Rotated Planar Surface Code"

        # initial Z-type parity checks
        p_check_weight_2 = (1 << self._dimension) + 1
        p_check_weight_4 = ((3 << self._dimension) + 3) << 1

        p_temp = 0
        for p in range(self._num_parity_check_qubits):

            # weight-2 Z-type parity checks
            if p_temp % self._dimension == 0:
                self._parity_checks["z"].append(p_check_weight_2)
                if p_temp == 0:
                    p_check_weight_2 = \
                        p_check_weight_2 << (2 * self._dimension) - 1
                    p_temp += 1
                elif p_temp == self._dimension:
                    p_check_weight_2 = p_check_weight_2 << 1
                    p_check_weight_4 = p_check_weight_4 << 2
                    p_temp = 0
            else:  # weight-4 Z-type parity checks
                self._parity_checks["z"].append(p_check_weight_4)
                p_check_weight_4 = p_check_weight_4 << 2
                p_temp += 1

        # initial X-type parity checks
        p_check_weight_2 = 6
        p_check_weight_4 = (3 << self._dimension) + 3

        p_temp = 0
        for p in range(self._num_parity_check_qubits):

            # weight-4 X-type parity checks
            if p >= (self._dimension - 1) / 2 and \
                    p < (self._dimension * (self._dimension - 1)) / 2:
                self._parity_checks["x"].append(p_check_weight_4)
                p_check_weight_4 = p_check_weight_4 << 2
                p_temp += 1
                if p_temp == (self._dimension - 1) / 2:
                    p_check_weight_4 = p_check_weight_4 << 2
                elif p_temp == self._dimension - 1:
                    p_temp = 0
            else:  # weight-2 X-type parity checks
                self._parity_checks["x"].append(p_check_weight_2)
                if p == ((self._dimension - 1) / 2) - 1:
                    p_check_weight_2 = \
                        p_check_weight_2 << ((self._dimension - 1)**2 + 1)
                else:
                    p_check_weight_2 = p_check_weight_2 << 2

    def draw(
        self,
        x_data_string: Union[int, List[int]] = 0,
        z_data_string: Union[int, List[int]] = 0,
        x_syndrome_string: Union[int, List[int]] = 0,
        z_syndrome_string: Union[int, List[int]] = 0,
        restrict_graph: str = None,
        **kwargs
    ) -> None:

        x_syndrome = 0
        z_syndrome = 0
        if x_data_string:
            if isinstance(x_data_string, List):
                x_data_string = convert_qubit_list_to_binary(x_data_string)
            z_syndrome = self.generate_syndrome(x_data_string, error_type="x")
        if z_data_string:
            if isinstance(z_data_string, List):
                z_data_string = convert_qubit_list_to_binary(z_data_string)
            x_syndrome = self.generate_syndrome(z_data_string, error_type="z")

        if z_syndrome_string:
            if isinstance(z_syndrome_string, List):
                z_syndrome_string = convert_qubit_list_to_binary(
                    z_syndrome_string
                )
            z_syndrome = z_syndrome_string

        x = 0
        z = 0
        str_code = "{:>16}".format("")

        # first X qubits
        for i in range(math.floor(self._dimension / 2)):
            if (1 << x) & x_syndrome:
                str_code += "\033[92m"
            if not restrict_graph == "z":
                str_code += "X" + "{:<15}".format(i)
            else:
                str_code += "{:<16}".format("")
            str_code += "\033[0m"
            x += 1
        str_code += "\n" + "{:>4}".format("")

        # D, Z and remaining X qubits
        even_row = True  # toggle flag between rows
        for d in range(self._dimension**2):
            if d > 0 and d % self._dimension == 0:
                str_code += "\n"
                if z < self._num_parity_check_qubits:
                    if not even_row:
                        str_code += "{:>8}".format("")

                    for i in range(math.ceil(self._dimension / 2)):

                        if (1 << z) & z_syndrome:
                            str_code += "\033[93m"
                        if not restrict_graph == "x":
                            str_code += "Z" + "{:<7}".format(z)
                        else:
                            str_code += "{:<8}".format("")
                        str_code += "\033[0m"
                        z += 1

                        if i != math.ceil(self._dimension / 2) - 1:

                            if (1 << x) & x_syndrome:
                                str_code += "\033[92m"
                            if not restrict_graph == "z":
                                str_code += "X" + "{:<7}".format(x)
                            else:
                                str_code += "{:<8}".format("")
                            str_code += "\033[0m"
                            x += 1

                    str_code += "\n" + "{:>4}".format("")
                    even_row = not even_row

            has_x_error = (1 << d) & x_data_string
            has_z_error = (1 << d) & z_data_string

            if has_x_error and has_z_error:
                str_code += "\033[95m"
            elif has_x_error:
                str_code += "\033[91m"
            elif has_z_error:
                str_code += "\033[94m"

            str_code += "D" + "{:<7}".format(d) + "\033[0m"
        str_code += "\n" + "{:>8}".format("")

        # last X qubits
        for i in range(x, self._num_parity_check_qubits):
            if (1 << i) & x_syndrome:
                str_code += "\033[92m"
            if not restrict_graph == "z":
                str_code += "X" + "{:<15}".format(i)
            else:
                str_code += "{:<16}".format("")
            str_code += "\033[0m"
        str_code += "\n"

        print()
        print(str_code)
        print()
        print(self._name)
        print()
        print("\033[91mX errors")
        print("\033[94mZ errors")
        print("\033[95mXZ errors")
        print("\033[92mX syndrome")
        print("\033[93mZ syndrome\033[0m")
        print()
