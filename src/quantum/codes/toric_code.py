from typing import List, Union

from src.classical.helpers import convert_qubit_list_to_binary
from src.quantum.codes.abstract_surface_code import AbstractSurfaceCode


class ToricCode(AbstractSurfaceCode):
    """Rotated Planar surface code, defined as coordinates over:
        D: data qubits
        X: X-type stabilizer qubits
        Z: Z-type stabilizer qubits
    Indices run from left-to-right, top-to-bottom, for all qubit types.

    Example for dimension = 3

    Z stabilizers:
        D0      D1      D2
    D3  Z0  D4  Z1  D5  Z2  D3
        D6      D7      D8
    D9  Z3  D10 Z4  D11 Z5  D9
        D12     D13     D14
    D15 Z6  D16 Z7  D17 Z8  D15
        D0      D1      D2

    X stabilizers:
    X0  D0  X1  D1  X2  D2  X0
    D3      D4      D5      D3
    X3  D6  X4  D7  X5  D8  X3
    D9      D10     D11     D9
    X6  D12 X7  D13 X8  D14 X6
    D15     D16     D17     D15
    X0  D0  X1  D1  X2  D2  X0

    combined X and Z stabilizers:
    X0      D0      X1      D1      X2      D2      X0
    D3      Z0      D4      Z1      D5      Z2      D3
    X3      D6      X4      D7      X5      D8      X3
    D9      Z3      D10     Z4      D11     Z5      D9
    X6      D12     X7      D13     X8      D14     X6
    D15     Z6      D16     Z7      D17     Z8      D15
    X0      D0      X1      D1      X2      D2      X0


    Example for dimension = 5

    Z stabilizers:
        D0      D1      D2      D3      D4
    D5  Z0  D6  Z1  D7  Z2  D8  Z3  D9  Z4  D5
        D10     D11     D12     D13     D14
    D15 Z5  D16 Z6  D17 Z7  D18 Z8  D19 Z9  D15
        D20     D21     D22     D23     D24
    D25 Z10 D26 Z11 D27 Z12 D28 Z13 D29 Z14 D25
        D30     D31     D32     D33     D34
    D35 Z15 D36 Z16 D37 Z17 D38 Z18 D39 Z19 D35
        D40     D41     D42     D43     D44
    D45 Z20 D46 Z21 D47 Z22 D48 Z23 D49 Z24 D45
        D0      D1      D2      D3      D4

    X stabilizers:
    X0  D0  X1  D1  X2  D2  X3  D3  X4  D4  X0
    D5      D6      D7      D8      D9      D5
    X5  D10 X6  D11 X7  D12 X8  D13 X9  D14 X5
    D15     D16     D17     D18     D19     D15
    X10 D20 X11 D21 X12 D22 X13 D23 X14 D24 X10
    D25     D26     D27     D28     D29     D25
    X15 D30 X16 D31 X17 D32 X18 D33 X19 D34 X15
    D35     D36     D37     D38     D39     D35
    X20 D40 X21 D41 X22 D42 X23 D43 X24 D44 X20
    D45     D46     D47     D48     D49     D45
    X0  D0  X1  D1  X2  D2  X3  D3  X4  D4  X0
    """

    def __init__(self, dimension: int) -> None:
        super().__init__(dimension)

        self._num_stabilizer_qubits = dimension**2
        self._name = f"D = {dimension} Toric Code"

        # initial Z-type stabilizers
        p_check = (1 << (2 * self._dimension)) + (3 << self._dimension) + 1

        row = 0
        for p in range(self._num_stabilizer_qubits):
            p_save = p_check << (row * self._dimension)

            if (p + 1) % self._dimension == 0:
                p_save ^= ((1 << (2 * self._dimension)) +
                           (1 << self._dimension)) << (row * 2 * self._dimension)

            if row == self._dimension - 1:
                mask = (1 << (2 * self._dimension ** 2)) - 1
                p_save &= mask
                p_save += 1 << (p % self._dimension)

            if (p + 1) % self._dimension == 0:
                row += 1

            self._stabilizers["z"].append(p_save)
            p_check <<= 1

    def draw(
        self,
        x_data_string: Union[int, List[int]] = 0,
        z_data_string: Union[int, List[int]] = 0,
        x_syndrome_string: Union[int, List[int]] = 0,
        z_syndrome_string: Union[int, List[int]] = 0,
        restrict_graph: str = None,
        **kwargs
    ) -> None:

        if "simplify" in kwargs:
            z_syndrome_label = "\033[1m" + "V"
            data_label = "e"
        else:
            z_syndrome_label = "Z"
            data_label = "D"

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
        str_code = ""
        final_row = ""
        str_code_d = ""

        row = 0
        for d in range(2 * self._dimension**2):

            has_x_error = (1 << d) & x_data_string
            has_z_error = (1 << d) & z_data_string

            if has_x_error and has_z_error:
                str_code_d = "\033[95m"
            elif has_x_error:
                str_code_d = "\033[91m"
            elif has_z_error:
                str_code_d = "\033[94m"
            else:
                str_code_d = ""

            if row % 2 == 0:
                if not restrict_graph == "z":
                    str_code += "X" + "{:<7}".format(x)
                else:
                    str_code += "{:<8}".format("")
                str_code += str_code_d + data_label + \
                    "{:<7}".format(d) + "\033[0m"

                x += 1
                if (d + 1) % self._dimension == 0:
                    if not restrict_graph == "z":
                        str_code += "X" + "{:<7}".format(x - self._dimension)
                    str_code += "\n\n\n"
                    if row == 0:
                        final_row = str_code

                    row += 1

            else:
                str_code += str_code_d + data_label + "{:<7}".format(d)
                str_code += "\033[0m"
                if not restrict_graph == "x":

                    if (1 << z) & z_syndrome:
                        str_code += "\033[93m"

                    str_code += z_syndrome_label + \
                        "{:<7}".format(z) + "\033[0m"
                    z += 1

                if (d + 1) % self._dimension == 0:

                    has_x_error = (
                        1 << (d - self._dimension + 1)) & x_data_string
                    has_z_error = (
                        1 << (d - self._dimension + 1)) & z_data_string

                    if has_x_error and has_z_error:
                        str_code_d = "\033[95m"
                    elif has_x_error:
                        str_code_d = "\033[91m"
                    elif has_z_error:
                        str_code_d = "\033[94m"
                    else:
                        str_code_d = ""

                    str_code += str_code_d + data_label + \
                        "{:<7}".format(d - self._dimension + 1) + \
                        "\033[0m" + "\n\n\n"
                    row += 1

        str_code += final_row

        print()
        print(str_code)
        print()
        print(self._name)
        print()
        if "simplify" not in kwargs:
            print("\033[91mX errors")
            print("\033[94mZ errors")
            print("\033[95mXZ errors")
            print("\033[92mX syndrome")
            print("\033[93mZ syndrome\033[0m")
            print()
