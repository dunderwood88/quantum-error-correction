from typing import List, Union

import numpy as np
from termcolor import colored

from src.classical.helpers import convert_qubit_list_to_binary
from src.quantum.codes.abstract_surface_code import AbstractSurfaceCode


class ToricCode(AbstractSurfaceCode):
    """Rotated Planar surface code, defined as coordinates over:
        D: data qubits
        X: X-type parity check qubits
        Z: Z-type parity check qubits
    Indices run from left-to-right, top-to-bottom, for all qubit types.

    Example for dimension = 3x3

    Z parity checks:
        D0      D1      D2
    D3  Z0  D4  Z1  D5  Z2  D3
        D6      D7      D8
    D9  Z3  D10 Z4  D11 Z5  D9
        D12     D13     D14
    D15 Z6  D16 Z7  D17 Z8  D15
        D0      D1      D2

    X parity checks:
    X0  D0  X1  D1  X2  D2  X0
    D3      D4      D5      D3
    X3  D6  X4  D7  X5  D8  X3
    D9      D10     D11     D9
    X6  D12 X7  D13 X8  D14 X6
    D15     D16     D17     D15
    X0  D0  X1  D1  X2  D2  X0

    combined X and Z parity checks:
    X0      D0      X1      D1      X2      D2      X0
    D3      Z0      D4      Z1      D5      Z2      D3
    X3      D6      X4      D7      X5      D8      X3
    D9      Z3      D10     Z4      D11     Z5      D9
    X6      D12     X7      D13     X8      D14     X6
    D15     Z6      D16     Z7      D17     Z8      D15
    X0      D0      X1      D1      X2      D2      X0


    Example for dimension = 5x5

    Z parity checks:
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

    X parity checks:
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

    def __init__(self, width: int, length: int) -> None:
        super().__init__(width, length)

        self._name = f"{self._width}x{self._length} Toric Code"

        # initial Z-type parity check
        p_check = (1 << (2 * self._width)) + (3 << self._width) + 1

        for p in range(self._num_parity_check_qubits):
            row = p // self._width
            p_save = p_check << (row * self._width)

            if (p + 1) % self._width == 0:
                p_save ^= ((1 << (2 * self._width)) +
                           (1 << self._width)) << (row * 2 * self._width)

            if row == self._length - 1:
                mask = (1 << (2 * self._num_parity_check_qubits)) - 1
                p_save &= mask
                p_save += 1 << (p % self._width)

            self._parity_checks["z"].append(p_save)
            p_check <<= 1

        # final X-type parity check (generate in reverse)
        final_index = (2 * self._num_parity_check_qubits) - 1
        p_check = (1 << final_index) + \
                  ((1 << final_index) >> (2 * self._width)) + \
                  ((3 << final_index - 1) >> self._width) \

        for p in range(self._num_parity_check_qubits):
            row = p // self._width
            p_save = p_check >> (row * self._width)

            if (p + 1) % self._width == 0:
                p_save ^= (1 << (final_index - self._width)) + \
                          ((1 << final_index) >> (2 * self._width))

                if p != ((self._num_parity_check_qubits) - 1):
                    final_index -= (2 * self._width)

            if row == self._length - 1:
                p_save |= (
                    1 << ((2 * self._num_parity_check_qubits) -
                          1 - (p % self._width))
                    )

            self._parity_checks["x"].insert(0, p_save)
            p_check >>= 1

        self._parity_check_matrices = {
            "x": None,
            "z": None
        }

    def get_parity_check_matrices(self, check_type: str = "x") -> np.ndarray:
        """Returns either the x- or z-type parity check matrices for the code
        in the form of a numpy array.

        Each row represents a parity check qubit:
            - x-type = measure-X qubit / Z-stabilizer / face
            - z-type = measure-Z qubit / X-stabilizer / vertex
        Each column represents a data qubit.

        Parameters
        ----------
        check_type : str, optional
            the type of parity check matrix to return (x or z), by default "x"

        Returns
        -------
        np.ndarray
            the parity check matrix for the given parity check type (x or z)
        """

        if self._parity_check_matrices[check_type]:
            return self._parity_check_matrices[check_type]
        else:
            self._parity_check_matrices[check_type] = np.zeros(
                (self._num_parity_check_qubits, self._num_data_qubits),
                dtype=np.uint8
            )
            for row, parity_check in enumerate(self._parity_checks[check_type]):
                p = parity_check
                for column in range(self._num_data_qubits):
                    if p & 1:
                        self._parity_check_matrices[check_type][row, column] = 1
                    p >>= 1

        return self._parity_check_matrices[check_type]

    def _set_num_parity_check_qubits(self):
        self._num_parity_check_qubits = self._width * self._length

    def _set_num_data_qubits(self):
        self._num_data_qubits = 2 * (self._width * self._length)

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

        if x_syndrome_string:
            if isinstance(x_syndrome_string, List):
                x_syndrome_string = convert_qubit_list_to_binary(
                    x_syndrome_string
                )
            x_syndrome = x_syndrome_string

        # set up styling
        if "simplify" in kwargs:
            z_syndrome_label = "V"
            data_label = "E"
        else:
            z_syndrome_label = "Z"
            x_syndrome_label = "X"
            data_label = "D"

        str_code = colored("{:<7}".format(""), on_color="on_black")
        str_padding = str_code + colored(
            ((2 * self._width) + 1) * "{:<8}".format(""), on_color="on_black"
        ) + "\n"
        str_final_row = ""

        # draw the graph
        x = 0
        z = 0
        row = 0
        for d in range(2 * self._num_parity_check_qubits):

            has_x_error = (1 << d) & x_data_string
            has_z_error = (1 << d) & z_data_string

            if has_x_error and has_z_error:
                d_color = "magenta"
            elif has_x_error:
                d_color = "red"
            elif has_z_error:
                d_color = "blue"
            else:
                d_color = "light_grey"

            if row % 2 == 0:  # X-syndromes only (even rows)

                # 1) X-syndrome qubit
                x_color = "dark_grey"
                x_style = []
                if not restrict_graph == "z":  # show X-syndromes if wanted
                    if (1 << x) & x_syndrome:
                        x_color = "green"
                        x_style.append("bold")

                    str_qubit = colored(
                        x_syndrome_label + "{:<7}".format(x),
                        color=x_color,
                        attrs=x_style,
                        on_color="on_black"
                    )

                else:
                    str_qubit = colored(
                        "{:<8}".format(""), on_color="on_black"
                    )

                str_code += str_qubit
                if d == (row * self._width):
                    str_final_qubit = str_qubit

                # 2) data qubit
                str_code += colored(
                    data_label + "{:<7}".format(d),
                    color=d_color,
                    on_color="on_black"
                )

                # 3) deal with final column
                if (d + 1) % self._width == 0:  # if final column

                    str_code += str_final_qubit

                    if row == 0:
                        str_final_row = str_code  # final row == first row
                    row += 1

                    str_code += "\n" + 2 * str_padding
                    str_code += colored(
                        "{:<7}".format(""), on_color="on_black"
                    )

                x += 1

            else:  # Z-syndromes only (odd rows)

                # 1) data qubit
                str_qubit = colored(
                    data_label + "{:<7}".format(d),
                    color=d_color,
                    on_color="on_black"
                )

                str_code += str_qubit

                if d == (row * self._width):
                    str_final_qubit = str_qubit

                # 2) Z-syndrome qubit
                z_color = "dark_grey"
                z_style = []
                if not restrict_graph == "x":  # show Z-syndromes if wanted
                    if (1 << z) & z_syndrome:
                        z_color = "yellow"
                        z_style.append("bold")

                    str_code += colored(
                        z_syndrome_label + "{:<7}".format(z),
                        color=z_color,
                        attrs=z_style,
                        on_color="on_black"
                    )

                else:
                    str_code += colored(
                        "{:<8}".format(""), on_color="on_black"
                    )

                # 3) deal with final column
                if (d + 1) % self._width == 0:  # if final column

                    str_code += str_final_qubit
                    str_code += "\n" + 2 * str_padding

                    if row != ((2 * self._length) - 1):
                        str_code += colored(
                            "{:<7}".format(""), on_color="on_black"
                        )
                    row += 1

                z += 1

        str_code += str_final_row

        print()
        print(colored(self._name, attrs=["bold"]))
        print()
        print(2 * str_padding + str_code + "\n" + 2 * str_padding)

        if "simplify" not in kwargs:
            print(colored("X errors", "red", attrs=["bold"]) + " - " + colored("Z syndrome", "yellow", attrs=["bold"]))
            print(colored("Z errors", "blue", attrs=["bold"]) + " - " + colored("X syndrome", "green", attrs=["bold"]))
            print(colored("XZ errors", "magenta", attrs=["bold"]))
            print()
