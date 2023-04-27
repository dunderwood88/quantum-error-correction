import math
from typing import List, Tuple, Union

from quantum.error_correction.helpers import convert_qubit_list_to_binary

class RPlanarCode:
    """Rotated Planar surface code, defined as coordinates over:
        D: data qubits
        X: X-type stabilizer qubits
        Z: Z-type stabilizer qubits
    Indices run from left-to-right, top-to-bottom, for all qubit types.

    Example for dimension = 3
    (a.k.a. Surface-17 https://arxiv.org/pdf/1612.08208.pdf)
    
                    X0
        D0      D1      D2
    Z0      X1      Z1
        D3      D4      D5
            Z2      X2      Z3
        D6      D7      D8
            X3

    Example for dimension = 5
    
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

        if dimension & 0 or dimension == 1:
            raise ValueError("Dimension must be an odd number greater than 1!")

        self._dimension = dimension
        self._num_stabilizer_qubits = int((dimension**2 - 1) / 2)
        self._name = f"D = {dimension} Rotated Planar Surface Code"

        self._stabilizers = {
            "x": [],
            "z": []
        }

        # initial Z-type stabilizers
        p_check_weight_2 = (1 << self._dimension) + 1
        p_check_weight_4 = ((3 << self._dimension) + 3) << 1
        
        p_temp = 0
        for p in range(self._num_stabilizer_qubits):

            # weight-2 Z-type stabilizers
            if p_temp % self._dimension == 0:
                self._stabilizers["z"].append(p_check_weight_2)
                if p_temp == 0:
                    p_check_weight_2 = \
                        p_check_weight_2 << (2 * self._dimension) - 1
                    p_temp += 1
                elif p_temp == self._dimension:
                    p_check_weight_2 = p_check_weight_2 << 1
                    p_check_weight_4 = p_check_weight_4 << 2
                    p_temp = 0
            else:  # weight-4 Z-type stabilizers
                self._stabilizers["z"].append(p_check_weight_4)
                p_check_weight_4 = p_check_weight_4 << 2
                p_temp += 1

        # initial X-type stabilizers
        p_check_weight_2 = 6
        p_check_weight_4 = (3 << self._dimension) + 3

        p_temp = 0
        for p in range(self._num_stabilizer_qubits):

            # weight-4 X-type stabilizers
            if p >= (self._dimension - 1) / 2 and \
                p < (self._dimension * (self._dimension - 1)) / 2:
                self._stabilizers["x"].append(p_check_weight_4)
                p_check_weight_4 = p_check_weight_4 << 2
                p_temp += 1
                if p_temp == (self._dimension - 1) / 2:
                    p_check_weight_4 = p_check_weight_4 << 2
                elif p_temp == self._dimension - 1:
                    p_temp = 0
            else:  # weight-2 X-type stabilizers
                self._stabilizers["x"].append(p_check_weight_2)
                if p == ((self._dimension - 1) / 2) - 1:
                    p_check_weight_2 = \
                        p_check_weight_2 << ((self._dimension - 1)**2 + 1)
                else:
                    p_check_weight_2 = p_check_weight_2 << 2

    def _get_qubit_row_column(
        self, index: int, type: str="z"
    ) -> Tuple[int, int]:
        
        if (type == "x" or type == "z") and \
            (index > self._num_stabilizer_qubits - 1):
            raise ValueError("Index greater than number of stabilizer qubits!")


        if type == "z":
            dim = math.ceil(self._dimension / 2)
        elif type == "x":
            dim = math.floor(self._dimension / 2)
        elif type == "d":
            dim = self._dimension

        else:
            raise ValueError("Invalid qubit type!")

        return math.floor(index / dim), index % dim


    def generate_syndrome(
        self, error_string: Union[int, List[int]], error_type: str="x"
    ) -> int:
        """Given a data qubit error string in binary or listed qubit index 
        form, and the error type, generate a syndrome string in binary form for
        the corresponding syndrome qubits.
        
        error_type = "x" -> z-type syndrome string
        error_type = "z" -> x-type syndrome string

        Parameters
        ----------
        error_string : Union[int, List[int]]
            The data qubit error string, given in binary form or as a list of
            errored data qubit indices
        error_type : str, optional
            The type of error the error_string represents, by default "x"

        Returns
        -------
        int
            Corresponding syndrome string for the error_string provided
        """

        if isinstance(error_string, List):
            error_string = convert_qubit_list_to_binary(error_string)

        # syndrome type is opposite to error type
        syndrome_type = next(s for s in ["x", "z"] if s != error_type)

        # build the syndrome string
        syndrome = 0
        for i, stabilizer in enumerate(self._stabilizers[syndrome_type]):
            res = 0
            par = stabilizer & error_string
            while par:
                res ^= par & 1
                par >>= 1

            syndrome += res << i

        return syndrome 

    def draw(
        self,
        x_error_string: Union[int, List[int]]=0,
        z_error_string: Union[int, List[int]]=0,
        z_syndrome_string: Union[int, List[int]]=0,
        restrict_graph: str=None
    ) -> None:
        """Method to print diagrammatic string representation of the code to
        the terminal.
        """

        if restrict_graph not in [None, "x", "z"]:
            raise ValueError("restrict graph must be either 'x' or 'z'!")
        
        x_syndrome = 0
        z_syndrome = 0
        if x_error_string:
            if isinstance(x_error_string, List):
                x_error_string = convert_qubit_list_to_binary(x_error_string)
            z_syndrome = self.generate_syndrome(x_error_string, error_type="x")
        if z_error_string:
            if isinstance(z_error_string, List):
                z_error_string = convert_qubit_list_to_binary(z_error_string)
            x_syndrome = self.generate_syndrome(z_error_string, error_type="z")

        if z_syndrome_string:
            if isinstance(z_syndrome, List):
                z_syndrome_string = convert_qubit_list_to_binary(z_syndrome_string)
            z_syndrome |= z_syndrome_string

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
                if z < self._num_stabilizer_qubits:
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

            has_x_error = (1 << d) & x_error_string
            has_z_error = (1 << d) & z_error_string

            if has_x_error and has_z_error:
                str_code += "\033[95m"
            elif has_x_error:
                str_code += "\033[91m"
            elif has_z_error:
                str_code += "\033[94m"

            str_code += "D" + "{:<7}".format(d) + "\033[0m"
        str_code += "\n" + "{:>8}".format("")

        # last X qubits
        for i in range(x, self._num_stabilizer_qubits):
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
