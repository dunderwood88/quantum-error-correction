from abc import ABC, abstractmethod
from typing import List, Union

from src.classical.helpers import convert_qubit_list_to_binary


class AbstractSurfaceCode(ABC):

    def __init__(self, dimension: int) -> None:

        if dimension & 0 or dimension == 1:
            raise ValueError("Dimension must be an odd number greater than 1!")

        self._dimension = dimension

        self._stabilizers = {
            "x": [],
            "z": []
        }

    def get_stabilizers(self, index: int = None, stabilizer_type: str = "z") -> int:
        """Returns the stabilizer of specified type in terms of the
        adjacent data qubit indexes; the result is a binary string where the
        data qubit indexes are marked 1. If index not provided, provides list
        of all the stabilizers of the given type.

        Parameters
        ----------
        index : int, by default None
            stabilizer index as defined by the code layout
        stabilizer_type : str, optional
            the type of stabilizer required, by default "z"

        Returns
        -------
        int
            binary string specifying the data qubits of the stabilizer

        Raises
        ------
        ValueError
            raise for incorrectly given stabilizer_type
        """

        if stabilizer_type not in self._stabilizers.keys():
            raise ValueError("Stabilizer type must be either 'x' or 'z'!")

        if index is None:
            return self._stabilizers[stabilizer_type]
        return self._stabilizers[stabilizer_type][index]

    def generate_syndrome(
        self, error_string: Union[int, List[int]],
        error_type: str = "x",
        show_all_adjacent: bool = False
    ) -> int:
        """Given a data qubit error string in binary or listed qubit index
        form, and the error type, generates a syndrome string in binary form
        for the corresponding syndrome qubits.

        error_type = "x" -> z-type syndrome string
        error_type = "z" -> x-type syndrome string

        Parameters
        ----------
        error_string : Union[int, List[int]]
            data qubit error string, given in binary form or as a list of
            errored data qubit indexes
        error_type : str, optional
            type of error the error string represents, by default "x"
        show_all_adjacent : bool, optional
            flag to specify whether to override the modulo 2 addition of
            stabilizers so that adjacent errors do not cancel (useful for
            situations where knowing all adjacent stabilizer qubits are, e.g.
            for growing clusters), by default False

        Returns
        -------
        int
            Corresponding syndrome string resulting from the error string
            provided
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
                if show_all_adjacent:
                    res = 1
                else:
                    res ^= par & 1
                par >>= 1

            syndrome += res << i

        return syndrome

    @abstractmethod
    def draw(
        self,
        x_data_string: Union[int, List[int]] = 0,
        z_data_string: Union[int, List[int]] = 0,
        x_syndrome_string: Union[int, List[int]] = 0,
        z_syndrome_string: Union[int, List[int]] = 0,
        restrict_graph: str = None,
        **kwargs
    ) -> None:
        """Prints a visual representation of the code to the console, with the
        ability to mark qubits by providing binary strings for each qubit type.

        Parameters
        ----------
        x_data_string : Union[int, List[int]], optional
            binary representation of x data qubits to mark as x, by default 0
        z_data_string : Union[int, List[int]], optional
            binary representation of data qubits to mark as z, by default 0
        x_syndrome_string : Union[int, List[int]], optional
            binary representation of x syndrome qubits to mark, by default 0
        z_syndrome_string : Union[int, List[int]], optional
            binary representation of z syndrome qubits to mark, by default 0
        restrict_graph : str, optional
            restrict the diagram to x or z stabilizers, if required,
            by default None
        """
        pass
