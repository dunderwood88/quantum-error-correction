from abc import ABC, abstractmethod
from typing import List, Union

from src.classical.helpers import convert_qubit_list_to_binary


class AbstractSurfaceCode(ABC):

    def __init__(self, width: int, length: int) -> None:

        if width <= 1 or length <= 1:
            raise ValueError("Each dimension must be greater than 1!")

        self._width = width
        self._length = length
        self._set_num_parity_check_qubits()
        self._set_num_data_qubits()

        self._parity_checks = {
            "x": [],
            "z": []
        }

    def get_parity_checks(
        self,
        index: int = None,
        parity_check_type: str = "z"
    ) -> int:
        """Returns the parity check of specified type (x or z) in terms of the
        adjacent data qubit indexes; the result is a binary string where the
        data qubit indexes are marked 1. If index not provided, provides list
        of all the parity checks of the given type.

        Parameters
        ----------
        index : int, by default None
            parity check index as defined by the code layout
        parity_check_type : str, optional
            the type of parity check required ("x" or "z"), by default "z"

        Returns
        -------
        int
            binary string specifying data qubits involved in the parity check

        Raises
        ------
        ValueError
            raise for incorrectly given parity_check_type
        """

        if parity_check_type not in self._parity_checks.keys():
            raise ValueError("Parity check type must be either 'x' or 'z'!")

        if index is None:
            return self._parity_checks[parity_check_type]
        return self._parity_checks[parity_check_type][index]

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
            parity checks so that adjacent errors do not cancel (useful for
            situations where knowing all adjacent parity check qubits are, e.g.
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
        for i, parity_check in enumerate(self._parity_checks[syndrome_type]):
            res = 0
            par = parity_check & error_string
            while par:
                if show_all_adjacent:
                    res = 1
                else:
                    res ^= par & 1
                par >>= 1

            syndrome += res << i

        return syndrome

    @abstractmethod
    def _set_num_parity_check_qubits(self):
        """Method to initialise the number of parity check qubits for the given
        code.
        """
        pass

    @property
    def num_parity_check_qubits(self):
        return self._num_parity_check_qubits

    @abstractmethod
    def _set_num_data_qubits(self):
        """Method to initialise the number of data qubits for the given code.
        """
        pass

    @property
    def num_data_qubits(self):
        return self._num_data_qubits

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
            restrict the diagram to x- or z- type parity checks, if required,
            by default None
        """
        pass
