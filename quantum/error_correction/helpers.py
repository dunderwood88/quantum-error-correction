from typing import List

def convert_qubit_list_to_binary(error_string: List[int]) -> int:
    """Takes a list of marked qubit inidices and converts to a binary
    representation where unmarked qubit indices are set to 0, marked set
    to 1.

    Parameters
    ----------
    error_string : List[int]
        List of marked qubit indices

    Returns
    -------
    int
        Binary representation of unmarked and marked qubit indices
    """
    num = 0
    for i in error_string:
        num += (1 << i)
    return num