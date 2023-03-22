def syndrome_validation(dim: int, syndrome: list) -> dict, dict:
    
    odd_clusters = []

    while odd_clusters:
        
        # grow cluster by a half-edge



def unionfind_decoder(dim: int, syndrome: list) -> list:
    """Given a code dimension and a syndrome specification, performs the
    union-find algorithm to generate the most likely correction graph to apply

    Parameters
    ----------
    dim : int
        dimension of the code
    syndrome : list
        a list of syndrome vertices

    Returns
    -------
    list
        _description_
    """

    if len(syndrome) % 2 != 0:
        raise ValueError(
            "Requires an even number of syndrome vertices."
        )

    if dim < 3:
        raise ValueError(
            "Lattice dimension too small!"
        )

    # step 1: syndrome validation
    support, clusters = syndrome_validation(dim, syndrome)

    # step 2: spanning tree

    # step 3: peeling decoder

if __name__ == "__main__":
    print(unionfind_decoder(3, [1, 2, 3]))