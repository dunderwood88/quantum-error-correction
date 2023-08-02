from src.classical.helpers import convert_qubit_list_to_binary
from src.quantum.codes.toric_code import ToricCode

from src.classical.decoders.union_find.uf_functions import (
    generate_spanning_trees, syndrome_validation_naive, tree_peeler)


### --- D = 5x7  --- ###
code = ToricCode(5, 7)
error, error_type = ([3, 13, 43], "x")
### -------------- ###

### --- D = 3x4  --- ###
code = ToricCode(3, 4)
error, error_type = ([4, 6], "x")
### -------------- ###

# extract the syndrome
error_bin = convert_qubit_list_to_binary(error)
syndrome = code.generate_syndrome(error_bin, error_type)
syndrome_type = next(s for s in ["x", "z"] if s != error_type)

# perform union-find decoding
clusters, count = syndrome_validation_naive(syndrome, code, syndrome_type)
spanning_trees = generate_spanning_trees(
    clusters, code, syndrome, syndrome_type
)
corrections = tree_peeler(spanning_trees, syndrome)
total_correction = sum(corrections.values(), [])
overall_result = error_bin ^ convert_qubit_list_to_binary(total_correction)

print(code.get_parity_checks(parity_check_type="x"))

# print the code with error + syndrome
print("BEFORE UNION-FIND DECODING")
code.draw(
    x_data_string=error_bin,
    restrict_graph="z"
    # z_data_string=error_bin,
    # restrict_graph="x"
)

# print the code with applied correction + syndrome
print("AFTER UNION-FIND DECODING")
code.draw(
    x_data_string=overall_result,
    z_syndrome_string=syndrome,
    restrict_graph="z"
    # z_data_string=overall_result,
    # x_syndrome_string=syndrome,
    # restrict_graph="x"
)

print()
print(f"UNION-FIND DECODING COMPLETED IN {count} ROUNDS")
print(f"SUGGESTED DATA QUBIT CORRECTIONS: {total_correction}")
print()
