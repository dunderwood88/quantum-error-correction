import matplotlib.pyplot as plt
import numpy as np

from src.classical.helpers import convert_qubit_list_to_binary
from src.quantum.codes.toric_code import ToricCode

from src.classical.decoders.union_find.uf_functions import (
    generate_spanning_trees, syndrome_validation_naive, tree_peeler)


def is_logical_error(code: ToricCode, error):
    



# ### --- D = 5  --- ###
# code = ToricCode(5)
# x_error = [1, 14, 17]
# x_error = [1, 15]
# x_error = [1, 7, 12]
# x_error = [3, 13, 43]
# ### -------------- ###

dimension = 7
num_data_qubits = dimension**2 - dimension
max_cycles = 5000

results = {}

for dimension in [4, 8]:
    dim_results = {}

    # for prob in [
    #     0.5, 0.4, 0.3, 0.2, 0.1,
    #     0.05, 0.04, 0.03, 0.02, 0.01,
    #     0.005, 0.004, 0.003, 0.002, 0.001,
    #     0.0005, 0.0004, 0.0003, 0.0002, 0.0001
    # ]:
    for prob in np.linspace(0.01, 0.2, 9):

        logical_fail = 0

        cycle = 0
        while cycle < max_cycles:

            error_string = 0
            for i in range(num_data_qubits):
                flip = np.random.binomial(1, prob,)
                error_string |= flip
                error_string <<= 1

            code = ToricCode(dimension)
            error_syn = code.generate_syndrome(error_string, error_type="x")

            # ignore if no errors occurred or were detected...
            if not error_syn:
                continue

            # ...otherwise try and correct them
            clusters, count = syndrome_validation_naive(error_syn, code)
            spanning_trees = generate_spanning_trees(clusters, code, error_syn)
            corrections = tree_peeler(spanning_trees, error_syn)

            # check if the overall correction was successful
            total_correction = convert_qubit_list_to_binary(
                sum(corrections.values(), [])
            )
            result = code.generate_syndrome(
                error_string ^ total_correction,
                error_type="x"
            )

            if result or is_logical_error(code, error_string):
                logical_fail += 1

        dim_results[prob] = logical_fail / max_cycles

    results[dimension] = dim_results

print(results)


for d, res in results.items():
    plt.plot(list(res.keys()), list(res.values()), label=f"dimension: {d}")
plt.xlabel("Physical error rate")
plt.ylabel("Logical error rate")
plt.legend(loc=0)

plt.savefig("threshold.png")
