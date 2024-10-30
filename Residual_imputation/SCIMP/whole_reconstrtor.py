import numpy as np
import pandas as pd

def balance_data(log_data, step_encode):
    log_matrix = log_data.values
    step_matrix = step_encode.values
    total_changed_count = 0
    changed = True

    while changed:
        changed = False
        correlation_matrix = np.corrcoef(log_matrix, rowvar=False)
        top_2_indices = np.argsort(-correlation_matrix, axis=1)[:, 1:3]

        for i in range(step_matrix.shape[1]):
            relevant_columns = top_2_indices[i]
            mask = (step_matrix[:, relevant_columns] == 0).all(axis=1)
            to_change = mask & (step_matrix[:, i] != 0)

            if np.any(to_change):
                changed = True
                changed_count = np.sum(to_change)
                total_changed_count += changed_count
                step_matrix[to_change, i] = 0

        if not changed:
            break

    print(f"Total number of elements changed: {total_changed_count}")
    modified_step_encode = pd.DataFrame(step_matrix, columns=step_encode.columns)
    return modified_step_encode