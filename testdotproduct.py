import numpy as np
import numpy as np
def dot(vector, coords, data):
    result = np.zeros(len(vector))
    i = 0
    for coord in coords:
        [row, col] = coord
        result[row] += data[i] * vector[col]
        i += 1
    return result
vector = np.array([1, 2, 3])
coords = [[0, 0], [1, 0], [1, 2]]
data = [2, 3, 4]

result = dot(vector, coords, data)
print(result)
print(type(result))
