import pandas
import pathlib
import re
# from tensorflow.keras.models import Sequential, load_module
# from tensorflow.keras.layers import Conv1D
# from sklearn.model_selection import train_test_split
# from sklearn.metrics import accuracy_score
# print(tf.__version__)
# exit()
import random
import numpy as np
import plotly.graph_objects as go
import datetime
# from matplotlib import pyplot as plt
# from contextlib import redirect_stdout


cwd = pathlib.Path.cwd()
orbit_data = pandas.read_csv(str(cwd / "Orbit_Table.csv"), header = 0).set_index("t")
# orbit_tran = orbit_data.T

sys = go.Figure()

rows, cols = orbit_data.shape

for c in range(int(cols / 2)):
    cndex = (2*c, 2*c+1)
    # print(cndex)
    col_of_interest = orbit_data.iloc[:, cndex[0]:cndex[1] + 1].copy()
    # print(col_of_interest)
    sys.add_trace(go.Scatter(x = col_of_interest.iloc[:, 0], y = col_of_interest.iloc[:, 1]))
    # break
# print(orbit_data.shape)
# print(orbit_tran)
# print(orbit_tran.index)

sys.show()
# # print(orbit_data)


# objects = list(orbit_data)

# bodies = list()
# for o in objects:
#     body = re.split("_", o)[0]
#     bodies.append(body)

# bodies = set(bodies)
