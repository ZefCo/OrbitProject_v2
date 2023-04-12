import pandas
import pathlib
import re
import random
import numpy as np
import plotly.graph_objects as go
import datetime
# from matplotlib import pyplot as plt
# from contextlib import redirect_stdout

# https://stackoverflow.com/questions/63328589/animated-plot-with-plotly


cwd = pathlib.Path.cwd()
orbit_data = pandas.read_csv(str(cwd / "Orbit_Table.csv"), header = 0).set_index("t")
# orbit_tran = orbit_data.T

sys = go.Figure(layout = go.Layout(updatemenus=[dict(type="buttons", buttons=[dict(label="Play", method="animate", args=[None])])]))

rows, cols = orbit_data.shape

for c in range(int(cols / 2)):
    cndex = (2*c, 2*c+1)
    # print(cndex)
    col_of_interest: pandas.DataFrame = orbit_data.iloc[:, cndex[0]:cndex[1] + 1].copy()
    # print(tuple(col_of_interest.columns)[0])
    name = re.split("_", tuple(col_of_interest.columns)[0])[0]
    # exit()
    xm, xM = col_of_interest.iloc[0, 0], col_of_interest.iloc[rows - 1, 0]
    ym, yM = col_of_interest.iloc[0, 1], col_of_interest.iloc[rows - 1, 1]
    sys.add_trace(go.Scatter(x = col_of_interest.iloc[:, 0], y = col_of_interest.iloc[:, 1], name = name, mode = "lines", line = dict(width = 2)))
    sys.update_layout(xaxis = dict(range = [xm, xM], autorange = False, zeroline = False),
                      yaxis = dict(range = [ym, yM], autorange = False, zeroline = False))
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
