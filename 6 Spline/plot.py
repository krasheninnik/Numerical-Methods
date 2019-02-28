import matplotlib.pyplot as plt
import numpy as np



x_data = []
y_data = []

with open("spline.txt") as data:
    for line in data:
        line = line.lstrip()
        line = line.split(" ")
        x_data.append(float(line[0]))
        y_data.append(float(line[len(line)-1]))

plt.plot(x_data, y_data)

plt.ylabel("f(x)")
plt.xlabel("x")

x_min = min( x_data)
x_max = max(x_data)
y_min = min(y_data)
y_max = max(y_data)

_max = max(x_max, y_max)
_min = min(x_min, y_min)

plt.axis([_min, _max, _min, _max])

plt.grid()
plt.show()