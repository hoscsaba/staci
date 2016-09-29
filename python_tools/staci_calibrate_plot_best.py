import matplotlib.pyplot as plt
import numpy as np
import sys

plt.close('all')

filename = str(sys.argv[1])

print("")
print("Datafile plotted: %s" % filename)

pipe_names = []
pipe_dia = []
pipe_dia_mul = []
titles = []
data = []

with open(filename, 'r') as f:
    pipe_names = f.readline().split(';')
    pipe_names.pop()

    pipe_dia = f.readline().split(';')
    pipe_dia.pop()
    [float(i) for i in pipe_dia]

    pipe_dia_mul = f.readline().split(';')
    pipe_dia_mul.pop()
    [float(i) for i in pipe_dia_mul]

    titles = f.readline().split(';')
    titles.pop()
    titles = map(str.lower, titles)

    for line in f:
        number_strings = line.split(';')
        number_strings.pop()
        data.append([float(i) for i in number_strings])

print("Size of pipe_names is   %d" % len(pipe_names))
print("Size of pipe_dia is     %d" % len(pipe_dia))
print("Size of pipe_dia_mul is %d" % len(pipe_dia_mul))
print("Size of titles is       %d" % len(titles))
print("Number of rows in data: %d" % len(data))
print("")

for i in range(0, len(titles)/2):
    print(" %2d - %s" % (i, titles[2*i+1]))

plot_this = input("Choose data to plot: ")

x = np.linspace(0, len(data), len(data))/2
data = np.array(data)

plt.figure(1)
plt.plot(x, data[:, 2*plot_this], x, data[:, 2*plot_this+1])
plt.xlim([0, 24])
plt.title(titles[2*plot_this+1])
plt.xlabel("time (hours)")
plt.ylabel("Pool water level (m)")
plt.grid(True)
plt.show()
# plt.savefig("pools.pdf")
