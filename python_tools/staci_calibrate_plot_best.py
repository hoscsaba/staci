import matplotlib.pyplot as plt
import numpy as np
import sys

filename = str(sys.argv[1])
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

x = np.linspace(0, len(data), len(data))/2
data = np.array(data)

plt.figure(1)
for i in range(0, 7):
    plt.subplot(3, 3, i+1)
    plt.plot(x, data[:, 2*i], x, data[:, 2*i+1])
    plt.xlim([0, 24])
    plt.title(titles[2*i+1])
    if i > 5:
        plt.xlabel("time (hours)")
    if (i == 0) or (i == 3) or (i == 6):
        plt.ylabel("Pool water level (m)")
    if i < 6:
        plt.gca().axes.xaxis.set_ticklabels([])
    plt.grid(True)

# plt.show()
plt.savefig("pools.pdf")

plt.figure(2)
sf = 8
for i in range(0, 9):
    plt.subplot(3, 3, i+1)
    plt.plot(x, data[:, sf+2*i], x, data[:, sf+2*i+1])
    plt.xlim([0, 24])
    plt.title(titles[sf+2*i+1])
    if i > 5:
        plt.xlabel("time (hours)")
    if (i == 0) or (i == 3) or (i == 6):
        plt.ylabel("pressure (mwc)")
    if i < 6:
        plt.gca().axes.xaxis.set_ticklabels([])
    plt.grid(True)

# plt.show()
plt.savefig("pressures1.pdf")

plt.figure(3)
sf = 2*9
for i in range(0, 9):
    plt.subplot(3, 3, i+1)
    plt.plot(x, data[:, sf+2*i], x, data[:, sf+2*i+1])
    plt.xlim([0, 24])
    plt.title(titles[sf+2*i+1])
    if i > 5:
        plt.xlabel("time (hours)")
    if (i == 0) or (i == 3) or (i == 6):
        plt.ylabel("pressure (mwc)")
    if i < 6:
        plt.gca().axes.xaxis.set_ticklabels([])
    plt.grid(True)

# plt.show()
plt.savefig("pressures2.pdf")

plt.figure(4)
sf = 3*9
for i in range(0, 9):
    plt.subplot(3, 3, i+1)
    plt.plot(x, data[:, sf+2*i], x, data[:, sf+2*i+1])
    plt.xlim([0, 24])
    plt.title(titles[sf+2*i+1])
    if i > 5:
        plt.xlabel("time (hours)")
    if (i == 0) or (i == 3) or (i == 6):
        plt.ylabel("pressure (mwc)")
    if i < 6:
        plt.gca().axes.xaxis.set_ticklabels([])
    plt.grid(True)

# plt.show()
plt.savefig("pressures3.pdf")

plt.figure(5)
sf = 4*9
for i in range(0, 9):
    plt.subplot(3, 3, i+1)
    plt.plot(x, data[:, sf+2*i], x, data[:, sf+2*i+1])
    plt.xlim([0, 24])
    plt.title(titles[sf+2*i+1])
    if i > 5:
        plt.xlabel("time (hours)")
    if (i == 0) or (i == 3) or (i == 6):
        plt.ylabel("pressure (mwc)")
    if i < 6:
        plt.gca().axes.xaxis.set_ticklabels([])
    plt.grid(True)

# plt.show()
plt.savefig("pressures4.pdf")
