import numpy as np
import matplotlib.pyplot as plt

n = 10

x, v1, v2, v3 = np.loadtxt(
    "../out/problem7_n" + str(n) + ".txt",
    usecols=(0, 1, 2, 3),
    unpack=True,
    skiprows=0,
)

# Figure size (inches)
figwidth = 5.5
figheight = figwidth / 1.33333


# Plot absolute error vs h on log-log axes
plt.figure(figsize=(figwidth, figheight))

plt.rc("axes", titlesize=12)  # fontsize of the axes title
plt.rc("axes", labelsize=12)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=12)  # fontsize of the tick labels
plt.rc("ytick", labelsize=12)  # fontsize of the tick labels

plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
plt.title("Plot for discretization of x^ with n=" + str(n) + " steps.", wrap=True)
plt.xlabel("x^")
plt.ylabel("v")


plt.plot(x, v1, "--", c="0.8", linewidth=1.5)
plt.plot(x, v1, ".", c="black", markersize=10)

plt.plot(x, v2, "--", c="0.8", linewidth=1.5)
plt.plot(x, v3, ".", c="black", markersize=10)

plt.plot(x, v3, "--", c="0.8", linewidth=1.5)
plt.plot(x, v3, ".", c="black", markersize=10)

plt.savefig("../out/problem7_n" + str(n) + ".pdf")
plt.show()
