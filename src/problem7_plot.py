import numpy as np
import matplotlib.pyplot as plt

n = 100

# v1, v2, v3 = np.loadtxt(
#     "../out/problem7_n" + str(n) + ".txt",
#     usecols=(0, 1, 2),
#     unpack=True,
# )
#
# x = np.loadtxt(
#     "../out/problem7_x_n" + str(n) + ".txt"
# )

x, v1, v2, v3 = np.loadtxt(
    "../out/problem7b_n" + str(n) + ".txt",
    usecols=(0, 1, 2, 3),
    unpack=True,
)

# Figure size (inches)
# figwidth = 5.5
# figheight = figwidth / 1.33333
#
# plt.figure(figsize=(figwidth, figheight))
#
# plt.tight_layout()

plt.rc("axes", titlesize=12)  # fontsize of the axes title
plt.rc("axes", labelsize=12)  # fontsize of the x and y labels
plt.rc("xtick", labelsize=12)  # fontsize of the tick labels
plt.rc("ytick", labelsize=12)  # fontsize of the tick labels

plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
plt.title("Three solutions v($\hat{x}$) with n=" + str(n) + " steps.", wrap=True)
plt.xlabel("$\hat{x}$")
plt.ylabel("v($\hat{x}$)")


plt.plot(x, v1, linewidth=1.5)
#plt.plot(x, v1, ".", c="black", markersize=10)

plt.plot(x, v2, linewidth=1.5)
#plt.plot(x, v3, ".", c="black", markersize=10)

plt.plot(x, v3, linewidth=1.5)
#plt.plot(x, v3, ".", c="black", markersize=10)

plt.savefig("../out/problem7b_n" + str(n) + ".pdf")


plt.show()
