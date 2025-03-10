import numpy as np
import pymsis
import matplotlib.pyplot as plt

# dates = np.arange(np.datetime64("2003-10-28T00:00"), np.datetime64("2004-01-04T00:00"), np.timedelta64(1, "h"))
# # geomagnetic_activity=-1 is a storm-time run
# data = pymsis.calculate(dates, 0, 0, 100, geomagnetic_activity=-1)
# print(data[:, 0, 0, 0, 0])
# plt.plot(dates, data[:, 0,0,0,0])
# plt.show()

date = np.datetime64("2025-02-24T00:00")
# geomagnetic_activity=-1 is a storm-time run
data = pymsis.calculate(date, 0, 0, 100)

# output format:
# [Total mass density (kg/m3),
# N2 # density (m-3),
# O2 # density (m-3),
# O # density (m-3),
# He # density (m-3),
# H # density (m-3),
# Ar # density (m-3),
# N # density (m-3),
# Anomalous oxygen # density (m-3),
# NO # density (m-3),
# Temperature (K)]
print(data)

# mass of a single N_2 molecule
m_N2 = 4.6518e-26

# assume that the air only consists of N_2
num_particles = data[0,0]/m_N2
print(num_particles)
# print(sum(data[0,~np.isnan(data[0,:])][1:-1]))

alts = np.linspace(0,1000, 10000)
data = pymsis.calculate(date, 0, 0, alts)

data = np.squeeze(data)
print(data.shape)

fig,ax = plt.subplots()
ax.plot(alts, data[:, 0]/m_N2)
# ax.set_xscale("log")
ax.set_yscale("log")
plt.show()
