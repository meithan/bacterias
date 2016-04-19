import matplotlib
import matplotlib.pyplot as plt
import time
import sys

if len(sys.argv) != 2:
  print("Must provide simulation data filename as argument")
  exit()
else:
  datafname = sys.argv[1]

cycles = []
total_population = []
num_food_sources = []
strain_populations = []

f = open(datafname)
for line in f:
  data = map(int, line.split())
  if len(strain_populations) == 0:
    num_strains = len(data) - 3
    strain_populations = [[] for i in range(num_strains)]
  cycles.append(data[0])
  total_population.append(data[1])
  num_food_sources.append(data[2])
  for i in range(num_strains):
    strain_populations[i].append(data[3+i])
f.close()

plt.ion()

plt.figure()
plt.plot(cycles, total_population, "-", color="blue")
plt.xlabel("Cycles")
plt.ylabel("Total population")
plt.grid()
ax2 = plt.gca().twinx()
ax2.plot(cycles, num_food_sources, "-", color="orange")
ax2.set_ylabel("Number of food sources")
#plt.savefig("test1.png")

plt.figure()
for i in range(num_strains):
  plt.plot(cycles, strain_populations[i], "-", lw=1.0)
plt.xlabel("Cycles")
plt.ylabel("Strain populations")
plt.grid()
#plt.savefig("test2.png")

raw_input("Press ENTER to exit")
