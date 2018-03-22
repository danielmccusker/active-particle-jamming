import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections

num = 5000
sizes = 0.2 * np.random.random(num)
xy = 50 * np.random.random((num, 2))

# Note that the patches won't be added to the axes, instead a collection will
patches = [plt.Circle(center, size) for center, size in zip(xy, sizes)]
print xy
print sizes
fig, ax = plt.subplots()

coll = matplotlib.collections.PatchCollection(patches, facecolors='black')
ax.add_collection(coll)

ax.margins(0.01)
plt.show()
