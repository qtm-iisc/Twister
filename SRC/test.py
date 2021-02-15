from matplotlib.path import Path
import numpy as np

tupVerts=[(86, 52), (85, 52), (81, 53), (80, 52), (79, 48), (81, 49), (86, 53),
 (85, 51), (82, 54), (84, 54), (83, 49), (81, 52), (80, 50), (81, 48),
 (85, 50), (86, 54), (85, 54), (80, 48), (79, 50), (85, 49), (80, 51),
 (85, 53), (82, 49), (83, 54), (82, 53), (84, 49), (79, 49)]


x, y = np.meshgrid(np.arange(300), np.arange(300)) # make a canvas with coordinates
x, y = x.flatten(), y.flatten()
points = np.vstack((x,y)).T 
print points

p = Path(tupVerts) # make a polygon
grid = p.contains_points(points)
print grid
mask = grid.reshape(300,300) # now you have a mask with points inside a polygon

