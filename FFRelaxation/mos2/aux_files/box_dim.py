import numpy as np


def box_dim(index):
  i = index
  latcon = [3.1, 3.11, 3.12, 3.13, 3.14]
  for j in range(len(latcon)):
    alat = latcon[j] 
    angle =  np.arccos((3*i**2 + 3*i + 0.5)/(3*i**2 + 3*i + 1))*180/np.pi
    a_moire = np.sqrt(3*i**2 + 3*i + 1)*float(alat)
    print("a_moire: %f"%(a_moire))
    print("Following are for \sqrt3\sqrt3 with alat:%f"%(alat))
    print("change_box all xy final %f x final 0.0 %f y final 0.0 %f remap units box"%(np.sqrt(3.)*a_moire/2, np.sqrt(3.)*a_moire, np.sqrt(3)*np.sqrt(3)*a_moire/2))

# run--
box_dim(21)
