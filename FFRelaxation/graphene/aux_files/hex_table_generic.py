import numpy as np

def hex_t_init(alat, n_max, m_max, ang_min, ang_max):
  '''
  alat : The unit-cell lattice constant of 
         the 2d material for moire lattice
  ang_min : minimum of angle (degree)
  ang_max : maximum of angle (degree)
  n_max, m_max : max index to search
  '''
  f = open('unsorted', 'w')
  f.write("%s \t %s \t %s \t %s \t %s\n"%('angle', 'n  m', 'p  q','angle(rad)', 'length'))
  for n in range(1, n_max):
    for m in range(n, m_max):
      n = float(n)
      m = float(m)
      angle =   np.arccos((m**2 + 4*n*m + n**2)/(2*(m**2 + n*m + n**2)))*180./np.pi   # angles in deg
      if angle > ang_min and angle < ang_max:
        f.write("%6.10f\t%d\t%d\t%d\t%d\t%6.10f\t%6.10f\n"%(angle,  n, m, -1*m, n+m, angle*np.pi/180, np.sqrt(n**2 + n*m + m**2)*float(alat)))
  f.close()

def sort_angles(n_max, m_max, ang_min, ang_max, sorted_file, cutoff):
  '''
  cutoff : Maximum moire superlattice length 
           you wish to simulate
  sorted_file : final hex.table
  '''
  alat = input("Enter the lattice constant: ")
  # unsorted_table 
  hex_t_init(alat, n_max, m_max, ang_min, ang_max)
  angle, n, m, p, q, angle_rad, length = np.loadtxt('unsorted', unpack=True, skiprows=1)
  f = open(sorted_file, 'w+')
  f.write("%s \t %s \t %s \t %s \t %s\n\n"%('angle', 'n  m', 'p  q','angle(rad)', 'length'))

  unq_angle, indx = np.unique(angle, return_index=True)
  angle_f = [angle[indx[i]] for i in range(len(indx))]
  n_f = [n[indx[i]] for i in range(len(indx))]
  m_f = [m[indx[i]] for i in range(len(indx))]
  p_f = [p[indx[i]] for i in range(len(indx))]
  q_f = [q[indx[i]] for i in range(len(indx))]
  angle_rad_f = [angle_rad[indx[i]] for i in range(len(indx))]
  length_f = [length[indx[i]] for i in range(len(indx))]

  # save the data
  counter = 1
  for i in range(len(indx)):
    if length_f[i] <= cutoff:
      f.write("%d %6.10f %d %d %d %d %6.10f %6.10f\n"%(counter, angle_f[i], n_f[i], m_f[i], p_f[i], q_f[i], angle_rad_f[i], length_f[i]))
      counter = counter + 1
   		
  f.close()	
 
#-----------run---------
# Change these values if you want larger/smaller systems
sort_angles(100, 100, 0, 30, 'hex.table', 200)
