import numpy as np

def commensurate_least(index, hex_table):
  alat = input()
  f = open(hex_table, 'w')
  f.write("\n%s\n\n"%("i       Angle(Deg)      l       m       p       q       Angle(Rad)      L(Ang)"))
  for i in range(index):
    angle =  np.arccos((3*i**2 + 3*i + 0.5)/(3*i**2 + 3*i + 1))*180/np.pi
    f.write("%d\t%6.10f\t%d\t%d\t%d\t%d\t%6.10f\t%6.10f\n"%(i, angle, i, i+1, -1*(i+1), 2*i + 1, angle*np.pi/180, np.sqrt(3*i**2 + 3*i + 1)*np.float(alat)))
  f.close()

# run--
commensurate_least(200, 'hex.table')
