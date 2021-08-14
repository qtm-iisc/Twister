
# -*- coding: utf-8 -*-
"""
Created on Wed May 30 17:07:57 2018

@author: saism
"""
import time
import numpy as np
import matplotlib.pyplot as plt
from math import *
from itertools import product,combinations
from mpi4py import MPI as mpi
import sys
import os
import shutil
from scipy.linalg import polar
comm = mpi.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

#Combines arrays
def join(A,B):
  if len(A)!=0:
    A=np.concatenate((A,B),axis=0)
  else:
    A=B
  return A


#Plots array of vectors on a graph
def plotvec(V,a1, a2,Ar,ag):
  Ar=np.round(Ar)
  Ar_list=np.ndarray.tolist(np.unique(Ar))
  colors = plt.cm.autumn(np.linspace(0,1,len(Ar_list)))
  for i in range(len(V)):
    ind=Ar_list.index(Ar[i])
    V_i=np.asarray(V[i]).reshape(2,2)
    v1 = V_i[0][0]*a1 + V_i[0][1]*a2
    v2 = V_i[1][0]*a1 + V_i[1][1]*a2
    plt.arrow(0,0,v1[0],v1[1],width=0.25,head_width=1.0,length_includes_head=True,alpha=0.75,edgecolor=colors[ind],facecolor=colors[ind])
    plt.arrow(0,0,v2[0],v2[1],width=0.25,head_width=1.0,length_includes_head=True,alpha=0.75,edgecolor=colors[ind],facecolor=colors[ind])
  plt.axis('equal') 
  plt.xlim(-max(max(v1)-2,max(v2)),max(max(v1),max(v2))+2)
  plt.ylim(-max(max(v1),max(v2))-2,max(max(v1),max(v2))+2)
  plt.title('Twist angle='+str(ag)+'\N{DEGREE SIGN}')


#Generates cell of lattice points
def gen_cell(n1,n2,a1,a2,b):
    nb=len(b)
    X=[]
    for i in range(-1*n1//2,n1//2):
        for j in range(-1*n2//2, n2//2):
            for k in range(nb):
                X.append([b[k][0]/n1+np.float(i)/n1,b[k][1]/n2+np.float(j)/n2])

    posit=np.array(X)

    return posit


#Crystal coordinated coverted to cartesian coordinates
def crys_to_ang(X,a1,a2):
    x=[]
    for i in range(len(X)):
        x.append(X[i][0]*a1 + X[i][1]*a2)
    p=np.array(x)
    return p

#Rotates vector x by theta
def rotate(x,theta):
    R=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
    for i in range(len(x)):
        x[i]=np.dot(R,x[i].T).T
    return x


#Main function
def  clt(a1,a2,b1,b2,theta_i,theta_f,dth,mismatch,strain_per,strain_layer,fix_ang,f_ang,range_nm_l,range_nm_u,th_er,alat_a,alat_b,nba_a,nba_b,plot,strain_tv,deepsearch):
  mism_ang={}
  new_mis_list={}
  min_str_list={}
  Sl_dict={}
  new_alat_a_dict={}
  new_alat_b_dict={}
  ind_list={}
  total_atom_list={}
  new_ag_list={}
  new_uv_list={}
  st_t = time.time()
  if strain_tv == 'True' and rank == 0:
  	print("Note: Since strain_tensor_vector is True, only the top layer will be strained regardless of strain_layer")
    
  range_nm=np.arange(range_nm_l,range_nm_u,1)  #Range of n,m to be scattered over prcoessors
  numDataPerRank_s1=int(len(range_nm)/size) #Array size for each rank
  data_s1=None
  
  ##Data scatter and receive 
  if numDataPerRank_s1==0:
    print('Too many processors used.Number of processors must be less than or equal to ',len(range_nm))
  else:
    if rank == 0:
      data_s1 = range_nm[:numDataPerRank_s1*size]  #Define array for recollecting after scatter
    recvbuf = np.empty(numDataPerRank_s1, dtype='int64') #Allocate space for buffer object for recieving 
    comm.Scatter(data_s1, recvbuf, root=0) #scatter numDataPerRank_s1*size to processors
    data_s2=None      
    
    if len(range_nm[numDataPerRank_s1*size:])!=0: 
      for i in range((len(range_nm[numDataPerRank_s1*size:])%numDataPerRank_s1)): #divide remainder of range_nm and send 
        if rank == 0:
          data_s2 = range_nm[numDataPerRank_s1*size+i]
          comm.send(data_s2, dest=i)  #Send data to rank i processor
        elif rank == i:
          data_s2 = comm.recv(source=0) #Receive data at node i processor
    if data_s2!=None:    
         data=np.insert(recvbuf,0,data_s2) 
    else:
      data=recvbuf
    ##Generate permutations of sets of (n1,n2,m1) from the data recieved with repetitions generated at each processor 
    perm_nm=np.empty((0))
    for i in data:  #permutations for each element in data 
      B=np.asarray(list(product(range_nm,repeat=3))) 
      A=i*np.ones((len(B),4)) #i is added as m2 to each permutation
      A[:,:-1]=B
      perm_nm=np.unique(join(perm_nm,A),axis=0) #Completed set of permutations of (n1,n2,m1,m2) for rank i processor
      
    ##Filter (n1,n2,m1,m2) permutaions by mismatch for each angle 
    for ag in np.arange(theta_i,theta_f,dth): 
      init_l1 = time.time()
      sl=[]
      tolr_vec1=[]
      tolr_vec2=[]
      tolr_avg=[]
      new_alat_a=[]
      new_alat_b=[]
      Strain_avg=[]
      sl_m=[]
      ang=ag*np.pi/180
      SL=[]
      SL_m=[]
      area=[]
      
      #Rotated unit cell vectors
      R=np.array([[np.cos(ang),-np.sin(ang)],[np.sin(ang),np.cos(ang)]]) #anti-clockwise rotation
      b2_t=b2.T
      b2_r=(np.dot(R,b2_t)).T 
      b1_t=b1.T
      b1_r=(np.dot(R,b1_t)).T
      
      Area_avg=[]
      A =np.array([a1,a2]) #Lattice 1 unit cell vectors
      B=np.array([b1_r,b2_r]) #Rotated latiice 2 unit cell vectors
      st1 = time.time()
      strain=[]
      cos_angle=[]
      new_alat_a=[]
      new_alat_b=[]
  
      a_nm=np.matmul(perm_nm[:,:2],A) #Lattice 1 superlattice vector (n1*a1 + n2*a2)
      b_nm=np.matmul(perm_nm[:,2:],B) #Lattice 2 supperlattice vector (m1*b1_r + m2*b2_r)
      A_B=a_nm-b_nm

      SL_ang=[]
      total_atoms=0
      new_mis_vec1=[]
      new_mis_vec2=[]
      strd_uv=[]
      p=[]
      u=[]
      new_ag=[]
      new_uv=[]
      tem_mis=[]
  
      with np.errstate(divide='ignore',invalid='ignore'): 
        a_a=np.linalg.norm(a_nm,axis=1)
        AB_AB=np.linalg.norm(A_B,axis=1)
        str_nm=AB_AB                #mismatch on overlap (|A-B|)
        
        #filter vectors by mismatch tolerance
        a_nm=a_nm[str_nm<mismatch]  
        b_nm=b_nm[str_nm<mismatch]
        strain=str_nm[str_nm<mismatch]
        sl=perm_nm[:,:2][str_nm<mismatch]
        sl_m=perm_nm[:,2:][str_nm<mismatch] 
      sl_minAr = []
      Area = []
      
      #send arrays back to rank 0 processor
      if rank!=0: 
         size_arr=np.array(len(a_nm)) 
         comm.Send(size_arr,dest=0,tag=5)
         comm.Send(a_nm,dest=0,tag=0)
         comm.Send(b_nm,dest=0,tag=1)
         comm.Send(strain,dest=0,tag=2)
         comm.Send(sl,dest=0,tag=3)
         comm.Send(sl_m,dest=0,tag=4)
      #concatenate recieved data 
      elif rank==0:
        for j in range(1,size):  
            size_arr=np.empty([1,1],dtype='int64')
            comm.Recv(size_arr,source=j,tag=5)
            data_a_nm=np.empty([size_arr[0][0],2])
            comm.Recv(data_a_nm,source=j,tag=0)
            if len(a_nm)!=0:
              a_nm=np.concatenate((a_nm,data_a_nm),axis=0)
            else:
              a_nm=data_a_nm
            data_b_nm=np.empty([size_arr[0][0],2])
            comm.Recv(data_b_nm,source=j,tag=1)
            if len(b_nm)!=0:
              b_nm=np.concatenate((b_nm,data_b_nm),axis=0)
            else:
              b_nm=data_b_nm
            data_strain=np.empty(size_arr[0][0])
            comm.Recv(data_strain,source=j,tag=2)
            if len(strain)!=0:
              strain=np.concatenate((strain,data_strain),axis=0)
            else:
              strain=data_strain
            data_sl=np.empty([size_arr[0][0],2])
            comm.Recv(data_sl,source=j,tag=3)
            if len(sl)!=0:
              sl=np.concatenate((sl,data_sl),axis=0)
            else:
              sl=data_sl
            data_sl_m=np.empty([size_arr[0][0],2]) 
            comm.Recv(data_sl_m,source=j,tag=4)
            if len(sl_m)!=0:
              sl_m=np.concatenate((sl_m,data_sl_m),axis=0)
            else:
              sl_m=data_sl_m

      if rank==0:
          
          if len(a_nm)>1:
            
            #generate all possible pairs of superlattice vectors  
            
            a_nm_comb=np.array(list(combinations(a_nm,2))) 
            a_nm_comb=a_nm_comb.reshape((len(a_nm_comb),-1))
            a_x_1 = np.hstack((a_nm_comb[:,:2], np.zeros((a_nm_comb[:,:2].shape[0], 1), dtype=a_nm.dtype)))
            a_x_2 = np.hstack((a_nm_comb[:,2:], np.zeros((a_nm_comb[:,2:].shape[0], 1), dtype=a_nm.dtype)))
            ar_sl=np.linalg.norm(np.cross(a_x_1,a_x_2,axis=1),axis=1)
            a_unit1=(np.linalg.norm(np.cross(a1,a2)))
            a_nm_comb=a_nm_comb[ar_sl>a_unit1]
            
       
            b_nm_comb=np.array(list(combinations(b_nm,2)))
            b_nm_comb=b_nm_comb.reshape(len(b_nm_comb),-1)
            b_nm_comb=b_nm_comb[ar_sl>a_unit1]
       
            sl_comb=np.array(list(combinations(sl,2)))
            sl_comb=sl_comb.reshape(len(sl_comb),-1)
            sl_comb=sl_comb[ar_sl>a_unit1]
       
            sl_m_comb=np.array(list(combinations(sl_m,2)))
            sl_m_comb=sl_m_comb.reshape(len(sl_m_comb),-1)
            sl_m_comb=sl_m_comb[ar_sl>a_unit1]
       
            strain_comb=np.array(list(combinations(strain,2)))
            strain_comb=strain_comb.reshape(len(strain_comb),-1)
            strain_comb=strain_comb[ar_sl>a_unit1]
           
             
       
            ar_sl=ar_sl[ar_sl>a_unit1]
            
       
            magn_A=np.linalg.norm(a_nm_comb[:,:2],axis=1)
            magn_B=np.linalg.norm(a_nm_comb[:,2:],axis=1)
          
            #check for angle between vectors if specified
            if fix_ang=='True':
              ij_dot_cos=np.divide(np.einsum('ij, ij->i', a_nm_comb[:,:2], a_nm_comb[:,2:]),magn_A*magn_B)
              
              ij_dot=abs(ij_dot_cos-np.cos(f_ang))
              
              a_nm_comb=a_nm_comb[ij_dot<th_er]
              reshape_anm=[]
              reshape_bnm=[]
              strd_uv=[]
              symm_def=[]
              symm_def_uv=[]
              sl_comb=sl_comb[ij_dot<th_er]
              sl_m_comb=sl_m_comb[ij_dot<th_er]
              b_nm_comb=b_nm_comb[ij_dot<th_er]
              strain_comb=strain_comb[ij_dot<th_er]
              
              ar_sl=ar_sl[ij_dot<th_er]
              ij_dot_cos=ij_dot_cos[ij_dot<th_er]
                            

                      
                  
     
            else:
              ij_dot=np.divide(np.einsum('ij, ij->i', a_nm_comb[:,:2], a_nm_comb[:,2:]),magn_A*magn_B)
       
     
            a_nm_i=np.linalg.norm(a_nm_comb[:,:2],axis=1) 
            b_nm_i=np.linalg.norm(b_nm_comb[:,:2],axis=1)


            a_nm_j=np.linalg.norm(a_nm_comb[:,2:],axis=1)
            b_nm_j=np.linalg.norm(b_nm_comb[:,2:],axis=1)
            
           #Straining lattice parameters
            if strain_tv=='False': 
              if strain_layer=='Both': #Strain lattice parameters of both layers                
                tolr_1=np.divide(b_nm_i-a_nm_i,b_nm_i+a_nm_i)
                tolr_2=np.divide(b_nm_j-a_nm_j,b_nm_j+a_nm_j)

                avg=((tolr_1+tolr_2)/2)
                bol_str=avg<strain_per #Filter by tolerated strain percentage
                avg=avg[bol_str]
                new_alat_a=np.dstack((alat_a[0]*(1+avg),alat_a[1]*(1+avg)))[0] #New lattice parameters
                new_alat_b=np.dstack((alat_b[0]*(1-avg),alat_b[1]*(1-avg)))[0]
            
              elif strain_layer=='Bottom':  #Strain lattice parameters of bottom layer
                tolr_1=np.divide(b_nm_i-a_nm_i,a_nm_i)
                tolr_2=np.divide(b_nm_j-a_nm_j,a_nm_j) 

                avg=((tolr_1+tolr_2)/2)
                bol_str=avg<strain_per
                avg=avg[bol_str]            #Filter by tolerated strain percentage
                new_alat_a=np.dstack((alat_a[0]*(1+avg),alat_a[1]*(1+avg)))[0] #New lattice parameters
                new_alat_b=np.dstack((alat_b[0]*(1-0*avg),alat_b[1]*(1-0*avg)))[0]
              else :                       #Strain lattice parameters of top layer
                tolr_1=np.divide(a_nm_i-b_nm_i,b_nm_i)
                tolr_2=np.divide(a_nm_j-b_nm_j,b_nm_j)  
                avg=((tolr_1+tolr_2)/2)
                bol_str=avg<strain_per
                avg=avg[bol_str]           #Filter by tolerated strain percentage
                new_alat_a=np.dstack((alat_a[0]*(1+0*avg),alat_a[1]*(1+0*avg)))[0]
                new_alat_b=np.dstack((alat_b[0]*(1+avg),alat_b[1]*(1+avg)))[0] #New lattice parameters
              tolr_vec1=(tolr_1[bol_str])
              tolr_vec2=(tolr_2[bol_str])
              a_nm_comb=a_nm_comb[bol_str]
              sl_comb=sl_comb[bol_str]
              b_nm_comb=b_nm_comb[bol_str]
              cos_angle=ij_dot_cos[bol_str]
              Area_avg=ar_sl[bol_str]
              SL=np.hstack((sl_comb[:,:2],sl_comb[:,2:]))
              SL_m=np.hstack((sl_m_comb[:,:2],sl_m_comb[:,2:]))
              Strain_avg=strain_comb[bol_str] 
              
            #Straining unit vectors
            elif strain_tv=='True':
                reshape_anm=[]
                reshape_bnm=[]
                reshape_slm_nm=[]
                reshape_sl_nm=[]
                tem_mis=[]
                strd_uv=[]
                bol_str=np.empty(len(a_nm_comb),dtype=bool)
                max_str=[]
                for k in range(len(b_nm_comb)):
                   reshape_anm.append(a_nm_comb[k].reshape(2,2))
                   reshape_bnm.append(b_nm_comb[k].reshape(2,2))
                   reshape_slm_nm.append(sl_m_comb[k].reshape(2,2))
                   reshape_sl_nm.append(sl_comb[k].reshape(2,2))
                        
                   strd_uv.append(np.dot(np.linalg.inv(np.array([b1,b2])),np.dot(np.linalg.inv(reshape_slm_nm[k]),np.dot(reshape_anm[k],np.linalg.inv(R.T)))))
                   bol_str[k]=(np.all(abs(abs(strd_uv[k])-np.identity(2))<strain_per)) #Filter by tolerated strain percentage
                   max_str.append(np.max(abs(abs(strd_uv[k])-np.identity(2))))
                   tem_mis.append(np.dot(reshape_slm_nm[k],np.dot(np.array([b1,b2]),np.dot(strd_uv[k],R.T)))-reshape_anm[k])
                   u_strd,p_strd=polar(strd_uv[-1],side='right')
                   symm_def.append(np.dot(np.array([b1_r,b2_r]),p_strd))
                   symm_def_uv.append(np.dot(np.array([b1,b2]),p_strd))
                symm_def=np.asarray(symm_def)
                strd_uv=np.asarray(strd_uv)
                symm_def_uv=np.asarray(symm_def_uv)
                strd_uv=strd_uv[bol_str]  #Strained unit vectors

                if len(strd_uv)>0:
                  new_uv_r=symm_def[bol_str]
                  new_uv=symm_def_uv[bol_str]
                  tem_mis=np.asarray(tem_mis)
                  reshape_anm=np.asarray(reshape_anm)
                  reshape_bnm=np.asarray(reshape_bnm)
                  reshape_sl_nm=np.asarray(reshape_sl_nm)
                  reshape_slm_nm=np.asarray(reshape_slm_nm)
                  tem_mis=tem_mis[bol_str]
                  reshape_anm=reshape_anm[bol_str]
                  reshape_bnm=reshape_bnm[bol_str]
                  SL=reshape_sl_nm[bol_str]
                  SL_m=reshape_slm_nm[bol_str]
                  Area_avg=ar_sl[bol_str]
                  cos_angle=ij_dot_cos[bol_str]
                  Strain_avg=strain_comb[bol_str]   
                  max_str=np.asarray(max_str)[bol_str]

          
          else:

            continue                      
          if len(Area_avg)==0:    
            continue
          else:
            #Output file
            if os.path.exists("solutions")==False:
              os.system("mkdir solutions")
            ag=round(ag,abs(int(np.log10(float(dth)))))
            f=open("solutions/Angle_"+str(ag)+".out","w")
            f.write("Twist angle:")
            f.write("%f\n" %ag)
          
            for ind in range(len(Area_avg)):
              da=min(Area_avg)
              
              #minimum area filter
              if ((Area_avg[ind])-min(Area_avg)) < da :  
                sl_minAr.append(SL[ind])
                Area.append(Area_avg[ind])
                #Sl_dict[ag]=SL
                if (ag in Sl_dict) == False:
                    Sl_dict[ag]=[SL[ind]]
                else:
                   Sl_dict[ag].append(SL[ind])
                new_alat_a_dict[ag]=new_alat_a
                new_alat_b_dict[ag]=new_alat_b
                if (ag in ind_list) == False:
                    ind_list[ag]=[ind]
                else:
                   ind_list[ag].append(ind)    
                #Strained lattice paramters
                if strain_tv=='False':
                  
                   #Strained unit cell vectors
                   new_a1=a1*new_alat_a[ind][0]/(alat_a[0])  
                   new_a2=a2*new_alat_a[ind][1]/alat_a[1]
                   new_b1_r=b1_r*new_alat_b[ind][0]/(alat_b[0])
                   new_b2_r=b2_r*new_alat_b[ind][1]/alat_b[1]
                   
                   #Area of unit cell vectors 
                   a_unit1=(np.linalg.norm(np.cross(new_a1,new_a2)))
                   a_unit2=(np.linalg.norm(np.cross(new_b1_r,new_b2_r)))
                    
                   #Strained superlattice vectors 
                  
                   #Lattice 1
                   new_v1=SL[:,:2][ind][0]*new_a1+SL[:,:2][ind][1]*new_a2
                   new_v2=SL[:,2:][ind][0]*new_a1+SL[:,2:][ind][1]*new_a2
                   #Lattice 2
                   new_v1_m=SL_m[:,:2][ind][0]*new_b1_r+SL_m[:,:2][ind][1]*new_b2_r
                   new_v2_m=SL_m[:,2:][ind][0]*new_b1_r+SL_m[:,2:][ind][1]*new_b2_r
         
                   #Number of atoms
                   A_fac1=round(nba_a[0]*np.linalg.norm(np.cross(new_v1,new_v2))/a_unit1) #Number of atoms in superlattice from lattice 1
                   A_fac2=round(nba_b[0]*np.linalg.norm(np.cross(new_v1_m,new_v2_m))/a_unit2) #Number of atoms in superlattice from lattice 2               
                   total_atoms=A_fac1+A_fac2

                   if (ag in total_atom_list) == False:
                      total_atom_list[ag]=[total_atoms]
                   else:
                      total_atom_list[ag].append(total_atoms)
                
                   #Angle between superllatice vectors
                   SL_ang=(acos(cos_angle[ind])*180/np.pi)
                
                   #Mismtatch between strained layers
                   new_mis_vec1= new_v1-(new_v1_m)
                   new_mis_vec2= new_v2-(new_v2_m)
                   if (ag in new_mis_list)== False :
                      new_mis_list[ag]=[min(np.linalg.norm(new_mis_vec1),np.linalg.norm(new_mis_vec2))]
                   else:
                      new_mis_list[ag].append(min(np.linalg.norm(new_mis_vec1),np.linalg.norm(new_mis_vec2)))

                #Strained unit vectors  
                elif strain_tv=='True':  
                     new_b1_r=new_uv_r[ind][0]
                     new_b2_r=new_uv_r[ind][1]
                     
                     #Area of unit cell vectors 
                     a_unit1=(np.linalg.norm(np.cross(a1,a2))) 
                     a_unit2=(np.linalg.norm(np.cross(new_b1_r,new_b2_r)))
                     #Area of superlattice
                     Ar_sl=np.linalg.norm(np.cross(reshape_anm[ind][0],reshape_anm[ind][1]))
                     
                     #Number of atoms
                     A_fac1=round((Ar_sl/a_unit1)*nba_a[0]) #Number of atoms in superlattice from lattice 1
                     A_fac2=round((Ar_sl/a_unit2)*nba_b[0]) #Number of atoms in superlattice from lattice 2  
                     total_atoms=A_fac1+A_fac2
                     if (ag in total_atom_list) == False:
                      total_atom_list[ag]=[total_atoms]
                     else:
                      total_atom_list[ag].append(total_atoms)
                     

                     #Angle between superllatice vectors
                     SL_ang=(acos(np.dot(new_uv[ind][0],new_uv[ind][1])/(np.linalg.norm(new_b1_r)*np.linalg.norm(new_b2_r)))*180/np.pi)
                     if len(min_str_list)==0:
                        min_str=min(max_str)
                        min_str_list.update({ag:min_str})

                     elif min(max_str)< min_str:
                        min_new_mis=min(max_str)
                        min_str_list.update({ag:min_str})
 
    
                f.write("\n")

                #minimum mismatch/strain for each angle
                
                if len(mism_ang)==0:
                     mism_ang[ag]=(min(np.max(abs(Strain_avg),axis=1)))
                elif list(mism_ang.values())[-1]!= min(np.max(abs(Strain_avg),axis=1)) :
                     mism_ang.update({ag:min(np.max(abs(Strain_avg),axis=1))})

                #Output for straining lattice parameters
                if strain_tv=='False':
                    strain_tv_false_write(f,strain_layer,ind,SL,SL_m,total_atoms,SL_ang,Strain_avg,tolr_vec1,tolr_vec2, new_alat_a, new_alat_b,alat_a,alat_b,new_mis_vec1,new_mis_vec2)
                #Output for straining unit vectors
                elif strain_tv=='True':
                    u,p=polar(strd_uv[ind],side='right')
                    new_ag=np.arcsin(np.dot(u,R.T)[0,1])*180/np.pi
                    if (ag in new_ag_list)== False :
                      new_ag_list[ag]=[new_ag]
                    else:
                      new_ag_list[ag].append(new_ag)
                    
                    #new_uv
                    if (ag in new_uv_list)== False :
                        new_uv_list[ag]=[new_uv[ind]]
                    else:
                        new_uv_list[ag].append(new_uv[ind])


                    new_ag_rad=np.arcsin(np.dot(u,R.T)[0][1])
                    strain_tv_true_write(f,strain_layer,ind,SL,SL_m,total_atoms,SL_ang,Strain_avg,strd_uv,p,u,new_ag,new_uv,tem_mis,alat_a,alat_b)
            
            #plotting result vectors
            if plot=='True' :
    
                        
                  n1=int(max(range_nm_u,range_nm_l)*max(alat_a[0],alat_a[1],alat_b[0],alat_b[1]))
                  n2=n1
                  n1=60
                  n2=60
                  axis=[0,0]
                  basis1=[[0,0]]
                  basis2=[[0,0]]

                  l1=gen_cell(n1,n2,a1,a2,basis1)
                  l2=gen_cell(n1,n2,b1,b2,basis2)
                  
                  
                  # New supercell lattice vectors
                  A1 = n1*a1
                  A2 = n2*a2
                  B1 = n1*b1
                  B2 = n2*b2
                  L1=crys_to_ang(l1,A1,A2)
                  L2=crys_to_ang(l2,B1,B2)

                  # Define layer 1 and layer 2
                  L1_c= np.copy(L1)
                  L2_c=np.copy(L2)
                  L1_t=L1_c-axis[0]*a1-axis[1]*a2
                  L2_t=L2_c-axis[0]*b1-axis[1]*b2
                  L2_r=rotate(L2_t,ang)


                  tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
                          (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
                          (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
                          (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
                          (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]

                  for i in range(len(tableau20)):
                      r, g, b = tableau20[i]
                      tableau20[i] = (r / 255., g / 255., b / 255.)


                  leg=['Bottom Layer', 'Top Layer']

                  fig=plt.figure(figsize=(12,9))

                  plt.scatter(L1_t[:,0],L1_t[:,1],s = 40, color = 'darkorange')
                  plt.scatter(L2_r[:,0],L2_r[:,1],s = 40, color ='k')
                  plt.axis('off')
                  plt.legend(leg)
                  plotvec(np.array(sl_minAr),a1,a2,Area,ag)
                  plt.savefig('twist_angle_'+str(ag)+'.png',dpi=300)
  
  mism_ang=comm.bcast(mism_ang,root=0)
  Sl_dict=comm.bcast(Sl_dict,root=0)
  SL_m=comm.bcast(SL_m,root=0)
  total_atoms=comm.bcast(total_atoms,root=0)
  SL_ang=comm.bcast(SL_ang,root=0)
  Strain_avg=comm.bcast(Strain_avg,root=0)
  ind_list=comm.bcast(ind_list,root=0)
  total_atom_list=comm.bcast(total_atom_list,root=0)
  if len(mism_ang)!=0 and deepsearch=='True' and strain_tv=='False':
         new_mis_list=comm.bcast(new_mis_list,root=0)
         tolr_vec1=comm.bcast(tolr_vec1,root=0)
         tolr_vec2=comm.bcast(tolr_vec2,root=0)
         new_alat_a_dict=comm.bcast(new_alat_a_dict,root=0)
         new_alat_b_dict=comm.bcast(new_alat_b_dict,root=0)
         alat_a=comm.bcast(alat_a,root=0)
         alat_b=comm.bcast(alat_b,root=0)
         new_mis_vec1=comm.bcast(new_mis_vec1,root=0)
         new_mis_vec2=comm.bcast(new_mis_vec2,root=0)
         return mism_ang,new_mis_list,Sl_dict,SL_m,total_atoms,SL_ang,Strain_avg,tolr_vec1,tolr_vec2, new_alat_a_dict, new_alat_b_dict,alat_a,alat_b,new_mis_vec1,new_mis_vec2

  if len(mism_ang)==0:
      if rank==0:
        print('No solution found.')
        sys.exit()
      else:
        sys.exit()

  if strain_tv=='True' and len(mism_ang)!=0 : 
       if rank==0:
            min_mis=min(mism_ang.values())
            min_mis_pos=list(mism_ang.values()).index(min(mism_ang.values()))
            min_angle= list(mism_ang.keys())[min_mis_pos]
            min_angle=round(min_angle,abs(int(np.log10(float(dth)))))
            f=open("Angle_"+str(min_angle)+".out","w")
            f.write("Twist angle:")
            f.write("%f\n" %min_angle)
            shutil.copyfile('solutions/Angle_'+str(min_angle)+'.out','Angle_'+str(min_angle)+'.out')
            plot_lat=str("False")
            new_alat_a=0
            new_alat_b=0
            new_strain_avg=0
            twist_inp(strain_tv,strain_layer,new_uv_list[min_angle],new_ag_list,strd_uv,ind_list,total_atom_list,a1,a2,b1,b2,alat_a,alat_b,new_alat_a,new_alat_b,min_angle,Sl_dict,new_strain_avg,plot_lat)


  if len(new_mis_list.values())!=0  and deepsearch=='False' and strain_tv=='False':
         min_mis=(min(new_mis_list.values()))
         res = [key for key in new_mis_list if new_mis_list[key] == (min_mis)]
         min_angle=res[0]
         min_angle=round(min_angle,abs(int(np.log10(float(dth)))))
         f=open("Angle_"+str(min_angle)+".out","w")
         f.write("Twist angle:")
         f.write("%f\n" %min_angle)
         shutil.copyfile('solutions/Angle_'+str(min_angle)+'.out','Angle_'+str(min_angle)+'.out')
         plot_lat=str("False")
         new_strain_avg=new_mis_list[min_angle]
         new_uv=0
         strd_uv=0
         new_ag=0
         ind_list=0
         total_atom_list=0
         twist_inp(strain_tv,strain_layer,new_uv,new_ag,strd_uv,ind_list,total_atom_list,a1,a2,b1,b2,alat_a,alat_b,new_alat_a_dict,new_alat_b_dict,min_angle,Sl_dict,new_strain_avg,plot_lat)


      #   return mism_ang,new_mis_list,Sl_dict,SL_m,total_atoms,SL_ang,Strain_avg,tolr_vec1,tolr_vec2, new_alat_a_dict, new_alat_b_dict,alat_a,alat_b,new_mis_vec1,new_mis_vec2



def strain_tv_true_write(f,strain_layer,ind,SL,SL_m,total_atoms,SL_ang,Strain_avg,strd_uv,p,u,new_ag,new_uv,tem_mis,alat_a,alat_b):
    f.write("Superlattice (vector1 (n1,n2),vector2 (n1',n2')):\n")
    f.write("%3.9f %3.9f\n" %(SL[ind][0][0],SL[ind][0][1]))
    f.write("%3.9f %3.9f\n" %(SL[ind][1][0],SL[ind][1][1]))
    f.write("Superlattice (vector1 (m1,m2),vector2 (m1',m2')):\n")
    f.write("%3.9f %3.9f\n" %(SL_m[ind][0][0],SL_m[ind][0][1]))
    f.write("%3.9f %3.9f\n" %(SL_m[ind][1][0],SL_m[ind][1][1]))
    f.write("Total number of atoms:")
    f.write("%d\n" %total_atoms)
    f.write("Angle between super lattice vectors (in degree):")
    f.write("%f\n" %SL_ang)
    f.write("Mismatch (vector1,vector2) (in Angstrom):")
    f.write("%3.9f %3.9f\n" %(Strain_avg[ind][0],Strain_avg[ind][1]))
    f.write("Deformation tensor:\n ")
    f.write("%3.9f %3.9f\n" %((((strd_uv[ind]))[0][0]),((((strd_uv[ind]))[0][1]))))
    f.write("%3.9f %3.9f\n" %((((strd_uv[ind]))[1][0]),((((strd_uv[ind]))[1][1]))))
    f.write("Symmetric factor of deformation tensor:\n ")
    f.write("%3.9f %3.9f\n" %((((p))[0][0]),((((p))[0][1]))))
    f.write("%3.9f %3.9f\n" %((((p))[1][0]),((((p))[1][1]))))
    f.write("Rotation factor of deformation tensor:\n ")
    f.write("%3.9f %3.9f\n" %((((u))[0][0]),((((u))[0][1]))))
    f.write("%3.9f %3.9f\n" %((((u))[1][0]),((((u))[1][1]))))
    f.write("New angle of rotation (in degrees):")
    f.write("%3.12f\n" %new_ag)
    f.write("New unit vector of top layer (in Angstrom):\n ")
    f.write("%3.15f %3.15f\n" %(new_uv[ind][0][0]/alat_b[0],new_uv[ind][0][1]/alat_b[0]))
    f.write("%3.15f %3.15f\n" %(new_uv[ind][1][0]/alat_b[1],new_uv[ind][1][1]/alat_b[1]))
    f.write("Mismatch with new unit vectors (vector1,vector2) (in Angstrom):")
    f.write("%3.15f %3.15f\n" %(np.linalg.norm(tem_mis[ind][0]),np.linalg.norm(tem_mis[ind][1])))


def strain_tv_false_write(f,strain_layer,ind,SL,SL_m,total_atoms,SL_ang,Strain_avg,tolr_vec1,tolr_vec2, new_alat_a, new_alat_b,alat_a,alat_b,new_mis_vec1,new_mis_vec2):
    f.write("Superlattice (vector1 (n1,n2), vector2 (n1',n2')):\n")
    f.write("%3.9f %3.9f\n" %(SL[:,:2][ind][0],SL[:,:2][ind][1]))
    f.write("%3.9f %3.9f\n" %(SL[:,2:][ind][0],SL[:,2:][ind][1]))
    f.write("Superlattice (vector1 (m1,m2),vector2 (m1',m2')):\n")
    f.write("%3.9f %3.9f\n" %(SL_m[:,:2][ind][0],SL_m[:,:2][ind][1]))
    f.write("%3.9f %3.9f\n" %(SL_m[:,2:][ind][0],SL_m[:,2:][ind][1]))
    f.write("Total number of atoms: ")
    f.write("%d\n" %total_atoms)
    f.write("Angle between super lattice vectors (in degree): ")
    f.write("%f\n" %SL_ang)
    f.write("Mismatch (vector1,vector2) (in Angstrom): ")
    f.write("%3.9f %3.9f\n" %(Strain_avg[ind][0],Strain_avg[ind][1]))
  
    if strain_layer=='Both':
        f.write("Strain Percent (vector 1,vector2): ")
        f.write("%3.9f %3.9f\n" %(tolr_vec1[ind]*100,-1*tolr_vec2[ind]*100))
        f.write("New lattice parameter layer 1 (in Angstrom):")
        f.write("%3.9f %3.9f %3.9f\n" %(new_alat_a[ind][0],new_alat_a[ind][1],alat_a[2]))
        f.write("New lattice parameter layer 2 (in Angstrom):")
        f.write("%3.9f %3.9f %3.9f\n" %(new_alat_b[ind][0],new_alat_b[ind][1],alat_b[2]))
    elif strain_layer=='Bottom':
        f.write("Strain Percent in bottom layer (vector 1,vector2): ")
        f.write("%3.9f %3.9f\n" %(tolr_vec1[ind]*100,tolr_vec2[ind]*100))
        f.write("New lattice parameter layer 1 (in Angstrom):")
        f.write("%3.9f %3.9f %3.9f\n" %(new_alat_a[ind][0],new_alat_a[ind][1],alat_a[2]))
        f.write("Lattice parameter layer 2 (in Angstrom): ")
        f.write("%3.9f %3.9f %3.9f\n" %(new_alat_b[ind][0],new_alat_b[ind][1],alat_b[2]))
    else:
        f.write("Strain Percent in top layer (vector 1,vector2): ")
        f.write("%3.9f %3.9f\n" %(tolr_vec1[ind]*100,tolr_vec2[ind]*100))
        f.write("Lattice parameter layer 1 (in Angstrom):")
        f.write("%3.9f %3.9f %3.9f\n" %(new_alat_a[ind][0],new_alat_a[ind][1],alat_a[2]))

        f.write("New lattice parameter layer 2 (in Angstrom): ")
        f.write("%3.9f %3.9f %3.9f\n" %(new_alat_b[ind][0],new_alat_b[ind][1],alat_b[2]))
    f.write("Mismatch with new lattice parmaters (vector1,vector2) (in Angstrom): ")
    f.write("%3.9f %3.9f\n" %(np.linalg.norm(new_mis_vec1),np.linalg.norm(new_mis_vec2)))

def twist_inp(strain_tv,strain_layer,new_uv,new_ag_list,strd_uv,ind_list,total_atom_list,a1,a2,b1,b2,alat_a,alat_b,new_alat_a,new_alat_b,min_angle,SL,new_strain_avg,plot_lat):
    f=open("twist.inp","w")
    if strain_tv=='False':
       ind=np.argmin(new_strain_avg)
       SL=SL[min_angle]
       f.write("a1_l:\n")
       f.write("%3.9f %3.9f %3.9f\n" %(a1[0]/(alat_a[0]),a1[1]/(alat_a[1]),0))
       f.write("a2_l:\n")
       f.write("%3.9f %3.9f %3.9f \n" %(a2[0]/(alat_a[0]),a2[1]/(alat_a[1]),0))
       f.write("a3_l:\n")
       f.write("0.0 0.0 1.0\n")
       f.write("\n")
       f.write("celldm1_l, celldm2_l, celldm3_l: (Angstrom)\n")
       f.write("%3.9f %3.9f %3.9f\n" %(new_alat_a[min_angle][ind][0],new_alat_a[min_angle][ind][1],alat_a[2]))
       f.write("\n")
       f.write("a1_u:\n")
       f.write("%3.9f %3.9f %3.9f\n" %(b1[0]/alat_b[0],b1[1]/alat_b[1],0))
       f.write("a2_u:\n")
       f.write("%3.9f %3.9f %3.9f\n" %(b2[0]/alat_b[0],b2[1]/alat_b[1],0))
       f.write("a3_u:\n")
       f.write("0.0 0.0 1.0\n")
       f.write("\n")
       f.write("celldm1_u, celldm2_u, celldm3_u: (Angstrom)\n")
       f.write("%3.9f %3.9f %3.9f\n" %(new_alat_b[min_angle][ind][0],new_alat_b[min_angle][ind][1],alat_b[2]))
       f.write("\n")
       f.write("angle: (radians)\n")
       f.write("%3.9f\n" %(min_angle*np.pi/180))
       f.write("\n")
       f.write("layer2_from_file:\n")
       f.write("True basis_pos_crys_layer2\n")
       f.write("\n")
       f.write("translate_z: (Angstrom)\n")
       f.write("6.2\n")
       f.write("\n")
       f.write("Superlattice1: (m,n)\n")
       f.write("%3.9f %3.9f\n" %(SL[ind][:2][0],SL[ind][:2][1]))
       f.write("\n")
       f.write("Superlattice2: (p,q)\n")
       f.write("%3.9f %3.9f\n" %(SL[ind][2:][0],SL[ind][2:][1]))
       f.write("\n")
       f.write("Plot_lattice:\n")
       f.write(plot_lat)
    else:
           ind=np.argmin(total_atom_list[min_angle])
           new_ag=new_ag_list[min_angle][ind]
           new_uv=new_uv[ind]
           #ind=ind_list[min_angle][ind]
           SL=SL[min_angle]
           f.write("a1_l:\n")
           f.write("%3.9f %3.9f %3.9f\n" %(a1[0]/(alat_a[0]),a1[1]/(alat_a[0]),0))
           f.write("a2_l:\n")
           f.write("%3.9f %3.9f %3.9f \n" %(a2[0]/(alat_a[1]),a2[1]/(alat_a[1]),0))
           f.write("a3_l:\n")
           f.write("0.0 0.0 1.0\n")
           f.write("\n")
           f.write("celldm1_l, celldm2_l, celldm3_l: (Angstrom)\n")
           f.write("%3.9f %3.9f %3.9f\n" %(alat_a[0],alat_a[1],alat_a[2]))
           f.write("\n")
           f.write("a1_u:\n")
           f.write("%3.9f %3.9f %3.9f\n" %(new_uv[0][0]/alat_b[0],new_uv[0][1]/alat_b[0],0))
           f.write("a2_u:\n")
           f.write("%3.9f %3.9f %3.9f\n" %(new_uv[1][0]/alat_b[1],new_uv[1][1]/alat_b[1],0))
           f.write("a3_u:\n")
           f.write("0.0 0.0 1.0\n")
           f.write("\n")
           f.write("celldm1_u, celldm2_u, celldm3_u: (Angstrom)\n")
           f.write("%3.9f %3.9f %3.9f\n" %(alat_b[0],alat_b[1],alat_b[2]))
           f.write("\n")
           f.write("angle: (radians)\n")
           f.write("%3.12f\n" %(new_ag*np.pi/180))
           f.write("\n")
           f.write("layer2_from_file:\n")
           f.write("True basis_pos_crys_layer2\n")
           f.write("\n")
           f.write("translate_z: (Angstrom)\n")
           f.write("6.2\n")
           f.write("\n")
           f.write("Superlattice1: (m,n)\n")
           f.write("%3.9f %3.9f\n" %(SL[ind][0][0],SL[ind][0][1]))
           f.write("\n")
           f.write("Superlattice2: (p,q)\n")
           f.write("%3.9f %3.9f\n" %((SL[ind][1][0],SL[ind][1][1])))
           f.write("\n")
           f.write("Plot_lattice:\n")
           f.write(plot_lat)
    print("Solution found.")       


def deep_search(a1,a2,b1,b2,theta_i,theta_f,dth,mismatch,strain_per,strain_layer,fix_ang,f_ang,range_nm_l,range_nm_u,th_er,alat_a,alat_b,nba_a,nba_b,plot,strain_tv,deepsearch,f_thresh):
    if strain_tv=='True' :
      if rank==0:  
        print('Please enter DeepSearch as False.')
        sys.exit()
    else:
        mism_ang,new_mis_list,Sl_dict,SL_m,total_atoms,SL_ang,Strain_avg,tolr_vec1,tolr_vec2, new_alat_a_dict, new_alat_b_dict,alat_a,alat_b,new_mis_vec1,new_mis_vec2= clt(a1,a2,b1,b2,theta_i,theta_f,dth,mismatch,strain_per,strain_layer,fix_ang,f_ang,range_nm_l,range_nm_u,th_er,alat_a,alat_b,nba_a,nba_b,plot,strain_tv,deepsearch)    
        i=0
        dth=10**-(len(str(dth))-2)
        while min(min(list(new_mis_list.values())))>f_thresh:
          if i <10 :  
             i=1+i
             min_mis_list=[]
             for j in range(len(new_mis_list.values())):
                min_mis_list.append(min(list(new_mis_list.values())[j]))
             res = [key for key in new_mis_list if min(new_mis_list[key]) == min(min_mis_list)]
             min_angle=res[0]
             theta_i=min_angle - 10**-(str(dth).index('1') -1)
             theta_f=min_angle + 10**-(str(dth).index('1')-1)
             dth = np.float(10**-(str(dth).index('1')))
             dth=(f"{dth:.15f}")
             mismatch=mism_ang[min_angle]+float(dth)
             if float(dth)!=0:
                mism_ang,new_mis_list,Sl_dict,SL_m,total_atoms,SL_ang,Strain_avg,tolr_vec1,tolr_vec2, new_alat_a_dict, new_alat_b_dict,alat_a,alat_b,new_mis_vec1,new_mis_vec2= clt(a1,a2,b1,b2,theta_i,theta_f,float(dth),mismatch,strain_per,strain_layer,fix_ang,f_ang,range_nm_l,range_nm_u,th_er,alat_a,alat_b,nba_a,nba_b,plot,strain_tv,deepsearch)
             else:
                 print('No solutions found.')
                 sys.exit()

          else:
            print('No solutions found.')
            sys.exit()

        if rank==0:
         min_mis=(min(new_mis_list.values()))
         res = [key for key in new_mis_list if new_mis_list[key] == (min_mis)]
         min_angle=res[0] 
         min_angle=round(min_angle,abs(int(np.log10(float(dth)))))
         f=open("Angle_"+str(min_angle)+".out","w")
         f.write("Twist angle:")
         f.write("%s\n" %str(min_angle))
         shutil.copyfile('solutions/Angle_'+str(min_angle)+'.out','Angle_'+str(min_angle)+'.out')
         plot_lat=str("False")
         new_strain_avg=new_mis_list[min_angle]
         new_uv=0
         strd_uv=0
         new_ag=0
         ind_list=0
         total_atom_list=0
         twist_inp(strain_tv,strain_layer,new_uv,new_ag,strd_uv,ind_list,total_atom_list,a1,a2,b1,b2,alat_a,alat_b,new_alat_a_dict,new_alat_b_dict,min_angle,Sl_dict,new_strain_avg,plot_lat)

#reads user input from get_ang.inp
def read_input(input_file):
 fp = open(input_file,'r')
 lines = fp.readlines()
 fp.close()
 for i in range(len(lines)):
   if "range_nm:" in lines[i]:
      w = lines[i+1].split()
      range_nm = [eval(w[0]),eval(w[1])]

   if "Number_basis_atoms_a:" in lines[i]:
      w = lines[i+1].split()
      nba_a=[eval(w[0])]

   if "Number_basis_atoms_b:" in lines[i]:
      w = lines[i+1].split()
      nba_b=[eval(w[0])]

   if "theta_range:" in lines[i]:
     w = lines[i+1].split()
     theta_range= [eval(w[0]),eval(w[1]),eval(w[2])]

   if "a1:" in lines[i]:
     w = lines[i+1].split()
     a1 = [eval(w[0]),eval(w[1]),eval(w[2])]

   if "a2:" in lines[i]:
     w = lines[i+1].split()
     a2 = [eval(w[0]),eval(w[1]),eval(w[2])]
   if "a3:" in lines[i]:
     w = lines[i+1].split()
     a3 = [eval(w[0]),eval(w[1]),eval(w[2])]

   if "celldm_a:" in lines[i]:
     w = lines[i+1].split()
     celldm_a=[eval(w[0]),eval(w[1]), eval(w[2])]

   if "celldm_b:" in lines[i]:
     w = lines[i+1].split()
     celldm_b = [eval(w[0]),eval(w[1]), eval(w[2])]

   if "b1:" in lines[i]:
     w = lines[i+1].split()
     b1 = [eval(w[0]),eval(w[1]),eval(w[2])]

   if "b2:" in lines[i]:
     w = lines[i+1].split()
     b2= [eval(w[0]),eval(w[1]),eval(w[2])]

   if "b3:" in lines[i]:
     w = lines[i+1].split()
     b3 = [eval(w[0]),eval(w[1]),eval(w[2])]

   if "initial_mismatch_threshold (Angstrom):" in lines[i]:
     w=lines[i+1].split()
     mismatch=[eval(w[0])]
    
   if "final_mismatch_threshold (Angstrom):" in lines[i]:
     w=lines[i+1].split()
     f_thresh=[eval(w[0])] 

   if "strain_per:" in lines[i]:
     w=lines[i+1].split()
     strain_per=[eval(w[0])]

   if "strain_layer:" in lines[i]:
     w=lines[i+1].split()
     strain_layer=[eval(w[0])]
   
   if "strain_tensor_vector:" in lines[i]:
       w=lines[i+1].split()
       strain_tv=[eval(w[0])]

   if "fix_ang:" in lines[i]:
     w=lines[i+1].split()
     fix_ang=[eval(w[0])]

   if "f_ang:" in lines[i]:
     w=lines[i+1].split()
     f_ang=[eval(w[0])]

   if "plot:" in lines[i]:
     w=lines[i+1].split()
     plot=[eval(w[0])]

   if "DeepSearch:" in lines[i]:
     w=lines[i+1].split()
     deepsearch=[eval(w[0])]



 return np.array(range_nm),np.array(celldm_a),np.array(a1),np.array(a2),np.array(a3),np.array(celldm_b),np.array(b1),np.array(b2),np.array(b3),np.array(theta_range),np.array(mismatch),np.array(fix_ang),np.array(f_ang),np.array(nba_a),np.array(nba_b),np.array(strain_per),np.array(plot),np.array(strain_layer),np.array(strain_tv),np.array(deepsearch),np.array(f_thresh)


