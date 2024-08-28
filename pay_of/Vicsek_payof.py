import numpy as np 
from scipy.spatial import cKDTree
import scipy as sp
from scipy import sparse

#Some of the from 
#https://francescoturci.net/2020/06/19/minimal-vicsek-model-in-python/
#Has been used

def Vicsek_payof(v_0,r, N,L,iterations,payof_matrix):
        theta=np.zeros([iterations,N])
        pos=np.zeros([iterations,N,2])
       

        initial_pos= [ [np.random.uniform(0,L),np.random.uniform(0,L) ] for k in range(N) ]
        pos[0]=initial_pos
        initial_theta=np.random.uniform(-np.pi,np.pi,N)
        theta[0]=initial_theta
        
        #Creating the s
        initial_s=np.array([ [1,0] if k%2==0  else [0,1] for k in range(N) ])
        s=np.array([ initial_s for i in range(iterations)])
        

        for iteration in range(iterations-1):
                #Updating the postions 
                pos[iteration + 1] = (pos[iteration] + v_0 * np.stack([np.cos(theta[iteration]), np.sin(theta[iteration])])) 
               
                


                
                #updating the angles
                tree = cKDTree(pos[iteration])
                negibors=tree.query_ball_point(pos[iteration],r ) 
                pay_of_list=np.zeros(N)


                dist = tree.sparse_distance_matrix(tree, max_distance=r,output_type='coo_matrix')
 
               #important 3 lines: we evaluate a quantity for every column j
                data = np.exp(theta[iteration][dist.col]*1j)
               # construct  a new sparse marix with entries in the same places ij of the dist matrix
                neigh = sparse.coo_matrix((data,(dist.row,dist.col)), shape=dist.get_shape())
                # and sum along the columns (sum over j)
                complex_sum = np.squeeze(np.asarray(neigh.tocsr().sum(axis=1)))  
                
                theta[iteration+1] = np.angle(data+ complex_sum)
 


                for particle_index in range(N):
                      negibor_indexes=negibors[particle_index]
                      #The angles of the neigbors
                     # angles=theta[iteration][negibor_indexes]
                      #theta[iteration+1][particle_index]=np.angle( np.exp(theta[iteration][particle_index]*1j )+sum(np.exp(angles*1j) ) )
                             


                      #Making sure that there are neigbors  
                      if len(negibor_indexes)!=0:
                        #Cacluting the payof   with the sum s_iAs_j  
                           neighbor_s=s[iteration][negibor_indexes]                                                       #The number of neighbors                                                                       
                           P=sum(   [np.dot(s[iteration][particle_index] ,np.dot(payof_matrix,k) ) for k in neighbor_s ])/len(negibor_indexes)
                           #If the length is zero the pay of will be zero due to np zeros 
                           pay_of_list[particle_index]=P
                           
                    
                #Now updating strategies based on neigbor with highest payof
                for particle_index in range(N):
                     negibor_indexes=negibors[particle_index]
                     pay_of_among_neigbors=pay_of_list[negibor_indexes]
                     #The case where max pay of among neigbors is greater than its own payof
                     if len(negibor_indexes)!=0:
                       max_index = np.argmax(pay_of_among_neigbors)  # Use pay_of_among_neigbors directly
                       max_neighbor_pay_of = pay_of_among_neigbors[max_index] 
                       if max_neighbor_pay_of>pay_of_list[particle_index]:
                           s[iteration + 1][particle_index] = s[iteration][negibor_indexes[max_index]]
                           
                     
               
        return pos,theta,s                      





def Vicsek_PD(v_0,r, N,L,iterations,b):
     payof_matrix=np.array( [[1,0],
                              [b,0]]
     )
     return Vicsek_payof(v_0,r, N,L,iterations,payof_matrix)

def Vicsek_SD(v_0,r, N,L,iterations):
     payof_matrix=np.array([ [1,1-r],
                              [1+r,0]]
     )
     return Vicsek_payof(v_0,r, N,L,iterations,payof_matrix)



