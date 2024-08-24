import numpy as np
from scipy.spatial import cKDTree






def Vicsek_modified(v_0,eta,r, N,L,iterations,alpha,beta):
        theta=np.zeros([iterations,N])
        pos=np.zeros([iterations,N,2])
       

        initial_pos= [ [np.random.uniform(0,L),np.random.uniform(0,L) ] for k in range(N) ]
        pos[0]=initial_pos
        initial_theta=np.random.uniform(-np.pi,np.pi,N)
        theta[0]=initial_theta
        
        #Creating the s
        s=np.array([ np.random.choice([0,1],N ) for i in range(iterations)])
        

        for iteration in range(iterations-1):
                #Updating the postions 
                pos[iteration + 1] = (pos[iteration] + v_0 * np.stack([np.cos(theta[iteration]), np.sin(theta[iteration])], axis=-1)) % L
               
                


                
                #updating the angles
 
                tree = cKDTree(pos[iteration],boxsize=[L,L])
                negibors=tree.query_ball_point(pos[iteration],r ) 

                for particle_index in range(N):
                      negibor_indexes=negibors[particle_index]
                      #The angles of the neigbors
                      angles=theta[iteration][negibor_indexes]
                      if s[iteration][particle_index]==1:
                             angle_average=np.angle(  sum(np.exp(angles*1j) )    )
                             theta[iteration+1][particle_index]=angle_average+eta*np.random.uniform(-np.pi,np.pi)
                             
                      else:
                              #Choosing random angle
                              theta[iteration+1][particle_index]=np.random.uniform(0,2*np.pi)
                       
                     #Picking a random neigbor
                      random_neigbor=np.random.choice( negibor_indexes)
                      
                      #No need to update s if they are already the same. And also ensuring that there are neigbors
                      if s[iteration][particle_index]!=s[iteration][random_neigbor] and len(negibor_indexes)!=0 :
                        
                    
                        r1=sum(np.cos( 0.5*(theta[iteration][particle_index]-angles )) )
                      
                        #The cost
                        P1=r1-alpha*s[iteration][particle_index]*r/L

                        #Doing the same with its neigbor
                      
                        #Neigbor's neigbors
                        negibor_indexes2=negibors[random_neigbor]
                        angles2=theta[iteration][negibor_indexes2]
                        r2=sum(np.cos( 0.5*(theta[iteration][random_neigbor]-angles2 )) )
                        P2=r2-alpha*s[iteration][ random_neigbor]*r/L

                        #Making the probabilitstic choice
                        if np.random.uniform(0,1)<1/(1+np.exp((P1-P2)/beta)):
                             s[iteration+1][particle_index]=s[iteration][ random_neigbor]

               
        return pos,theta,s                      