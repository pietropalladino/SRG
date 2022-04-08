#----------------------------------------------------------------------------------------------------------------------
# SRG.py

# author: Pietro Palladino
# date: 07/04/2022
# tested with Python v3.0

# SRG evolution of a two dimensional simplified potential and a chiral nucleon-nucleon realistic interaction for the deuteron
# bound states, with a uniform mesh for the momenta discretization
#----------------------------------------------------------------------------------------------------------------------
import myfunctions 
from myfunctions import uniform_weights, read_mesh, read_interaction, plot_snapshots_delta, plot_snapshots_chiral, commutator
from myfunctions import np, array, dot, diag, reshape, sqrt, pi, eigvalsh, ode, interpolate, plot_3d_chiral, sys
#----------------------------------------------------------------------------------------------------------------------
# Define the SRG flow differential equation
#----------------------------------------------------------------------------------------------------------------------
def derivative(lam, y, T):

 dim=T.shape[0]

 V=reshape(y, (dim,dim)) 

 eta=commutator(T,V) 

 dV=-4.0/(lam**5)*commutator(eta, T+V) 

 dydt=reshape(dV, -1) 

 return dydt
#----------------------------------------------------------------------------------------------------------------------
# Define the simplified potential delta and solve the previous ODE
#----------------------------------------------------------------------------------------------------------------------

def Delta():
 
 print('The potential in momentum space is constant and of the form: alpha/2*pi\n')
 print('choose a value for alpha:\n')
 
 alpha=float(input())
 momenta= read_mesh('/home/utente/SRG/docs/momenta.txt') 
 weights=uniform_weights(momenta)
 dim= len(momenta)
 T=diag(momenta*momenta) 
 V_delta=np.ones_like(T)
 V_delta*=-(alpha)/(2*pi)
 weight_matrix=np.zeros_like(T)

# Define a matrix which depends on the weights and the momenta

 for i in range(dim):
  for j in range(dim):
   qiqj=momenta[i]*momenta[j]
   weight_matrix[i,j]=qiqj*sqrt(weights[i]*weights[j])

 V_delta*=weight_matrix
 y0=reshape(V_delta,-1) 
 
 lam_initial =10.0
 lam_final = 0.0

# Solve explicitly the differential equations using a backward differentiation method
 
 solver = ode(derivative,jac=None)
 solver.set_integrator('vode', method='bdf', order=5, nsteps=1000)
 solver.set_f_params(T)
 solver.set_initial_value(y0, lam_initial)

 flowparams=([lam_initial])
 V_deltas=([V_delta])

# Use an adaptive step size method

 while solver.successful() and solver.t > lam_final:
   
   if solver.t >= 5.0:
     ys = solver.integrate(solver.t-1.0)
   elif solver.t >= 2.5:
     ys = solver.integrate(solver.t-0.5)
   elif solver.t < 2.5 and solver.t >= lam_final:
     ys = solver.integrate(solver.t-0.1)
   flowparams.append(solver.t)
   V_matrix= reshape(ys,(dim,dim))
   V_deltas.append(V_matrix)
   
# Print all the eigenvalues at each step 
 
 for s in range(0,len(flowparams)) :
  
  sourcefile=open('/home/utente/SRG/results/simplified_potential/Delta_Eigenval/Eigenvalues at lambda=%s'%flowparams[s],'w')
  print(eigvalsh((T+V_deltas)[s]), file=sourcefile)
 sourcefile=open('/home/utente/SRG/results/simplified_potential/Delta_Eigenval/Eigenvalues differences', 'w')
 for s in range(0,len(flowparams)-1):
  a=eigvalsh((T+V_deltas)[s+1])
  b=eigvalsh((T+V_deltas)[s])
  r=a-b
  print(r,file=sourcefile)
# Print the results in a sequence of plots in a pdf file
 
 plot_snapshots_delta((V_deltas[0:]/weight_matrix), flowparams[0:], momenta, alpha)
#--------------------------------------------------------------------------------------------------------------------------
# Define the realistic Chiral potential for nuclear interaction
#--------------------------------------------------------------------------------------------------------------------------

def Chiral():

# Define a constant for scattering units (hbarm=(hbar)^{2}/m_N) where hbar=6.582E-22 MeV*s is the Planck constant and m_N=939 MeV/c^{2} is the nucleon mass
# hbar2_m is expressed in MeV*fm^{2}

 hbar2_m=41.47
 
 mom=read_mesh('/home/utente/SRG/docs/momenta.txt')
 momenta=np.concatenate([mom,mom])
 weights=uniform_weights(momenta)
 dim=len(momenta)
 T=diag(momenta*momenta)
  
 partial_waves=[]

# Define the potential matrix by composition of submatrices related to the channels

 for filename in ('/home/utente/SRG/docs/s1.txt','/home/utente/SRG/docs/d1.txt','/home/utente/SRG/docs/sd1.txt'):
  partial_waves.append(read_interaction(filename)) 

 V= np.vstack((np.hstack((partial_waves[0], partial_waves[2])), 
    np.hstack((np.transpose(partial_waves[2]), partial_waves[1])) 
    ))
# Define a matrix depending on the weights and the momenta 
 weight_matrix = np.zeros_like(T)
 for i in range(dim):
   for j in range(dim):
     
     qiqj = momenta[i]*momenta[j]
     weight_matrix[i,j] = qiqj*sqrt(weights[i]*weights[j])

 V *= weight_matrix/hbar2_m
 y0  = reshape(V, -1) 
 
               
 lam_initial = 20.0
 lam_final = 1.5

# Solve explicitly the differential equations using a backward differentiation method
 
 solver = ode(derivative,jac=None)
 solver.set_integrator('vode', method='bdf', order=5, nsteps=1000)
 solver.set_f_params(T)
 solver.set_initial_value(y0, lam_initial)


 flowparams=([lam_initial])
 Vs=([V])

 while solver.successful() and solver.t > lam_final:
   
   if solver.t >= 6.0:
     ys = solver.integrate(solver.t-1.0)
   elif solver.t >= 2.5:
     ys = solver.integrate(solver.t-0.5)
   elif solver.t < 2.5 and solver.t >= lam_final:
     ys = solver.integrate(solver.t-0.1)
        
   flowparams.append(solver.t)
   Vm = reshape(ys,(dim,dim))
   Vs.append(Vm)
 
# Print all the eigenvalues at each step 

 for s in range(0,len(flowparams)):
  sourcefile=open('/home/utente/SRG/results/realistic_potential/SRG_chiral_flow/Chiral_Eigenval/Eigenvalues at lambda=%s'%flowparams[s],'w')
  print(eigvalsh((T+Vs)[s])*hbar2_m, file=sourcefile)
 sourcefile=open('/home/utente/SRG/results/realistic_potential/SRG_chiral_flow/Chiral_Eigenval/Eigenvalues differences','w')
 for s in range(0,len(flowparams)-1):
  a=eigvalsh((T+Vs)[s+1])*hbar2_m
  b=eigvalsh((T+Vs)[s])*hbar2_m
  r=a-b
  print(r,file=sourcefile)
  
# Print the results in a sequence of 3D plots and 2D projections on the momenta plane

 plot_snapshots_chiral((Vs[0:]/weight_matrix), flowparams[0:], momenta, 4.0)
 plot_3d_chiral((Vs[0:]/weight_matrix), flowparams[0:], momenta, 4.0)
 
#----------------------------------------------------------------------------------------------------------------------------------
# Main program
#----------------------------------------------------------------------------------------------------------------------------------

def main():
 
 def question_block(): 
  print('Which potential do you choose?:(digit A\B)\n\n A) 2d delta function (simplified)  \n B) Chiral NN potential (realistic) \n')
  x=input()

  if x== 'A':
   Delta()

  elif x== 'B':
   Chiral()

 question_block()
 print('Do you want to finish?:(digit yes/no)')
 y=input()

 while y == 'no':
  question_block()
  print('Do you want to finish?:(digit yes/no)')
  y=input()

 if y == 'yes': 
  exit(0)
 

if __name__ == "__main__": 
 main()

 



