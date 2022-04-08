import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import sys
import os
import scicomap as sc
import glob
import natsort
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages
from PyPDF2 import PdfFileMerger
from fpdf import FPDF
from matplotlib.pyplot import figure
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import AxesGrid, make_axes_locatable
from mpl_toolkits import mplot3d
from collections import defaultdict
from numpy import array, dot, diag, reshape, sqrt
from math import sqrt, pi
from scipy.linalg import eigvalsh
from scipy.integrate import ode
from scipy import interpolate
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar
from PIL import Image
from natsort import natsorted



#---------------------------------------------------------------------------------------------------------------------
# Define the commutator and some helper functions to extract data from files
#---------------------------------------------------------------------------------------------------------------------
def commutator(a,b):
 
 return dot(a,b)-dot(b,a)

def uniform_weights(momenta):

 weights=np.ones_like(momenta)
 weights *=abs(momenta[1]-momenta[0])

 return weights

def read_mesh(filename):
 
 data=np.loadtxt(filename)
 dim=data.shape[0]
 momenta=data[0:dim]

 return momenta
 
def read_interaction(filename):
 
 data=np.loadtxt(filename)
 dim=data.shape[1]
 V=data[0:,:dim]

 return V
#----------------------------------------------------------------------------------------------------------------------
# Define a function that prints plots containing 3D graphs and 2D projections of delta potential for SRG flow
#----------------------------------------------------------------------------------------------------------------------
def plot_snapshots_delta(Hs, flowparams, momenta, alpha):

 
 
 dim=len(momenta)
# Distinguish the maximum and the minimum according if alpha is positive or negative
 if alpha >0. :
  a=-alpha/(2*pi)
  b=alpha/(2*pi)
  levels=np.arange(a,b,0.05)
 if alpha <=0. :
  a=alpha/(2*pi)
  b=-alpha/(2*pi)
  levels=np.arange(a,b,0.05)
 x=np.linspace(0.0,7.0,141)
 y=np.linspace(0.0,7.0,141)
 X,Y=np.meshgrid(x,y)
#fix the colormap using scicomap package

 mpl_cmap_obj=plt.get_cmap('gist_rainbow')
 div_map=sc.ScicoDiverging(cmap=mpl_cmap_obj)
 div_map.unif_sym_cmap(lift=None,
                        bitonic=False,
                        diffuse=True)
 fixed_cmap=div_map.get_mpl_color_map()

 for s in range(Hs.shape[0]):
  
  fig=plt.figure()
  ax1=fig.add_subplot(121)
  ax2=fig.add_subplot(122,projection='3d')
  
  
  h=Hs[s,0:dim,0:dim]
# Use bicubic interpolation with a linear normalization for color map  
  
  
  img=ax1.imshow(h, 
             cmap=fixed_cmap,
             #cmap=plt.get_cmap('gist_rainbow'), 
             interpolation='bicubic', 
             norm=mpl.colors.Normalize(vmin=a, vmax=b),
             extent=(0,141,141,0)
             )
  
  
  ax2.plot_surface(X,Y,h,cmap=fixed_cmap,rstride=1,cstride=1,norm=mpl.colors.Normalize(vmin=a, vmax=b))
  ax1.contour(h, levels, colors='black', linestyles='dashed', origin=None,linewidths=1.)
  ax1_divider=make_axes_locatable(ax1)
  cax1=ax1_divider.append_axes('right',size='5%', pad='1%')
  cb1=fig.colorbar(img, cax=cax1)
  cb1.ax.tick_params(labelsize=5)
  ax1.set_title("$V(k,k',\lambda=%s)$"%flowparams[s], fontsize=6)
  ax1.set_ylabel("$k'\,[\mathrm{fm}^{-1}]$",fontsize=6)
  ax1.set_xlabel("$k\,[\mathrm{fm}^{-1}]$", fontsize=6)
  ax2.tick_params(axis='x', labelsize=6, pad=0)
  ax2.tick_params(axis='y', labelsize=6, pad=0)
  ax2.tick_params(axis='z', labelsize=6, pad=0)
  ax1.tick_params(axis='x', labelsize=6, pad=1, length=3)
  ax1.tick_params(axis='y', labelsize=6, pad=1, length=3)
  ax1.set_xticklabels(['$0.0$','$1.0$','$2.0$','$3.0$','$4.0$','$5.0$','$6.0$','$7.0$'], fontsize=6)
  ax1.set_yticklabels(['$0.0$','$1.0$','$2.0$','$3.0$','$4.0$','$5.0$','$6.0$','$7.0$'], fontsize=6)
  ax2.set_xlabel("$k'\,[\mathrm{fm}^{-1}]$", fontsize=6)
  ax2.set_ylabel('$k\,[\mathrm{fm}^{-1}]$',fontsize=6)
  ax2.set_title("$V(k,k',\lambda=%s)$"%flowparams[s], fontsize=6)
  ax2.view_init(45)
  
# Save the plots in a sequence of pdf files and then merge them into a single one

  plt.savefig('/home/utente/SRG/results/simplified_potential/SRG_delta_flow_gif/%d.png'%s)
  path='/home/utente/SRG/results/simplified_potential/SRG_delta_flow/'
  pdf=PdfPages(path+str(s)+'_delta.pdf')
  pdf.savefig(fig, bbox_inches='tight',pad=0.0)
  pdf.close()
 
 pdfs=[]
 for i in range(Hs.shape[0]):
  pdfs.append('/home/utente/SRG/results/simplified_potential/SRG_delta_flow/%d_delta.pdf' %i)
                            
 merger=PdfFileMerger()
 for pdf in pdfs:
  merger.append(pdf)

 merger.write('/home/utente/SRG/results/simplified_potential/SRG_delta_flow/SRG_delta_flow.pdf')
 merger.close()

# Some lines to disable some warnings

 np.seterr(divide='ignore', invalid='ignore')
 mpl.rc('figure', max_open_warning = 0)

# Create the frames
 frames = []
 imgs = natsort.natsorted(glob.glob("/home/utente/SRG/results/simplified_potential/SRG_delta_flow_gif/*.png"))
 for i in imgs:
   new_frame = Image.open(i)
   frames.append(new_frame)
 
# Save into a GIF file that loops forever
 frames[0].save('/home/utente/SRG/results/simplified_potential/SRG_delta_flow_gif/SRG_delta.gif', format='GIF',
               append_images=frames[1:],
               save_all=True,
               duration=300, loop=0)
 
 return

#----------------------------------------------------------------------------------------------------------------------------
# Define a function to create and print 2D projections of chiral potential for SRG flow
#---------------------------------------------------------------------------------------------------------------------------- 

def plot_snapshots_chiral(Hs,flowparams, momenta, qMax):
  
  
  
  
# Define some indeces to print only a sector for each submatrix related to the channels

  dist=np.absolute(momenta-qMax)
  c=min(dist)
  semidist=int(len(dist)/2)
  for i in range(0,semidist):
   if dist[i]==c:
    indexmax=i
    
    break
  for j in range(semidist,int(len(dist))):
   if dist[j]==c:
    indexxmax=j
   
    break
  
  
  edge = int(len(momenta)/2)
  levels = np.arange(-1, 1, 0.12)

# Print for each value of the flow parameter a grid with 4 subplots each containing the submatrices up to a given maximum momentum

  for s in range(Hs.shape[0]):
   
   fig=plt.figure()
   grid=AxesGrid(fig, 111, 
                  nrows_ncols=(1,1),
                  axes_pad=0.,
                  share_all=True,
                  cbar_location="right",
                  cbar_mode="single",
                  cbar_size="7%",
                  cbar_pad=0.15,
                  )
   
   
   h = np.vstack((np.hstack((Hs[s,0:indexmax,0:indexmax], Hs[s,0:indexmax,edge:indexxmax])), 
                 np.hstack((Hs[s,edge:indexxmax,0:indexmax], Hs[s,edge:indexxmax,edge:indexxmax]))
                ))
 
   im = grid[0].imshow(h,
                     cmap=plt.get_cmap('Spectral_r'),                            
                     interpolation='bicubic',
                     norm=mpl.colors.Normalize(vmin=-0.5, vmax=+0.5),
                     )
             
   grid[0].contour(h, levels, colors='black', linestyles='dashed', origin=None,linewidths=0.5)
   grid[0].set_title("$V(k,k',\lambda=%s)$"%flowparams[s])
   grid[0].set_xticks([0,20,40,60,80,100,120,140,160])
   grid[0].set_yticks([0,20,40,60,80,100,120,140,160])
   grid[0].set_xticklabels(['$0.0$','$1.0$','$2.0$','$3.0$','$4.0$','$1.0$','$2.0$','$3.0$','$4.0$'])
   grid[0].set_yticklabels(['$0.0$','$1.0$','$2.0$','$3.0$','$4.0$','$1.0$','$2.0$','$3.0$','$4.0$'])
   grid[0].tick_params(axis='both',which='both',width=1,length=3)
   grid[0].xaxis.set_label_text("$k\,[\mathrm{fm}^{-1}]$")
   grid[0].yaxis.set_label_text("$k'\,[\mathrm{fm}^{-1}]$")
   grid[0].axvline(x=[80],ls='--',color='black',linewidth=0.8)
   grid[0].axhline(y=[80],ls='--', color='black',linewidth=0.8)
   cbar=grid.cbar_axes[0]
   plt.colorbar(im, cax=cbar)
   

# Print all the pdf files and then merge them together 
 
   plt.savefig('/home/utente/SRG/results/realistic_potential/SRG_chiral_flow/SRG_chiral_flowb_gif/%d.png'%s)
   path='/home/utente/SRG/results/realistic_potential/SRG_chiral_flow/SRG_chiral_flowb/'
   pdf=PdfPages(path+str(s)+'_chiral.pdf')
   pdf.savefig(fig, bbox_inches='tight')
   pdf.close()
   
  pdfs=[]
  
  for i in range(Hs.shape[0]):
    pdfs.append('/home/utente/SRG/results/realistic_potential/SRG_chiral_flow/SRG_chiral_flowb/%d_chiral.pdf' %i)
                            
  merger=PdfFileMerger()
  for pdf in pdfs:
   merger.append(pdf)

  merger.write('/home/utente/SRG/results/realistic_potential/SRG_chiral_flow/SRG_chiral_flowb/SRG_chiral_flow.pdf')
  merger.close()

  np.seterr(divide='ignore', invalid='ignore')
  mpl.rc('figure', max_open_warning = 0)

# Create the frames
  frames = []
  imgs = natsort.natsorted(glob.glob("/home/utente/SRG/results/realistic_potential/SRG_chiral_flow/SRG_chiral_flowb_gif/*.png"))
  for i in imgs:
    new_frame = Image.open(i)
    frames.append(new_frame)
 
# Save into a GIF file that loops forever
  frames[0].save('/home/utente/SRG/results/realistic_potential/SRG_chiral_flow/SRG_chiral_flowb_gif/SRG_chiral.gif', format='GIF',
               append_images=frames[1:],
               save_all=True,
               duration=300, loop=0)

                 
  return

#---------------------------------------------------------------------------------------------------------------------------
# Define a function to print 3D graphs for the SRG flow
#---------------------------------------------------------------------------------------------------------------------------

def plot_3d_chiral(Hs, flowparams, momenta, qMax):

 

 

# Find indeces to print the submatrices related to the channels up to a given value of the momentum

 dist=np.absolute(momenta-qMax)
 c=min(dist)
 semidist=int(len(dist)/2)
 for i in range(0,semidist):

  if dist[i]==c:
   indexmax=i
   
   break
 for j in range(semidist,int(len(dist))):
  if dist[j]==c:
   indexxmax=j
  
   
   break

 x=np.linspace(0.0,4.0,80)
 y=np.linspace(0.0,4.0,80)
 X,Y=np.meshgrid(x,y)
 edge=int(len(momenta)/2)
 levels = np.arange(-2, 1, 0.12)
 
 for s in range(Hs.shape[0]):
  
  fig=plt.figure()
  ax0=fig.add_subplot(221, projection='3d')
  ax1=fig.add_subplot(222, projection='3d')
  ax2=fig.add_subplot(223, projection='3d')
  ax3=fig.add_subplot(224, projection='3d')
  plt.tight_layout()

  h0 = Hs[s,0:indexmax,0:indexmax]
  h1= Hs[s,0:indexmax,edge:indexxmax]
  h2=Hs[s,edge:indexxmax,0:indexmax]
  h3= Hs[s,edge:indexxmax,edge:indexxmax]
                
# Draw the surfaces using the values of the submatrices as f(X,Y) bivariable function 
  
  surf1=ax0.plot_surface(X,Y,h0,cmap='Spectral_r',rstride=1,cstride=1,norm=mpl.colors.Normalize(vmin=-0.5, vmax=+0.5))
  surf2=ax1.plot_surface(X,Y,h1,cmap='Spectral_r',rstride=1,cstride=1,norm=mpl.colors.Normalize(vmin=-0.5, vmax=+0.5))
  surf3=ax2.plot_surface(X,Y,h2,cmap='Spectral_r',rstride=1,cstride=1,norm=mpl.colors.Normalize(vmin=-0.5, vmax=+0.5))
  surf4=ax3.plot_surface(X,Y,h3,cmap='Spectral_r',rstride=1,cstride=1,norm=mpl.colors.Normalize(vmin=-0.5, vmax=+0.5))
  cbar1=fig.colorbar(surf1, ax=ax0, shrink=0.5 )
  cbar2=fig.colorbar(surf2, ax=ax1, shrink=0.5)
  cbar3=fig.colorbar(surf3, ax=ax2, shrink=0.5) 
  cbar4=fig.colorbar(surf4, ax=ax3, shrink=0.5)
  cbar1.set_label("$^{3}S1$",size=8)
  cbar2.set_label("$^{3}S1-^{3}D1$",size=8)
  cbar3.set_label("$^{3}S1-^{3}D1$",size=8)
  cbar4.set_label("$^{3}D1$",size=8)
  cbar1.ax.tick_params(labelsize=5)
  cbar2.ax.tick_params(labelsize=5)
  cbar3.ax.tick_params(labelsize=5)
  cbar4.ax.tick_params(labelsize=5)
  ax0.tick_params(axis='x', labelsize=5)
  fig.suptitle("$V(k,k',\lambda=%s)$"%flowparams[s], fontsize=10)
  ax0.set_xlabel("$k\,[\mathrm{fm}^{-1}]$", fontsize=5)
  ax0.set_ylabel("$k'\,[\mathrm{fm}^{-1}]$", fontsize=5)
  ax1.set_xlabel("$k\,[\mathrm{fm}^{-1}]$", fontsize=5)
  ax1.set_ylabel("$k'\,[\mathrm{fm}^{-1}]$", fontsize=5)
  ax2.set_xlabel("$k\,[\mathrm{fm}^{-1}]$", fontsize=5)
  ax2.set_ylabel("$k'\,[\mathrm{fm}^{-1}]$", fontsize=5)
  ax3.set_xlabel("$k\,[\mathrm{fm}^{-1}]$", fontsize=5)
  ax3.set_ylabel("$k'\,[\mathrm{fm}^{-1}]$", fontsize=5)
  ax0.tick_params(axis='x', labelsize=5, pad=-3)
  ax0.tick_params(axis='y', labelsize=5, pad=-3)
  ax0.tick_params(axis='z', labelsize=5, pad=-3)
  ax1.tick_params(axis='x', labelsize=5, pad=-3)
  ax1.tick_params(axis='y', labelsize=5, pad=-3)
  ax1.tick_params(axis='z', labelsize=5, pad=-3)
  ax2.tick_params(axis='x', labelsize=5, pad=-3)
  ax2.tick_params(axis='y', labelsize=5, pad=-3)
  ax2.tick_params(axis='z', labelsize=5, pad=-3)
  ax3.tick_params(axis='x', labelsize=5, pad=-3)
  ax3.tick_params(axis='y', labelsize=5, pad=-3)
  ax3.tick_params(axis='z', labelsize=5, pad=-3)
  
# Print and save all the pdfs and then merge them together 
  
  plt.savefig('/home/utente/SRG/results/realistic_potential/SRG_chiral_flow/SRG_chiral3d_flow_gif/%d.png'%s)
  path='/home/utente/SRG/results/realistic_potential/SRG_chiral_flow/SRG_chiral3d_flow/'
  pdf=PdfPages(path+str(s)+'_chiral3d.pdf')
  pdf.savefig(fig, bbox_inches='tight')
  pdf.close()

 pdfs=[]
 for i in range(Hs.shape[0]):
  pdfs.append('/home/utente/SRG/results/realistic_potential/SRG_chiral_flow/SRG_chiral3d_flow/%d_chiral3d.pdf' %i)
                            
 merger=PdfFileMerger()
 for pdf in pdfs:
  merger.append(pdf)

 merger.write('/home/utente/SRG/results/realistic_potential/SRG_chiral_flow/SRG_chiral3d_flow/SRG_chiral3d_flow.pdf')
 merger.close()

# Some lines to disable the warnings

 np.seterr(divide='ignore', invalid='ignore')
 mpl.rc('figure', max_open_warning = 0)

# Create the frames
 frames = []
 imgs = natsort.natsorted(glob.glob("/home/utente/SRG/results/realistic_potential/SRG_chiral_flow/SRG_chiral3d_flow_gif/*.png"))
 for i in imgs:
   new_frame = Image.open(i)
   frames.append(new_frame)
 
# Save into a GIF file that loops forever
 frames[0].save('/home/utente/SRG/results/realistic_potential/SRG_chiral_flow/SRG_chiral3d_flow_gif/SRG_chiral.gif', format='GIF',
               append_images=frames[1:],
               save_all=True,
               duration=300, loop=0)
 
 return



