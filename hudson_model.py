from math import pi
from collections import namedtuple
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from numpy import dot, transpose,linalg, sort, cos,sin,sqrt

#Hudsons model for cracked media with "weak" inclusions is based on an assumed distribution of penny shaped
#ellipsoid cracks. The model is shown in detail in Mavko 1999
#NOTE:This module is modified from Abers 2016 TI model and MS_phasevelo functions designed for MatLAB
class TIcrackmodel():

    def __init__(self):
        self.name = 'model'

    #results global variables
    def loadvar(self,aspar,Kfl,Rfl,K,G,rho,Cagg,phiar,fdist,increment,data):
        self.fdist = fdist*pi/180.                                      #Converts degrees to radians
        self.aspar = 1./aspar
        self.phiar = phiar/1e2
        self.Kfl = Kfl
        self.Rfl = Rfl
        self.Cagg = Cagg
        self.Kiso = K
        self.Giso = G
        self.Ragg = rho
        self.increment = increment
        self.data = data
        #get ranges for loops
        self.n = aspar.size
        self.m = phiar.size

    #Calculates C matrix based on K,G,rho,phi and crack distribution.
    #Cagg used ONLY for anisotropic solid case. Here we use K,rho,G, and phi
    #to build transversely isotropic solid with cracks and porosity 
    def computeC(self,cd,asp,Kfl,Rfl,K,G,rho,phi,s):
        #lame parameters and hudson constants for weak inclusions
        lam = K-(2/3.)*G
        mu = G
        #Hudson constants for "weak" inclusions (p. 195-6 in Mavko) within ISOTROPIC solid
        kapa = Kfl*(lam+2.*mu)/(pi*asp*mu*(lam+mu))
        u3 = (4/3.)*(lam+2*mu)/((lam+mu)*(1+kapa))
        u1 = (16/3.)*(lam+2*mu)/(3*lam+4*mu)
        #Fischer stat terms:  here "s" is the std. dev. of angle (radian): Fischer term
        #is exp(cos(self.inc)/s^2)  where theta is angle from symmetry (3) axis. For small-ish theta, this is exp(-theta^2/2/s^2)
        #in Abers code s = fdist in calculation
        ex = np.exp(1/(s**2))
        e11 = (-1+2.*s**2.*ex-2*s**4.*(ex-1))/(2*(ex-1))
        e1111 = (3/8.)*(-1+4*s**4.*(2*ex+1)-24*s**6.*ex+24*s**8.*(ex-1))/(ex-1)
        #Terms for rotationally symmetric  cracks around theta=0.  SHOULD ONLY
        #NEED:  e11, e22, e33,  e1111, e1122, e1133, e3333, e2323 (to inforce TI)
        e22 = e11
        e33 = 1-2*e11
        e1122 = e1111/3.
        e1133 = e11-4/3.*e1111
        e3333 = 8/3.*e1111-4*e11+1
        e2323 = e1133
        #Build the minimum number of cijkl for TI anisotropy
        c1111 = 4*cd*mu*u1*(e1111-e11)-cd/mu*u3*(lam**2+4*lam*mu*e11+4*mu**2.*e1111)
        c1122 = 4*cd*mu*u1*e1122-cd/mu*u3*(lam**2+2.*lam*mu*(e22+e11)+4*mu**2.*e1122)
        c1133 = 4*cd*mu*u1*e1133-cd/mu*u3*(lam**2+2.*lam*mu*(e33+e11)+4*mu**2.*e1133)
        c3333 = 4*cd*mu*u1*(e3333-e33)-cd/mu*u3*(lam**2+4*lam*mu*e33 +4*mu**2.*e3333)
        c2323 = cd*mu*u1*(4*e2323-e22-e33)-4*cd*mu*u3*e2323   #Fixed
        #Add C1 to C0 (the isotropic original)
        #C=zeros(6);  Voigt Notation
        C = np.zeros([6,6])
        C[0,0] = lam+2*mu+c1111
        C[0,1] = lam+c1122
        C[0,2] = lam+c1133
        C[2,2] = lam+2*mu+c3333
        C[3,3] = mu+c2323
        #These hold for TI:
        C[1,1] = C[0,0]
        C[1,2] = C[0,2]
        C[4,4] = C[3,3]
        C[5,5] = (C[0,0]-C[0,1])/2.
        #Symmetry terms
        C[1,0] = C[0,1]
        C[2,0] = C[0,2]
        C[2,1] = C[1,2]
        #crack density and cracked rock density
        ####Crack porosity set equal to rock porosity itself...
        ####may need to adjust this as needed
        #phicrack = 4*pi*asp/(3.*cd)
        #compute density
        self.rho = (1-phi)*rho+phi*Rfl
        #create TI matrix
        self.Cti = C
        return self.Cti,self.rho

    #Create iterable cordinate transformation
    def computeROT(self):
        self.npts = np.int(90./self.increment)
        #create incidence angle array
        self.inc = np.linspace(0,90,self.npts)
        #create azimuth array
        self.azm = np.zeros([self.npts,])
        self.l = self.azm.size
        #convert from degrees to radians
        self.azmr = self.azm*(pi/180.0)
        self.incr = self.inc*(pi/180.0)
        caz = cos(self.azmr)
        saz = sin(self.azmr)
        cinc = cos(self.incr)
        sinc = sin(self.incr)
        #create X and r matricies
        self.X = np.zeros([self.l,3])
        self.r = np.zeros([self.l,])

        #Calculate X and r
        for i in range(0,self.l):
            self.X[i,0], self.X[i,1], self.X[i,2] = caz[i]*cinc[i], -saz[i]*cinc[i], sinc[i]
            #c normalise to direction cosines
            self.r[i] =(np.sqrt(self.X[i,0]*self.X[i,0]+self.X[i,1]*self.X[i,1]+self.X[i,2]*self.X[i,2]))
            self.X[i,0],self.X[i,1],self.X[i,2] =  self.X[i,0]/self.r[i],self.X[i,1]/self.r[i],self.X[i,2]/self.r[i]

    def phasevelo(self,X,C,rho,i):
        #Create gamma matrix
        gamma = np.array([[self.X[i,0], 0.0, 0.0, 0.0, self.X[i,2],self.X[i,1]],
                        [0.0, self.X[i,1], 0.0, self.X[i,2],0.0, self.X[i,0]],
                        [0.0, 0.0, self.X[i,2], self.X[i,1], self.X[i,0], 0.0 ]])
        #compute T matrix (3X3)
        T = dot(dot(gamma,C),transpose(gamma))
        #Calculate eigen values and eigenvectors
        [Eval,Evec] = linalg.eig(T)
        #divide by density
        V_r = sqrt([Eval/rho])
        return T,gamma,Eval,Evec,V_r

    def run_model(self):
        self.vel_agg = np.zeros([self.l,3])
        self.vel_mod = np.zeros([self.l,3,self.n])
        #compute anisotropic solid case
        for i in range(0,self.l):
            self.vel_agg[i,:] = self.phasevelo(self.X,self.Cagg*1e9,self.Ragg,i)[4]/1e3
        self.vel_agg = sort(self.vel_agg,axis=1)
        #Begin loop to compute velocities for each porosity and aspect ratio combination
        for j in range(0,self.m):
            #Crack density used in computing Cijkl elements
            cd = (3*self.phiar[j])/(4*pi*self.aspar[j]);
            #build TI stiffness matrix with Kiso, Giso, and agg density
            self.Cti,self.rho = self.computeC(cd,self.aspar[j],self.Kfl,self.Rfl,
                                            self.Kiso,self.Giso,self.Ragg,self.phiar[j],self.fdist)
            #Loop through all inclinations and save to a single array
            for i in range(0,self.inc.size):
                #compute phase velocities for each compbination, place into 3-D array
                self.vel_mod[i,:,j] = self.phasevelo(self.X,self.Cti,self.rho,i)[4]/1e3
            #Sort eigenvalues, smallest in column 0, largest (vp) in column 2
            self.vel_mod[:,:,j] = sort(self.vel_mod[:,:,j],axis=1)

    def plotfigure(self):
        #create array of grey values based on number of samples.
        color = np.linspace(0.2,0.8,self.m)
        #set up figure
        fig = plt.figure(figsize=(8,10))
        #create first subplot
        ax = fig.add_subplot(2,1,1)
        #create ghost line for legend headers
        ax.plot(0,0,label='   $\\phi$       $a$',c = '1')
        #loop over all samples
        for i in range(0,self.m):
            ax.plot(90-self.inc,self.vel_mod[:,2,i],label = '%.2f %.3f'%(self.phiar[i]*100,
                                                            self.aspar[i]),linewidth = 1.5,c = '%f' %color[i])
        #plot anisotropic solid case
        ax.plot(90-self.inc,self.vel_agg[:,2],label = 'Aniso solid', linewidth = 2.5,c = 'g')
        #set plot parameters and add peak stress cordinate
        ax.add_patch(patches.Rectangle((85, self.data[0]/1e3),5,.37,alpha = 0.5,edgecolor = 'NONE'))
        ax.add_patch(patches.Rectangle((0, self.data[1]/1e3),5,.7,alpha = 0.5,edgecolor = 'NONE'))
        ax.set_ylim([6.0,7])
        ax.set_ylabel('Vp [$km/s$]',fontsize = 18)
        ax.tick_params(axis='both', which='major', labelsize=18)
        plt.legend(loc='lower right', fontsize = 14)

        #Create second subplot
        ax = fig.add_subplot(2,1,2)
        #plot anisotropic solid example
        ax.plot(90-self.inc,self.vel_agg[:,1],label = '$Vs_1$', linewidth = 2.5,c = 'g')
        ax.plot(90-self.inc,self.vel_agg[:,0],label = '$Vs_2$', linewidth = 2.5,c = 'g',linestyle = '--')
        #loop over all samples
        for i in range(0,self.m):
            ax.plot(90-self.inc,self.vel_mod[:,1,i], linewidth = 1.5,c = '%f' %color[i])
            ax.plot(90-self.inc,self.vel_mod[:,0,i], linewidth = 1.5,c = '%f' %color[i],linestyle = '--')
        #set plot parameters and add data
        ax.add_patch(patches.Rectangle((85, self.data[2]/1e3),5,.4,alpha = 0.5,edgecolor = 'NONE'))
        ax.add_patch(patches.Rectangle((0, self.data[3]/1e3),5,.7,alpha = 0.5,edgecolor = 'NONE'))
        ax.set_ylim([1.5,4.5])
        ax.set_xlabel('Incidence [$\\theta$]', fontsize=18)
        ax.set_ylabel('Vs [$km/s$]',fontsize = 18)
        ax.tick_params(axis='both', which='major', labelsize=18)
        plt.legend(loc='lower right', fontsize = 14)
        plt.show()
      
    #copyright - Peter Miller 8/11/2021
