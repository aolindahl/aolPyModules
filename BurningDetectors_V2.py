# -*- coding: utf-8 -*-
"""
Created on Mon Oct 13 17:50:18 2014

@author: hartmang
"""
import numpy as np
#detectors' angles
phi=[i*22.5 for i in range(0,16)]
phir=[np.radians(i) for i in phi]
#basis vectors
Bt=[[1.,np.cos(2*phir[i]),np.sin(2*phir[i])] for i in range(0,16)]
Bs=[[1.,np.cos(2*phir[i])**2,np.sin(2*phir[i])**2] for i in range(0,16)]
B0=np.sum(Bs,axis=0)
B1=np.sqrt(B0)
B=Bt/B1
B=np.array([[1./16., np.cos(2*phir[i])/8., np.sin(2*phir[i])/8. ] for i in
    range(0,16)])

# angular distribution
# I set the random numbers to fixed parameters in order to compare it with the Matlab Script -fine!
a=np.random.normal(size=1)
b=1+np.random.normal(size=1)
I0 = (100+200*np.random.normal(size=1))
d=np.radians(360*np.random.normal(size=1))

a=0.5
b=2.
I0=100.
d=np.radians(22.5)

Itest0=[I0 * (1. + b/4. * (3. * a * np.cos( 2.*(i-d) ) + 1.)) for i in phir]
Inoise=[min(0,np.sqrt(Itest0[i])*np.random.normal(size=1)) for i in range(0,16)]
# check: Inoise=[0 for i in range(0,16)] -fine!

Itest=[Itest0[i]+float(Inoise[i]) for i in range(0,16)]

c=np.dot(B.transpose(),Itest0)

#operating status
hassig0=[1,0,1,1,1,1,1,1,1,0,1,0,0,0,0,1]

#matrix stuff
hassig0b=[0.5 for i in range(0,16)]
for i in range(0,16):
    if hassig0[i]==1:
        hassig0b[i]=0
    else:
        hassig0b[i]=1
Bn1=np.array([np.array([hassig0b[i]*B.transpose()[j][i] for i in range(0,16)]) for j in range(0,3)])
Bn2=Bn1.transpose()
Bfir=np.linalg.inv(np.eye(3)-np.dot(Bn1,Bn2))
Bsec=[sum([hassig0[i]*B.transpose()[j][i]*Itest0[i] for i in range(0,16)]) for j in range(0,3)]

cgab=np.dot(Bfir,Bsec)

print c
print cgab




