import zfit
import uproot
import numpy as np
import matplotlib.pyplot as plt
from zfit import z

import os
import sys
sys.path.append("/home/oscar/Documents/Thesis/Code/scripts")
import plot_tools

#Importing MC
file = uproot.open("/home/oscar/Documents/Thesis/RootFiles/Prueba2/ntuple_BuJpsiK_ForXGboost_Year2022_Sample3_v0.root")
tree = file["treeBu;1"]
data = tree["massB"].array(library="pd")
data_fit = data.sample(n=70000, random_state=20)
datanp = data_fit.to_numpy()

#starting the Fit
space = zfit.Space("x", [4.8,5.5])

param1G = zfit.Parameter("mediano1",5.2793)
param2G = zfit.Parameter("scale1",0.0394556)

param2DCB = zfit.Parameter("sigma",0.0279811)
param3DCB = zfit.Parameter("alphal",1.07502)
param4DCB = zfit.Parameter("nl",2.26597)
param5DCB = zfit.Parameter("alphar",1.04001)
param6DCB = zfit.Parameter("nr",7.13824)
fraction = zfit.Parameter("frac",0.4)

dcb = zfit.pdf.DoubleCB(param1G, param2DCB, param3DCB, param4DCB, param5DCB, param6DCB, space)
gauss = zfit.pdf.Gauss(param1G, param2G, space)

modelo = zfit.pdf.SumPDF([gauss, dcb], fraction)

dataZ = zfit.Data.from_numpy(obs=space, array=datanp)
nll = zfit.loss.UnbinnedNLL(modelo, dataZ)
Minuit = zfit.minimize.Minuit()
minimum = Minuit.minimize(nll)
print(minimum)
print(minimum.hesse())

print(modelo.get_params())


#three_sigma1 = 0.4*0.0465675 + 0.6*0.0280631
#print(modelo.integrate([5.27905 - three_sigma1,5.27905 + three_sigma1]).numpy()[0])

#print(three_sigma1)
#for i in range (10000):
#    value = modelo.integrate([5.27905 - three_sigma1 -(i*0.0001),5.27905 + three_sigma1+(i*0.0001)])
#    if (value.numpy()[0]>0.6827):
#        three_sigma2 = three_sigma1 +  i*0.0001
#        print(three_sigma2)
#        break  
#SBL2 = 5.27905-(3*three_sigma2)
#SBR2 = 5.27905+(3*three_sigma2)
#print("Left Side Band: ", SBL2)
#print("Right Side BandL: ",  SBR2)

figure = plt.figure()
axes = plot_tools.create_axes_for_pulls(figure)
plot_tools.plot_model(datanp, modelo,
                      bins=50,
                      axis=axes[0],
                      pulls=True,
                      level=2,
                      plot_components=True,
                      axis_pulls=axes[1], chi_x=0.05, chi_y=0.9,
                      SBL=0, SBR=0, SBR2=0, SBL2=0)
plt.savefig('/home/oscar/Documents/Thesis/Code/Fits/Imagenes/massB_FitGood_NonRes_Antiradveto_Check1.png' , format='png')
plt.close()