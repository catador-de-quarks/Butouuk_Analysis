import zfit
import uproot
import numpy as np
import matplotlib.pyplot as plt
from zfit import z
import PDFS

import os
import sys
sys.path.append("/home/oscar/Documents/Thesis/Code/scripts")
import plot_tools


#Importing MC
file = uproot.open("/home/oscar/Documents/Thesis/RootFiles/ntuple_BuJpsik_MiniAOD_Year2022_Sample1_Preselection_MC.root")
tree = file["treeBu;1"]
data = tree["massB"].array(library="pd")
data_fit = data.sample(n=70000)
datanp = data_fit.to_numpy()

#starting the Fit
space = zfit.Space("x", [4.8,5.5])

param1G = zfit.Parameter("mediano1",5.26836,5.2,5.3)
param2G = zfit.Parameter("scale1",0.0426488,0.04,0.05)
param1DCB = zfit.Parameter("mediano2",5.27651,5.2,5.3)
param2DCB = zfit.Parameter("sigma",0.0260087,0.01,0.05)
param3DCB = zfit.Parameter("alphal",1.29111,1.0,1.6)
param4DCB = zfit.Parameter("nl",1.0,0.9,1.6)
param5DCB = zfit.Parameter("alphar",2.0,1.5,2.2)
param6DCB = zfit.Parameter("nr",0.9,0.9,1.5)
lambda1 = zfit.Parameter("lambda1", -0.8,-1.2,-0.5)

"""
param1G = zfit.param.ConstantParameter("mediano1",5.26836)
param2G = zfit.param.ConstantParameter("scale1",0.0426488)
param1DCB = zfit.param.ConstantParameter("mediano2",5.27651)
param2DCB = zfit.param.ConstantParameter("sigma",0.0260087)
param3DCB = zfit.param.ConstantParameter("alphal",0.861356)
param4DCB = zfit.param.ConstantParameter("nl",5.90993)
param5DCB = zfit.param.ConstantParameter("alphar",8.69809)
param6DCB = zfit.param.ConstantParameter("nr",0.434962)
lambda1 = zfit.param.ConstantParameter("lambda1", 0.0685827)
"""

dcb = zfit.pdf.DoubleCB(param1DCB, param2DCB, param3DCB, param4DCB, param5DCB, param6DCB, space)

log = PDFS.Logistic(param1G,param2G, space)
exp = zfit.pdf.Exponential(lambda1, space)
gauss = zfit.pdf.Gauss(param1G, param2G, space)

modelo = zfit.pdf.SumPDF([exp,gauss,dcb],[0.1, 0.4])

dataZ = zfit.Data.from_numpy(obs=space, array=datanp)
nll = zfit.loss.UnbinnedNLL(modelo, dataZ)
Minuit = zfit.minimize.Minuit()
minimum = Minuit.minimize(nll)
print(minimum)
print(minimum.hesse())




figure = plt.figure()
axes = plot_tools.create_axes_for_pulls(figure)
plot_tools.plot_model(datanp, modelo,
                      bins=50,
                      axis=axes[0],
                      pulls=True,
                      level=2,
                      plot_components=True,
                      axis_pulls=axes[1], chi_x=0.05, chi_y=0.9
                      )
plt.savefig('/home/oscar/Documents/Thesis/Code/Fits/Imagenes/massB_FitGood.png' , format='png')
plt.close()