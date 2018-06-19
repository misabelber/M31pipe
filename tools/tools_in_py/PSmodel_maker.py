
#This program produces separated xml models for the Point Sources contained in the general model 
#file with all the components 
import gammalib
import math
import numpy as np

#PATHS

LMC_PATH = '/home/queenmab/GitHub/LMC'
PATH_MODELS=LMC_PATH+'/models/'  #Path to store the models created (.xml)

#Open model file
models = gammalib.GModels(PATH_MODELS+'LMC_closer_files.xml') #This modelfile contains the models
# of all the point sources in the LMC region, but the ones that are inside the 10X10 ROI are  
#"marked" (their parameteres are set free, while the ones outside this ROI have their 
#parameters fixed)

for model in models:
    modelcontainer = gammalib.GModels()
    modelname = model.name()
    srctype = model.type()
    par = model.spectral().at(0)
    if srctype=='PointSource' and par.is_free(): #Do only for PS with free parameters (inside ROI)
        modelcontainer.append(model)
        filename = 'LMC_'+modelname+'.xml'
        modelcontainer.save(PATH_MODELS+filename)
    
