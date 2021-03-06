
# This program simulates, bins and creates model cubes for several Dark Matter models using a Pointing pattern similar to the one propopsed in the CTA KSP for LMC.

import gammalib
import ctools
import cscripts
import numpy as np
import multiprocessing
import math
import shutil
import os
import random
import config

PATH_MODEL = "../../models/" #Path where models are stored
PATH_OBS = config.DATA_PATH+"/Obs_DM/" #Path to store resulting .fits files.
PATH_HERE = "../pipes_in_py/" #Path where we are running

# If the observation path doesn't exist, create it.
if not os.path.exists(PATH_OBS): 
    os.makedirs(PATH_OBS)

#Pointing
centerx = 10.76
centery = 41.19

r = 3.0
ra_list=np.zeros(7)
dec_list=np.zeros(7)
ra_list[0] = centerx;
dec_list[0] = centery;
for i in range(1,7):
    angle=(i-1)*2*math.pi/6
    ra_list[i] = r*math.cos(angle)+centerx
    dec_list[i] = r*math.sin(angle)+centery


#Observation variables

rad = 5 #Radius of ROI
emin = 0.03 #Minimum energy in TeV
emax = 100.0 #Maximum enery in TeV
tstart = 0.0 #Startind time
duration = 180000 #Ending time
deadc = 0.95 #Dead time

binsz = 0.5 #Spatial binning
nxpix = 20
nypix = 20

enumbins = 20

caldb = gammalib.GCaldb(config.CTOOLS_PATH+"/share/caldb/data/cta/1dc/bcf/North_z20_50h") #Calibration Files for gammalib class
irf = "irf_file.fits"

caldb_= "1dc"
irf_="North_z20_50h"

particle = 'Tau' #Final state particle
#masses = [0]
#masses = [0.100,0.200,0.500,1,5,10,50,100]
masses = [0.100,0.200,0.300,0.400,0.500,0.600,0.800,1,4,5,8,10,40,50,80,100]


for mass in masses:
    rndseed = random.randint(1,300000)
    #Build wisely the strings for the filenames
    #______________________________________________
    if mass < 1:
        masstr = str(int(mass*1000))+'GeV'
    else:
        masstr = str(mass)+'TeV'
    specname = 'flux'+particle+masstr+'.txt'
    modelname = 'dm_M31_'+particle+masstr
    suf = '_jfactorNFW'
    
    model = PATH_MODEL+modelname+suf+'.xml'
    time = str(int(duration/3600))
    outfile = PATH_OBS+'observations_'+'LMC_'+modelname+suf+'.xml' #List of Observations file that will be produced ('.xml')
    cntcube = PATH_OBS+"cntcube_"+modelname+suf+'.fits'
    modcube = PATH_OBS+"modcube_"+modelname+suf+'.fits'
    #_______________________________________________
    
    file = open(outfile,'w')
    Obs_list = gammalib.GObservations()
    xml = gammalib.GXml(outfile)
    
    # RUN SIMULATION
    number=1
    for ra in ra_list:
        dec = dec_list[number-1]
        filename = 'events_'+'M31_'+modelname+'_KSP'+'0'+str(number)+'.fits'
        eventfile = gammalib.GFilename(filename);

        sim = ctools.ctobssim()
        sim["inmodel"]=model
        sim["seed"] = number+rndseed
        sim["outevents"]=filename
        sim["caldb"]=caldb_
        sim["irf"]=irf_
        sim["ra"]=ra
        sim["dec"]=dec
        sim["rad"]=rad
        sim["tmin"]=gammalib.GTime(tstart)
        sim["tmax"]=gammalib.GTime(tstart+duration)
        sim["emin"]=emin
        sim["emax"]=emax
        sim["debug"]=True
        sim.execute()

        #STORE THE OBSERVATION IN GObservations Class and store it in the .xml output file.

        #Allocate CTA observation
        obs = gammalib.GCTAObservation()
        #Set pointing direction
        pntdir = gammalib.GSkyDir()
        pntdir.radec_deg(ra,dec)
        
        pnt = gammalib.GCTAPointing()
        pnt.dir(pntdir)
        obs.pointing(pnt)
    
        #Set ROI

        roi = gammalib.GCTARoi()
        instdir = gammalib.GCTAInstDir()
        instdir.dir(pntdir)
        roi.centre(instdir)
        roi.radius(rad)
        
        #Set GTI
    
        gti = gammalib.GGti()
        start = gammalib.GTime(tstart)
        stop = gammalib.GTime(tstart+duration)
        gti.append(start,stop)
        
        #Set Energy Boundaries
        ebounds = gammalib.GEbounds()
        e_min = gammalib.GEnergy(emin,'TeV')
        e_max = gammalib.GEnergy(emax,'TeV')
        ebounds.append(e_min,e_max)
        
        #Allocate event list
        events = gammalib.GCTAEventList(eventfile)
        obs.eventfile(eventfile)
        events.roi(roi)
        events.gti(gti)
        events.ebounds(ebounds)
        obs.events(events)

        #Set instrument response
        obs.response(irf,caldb)
        
        #Set ontime, livetime, and deadtime correction factor
        obs.ontime(duration)
        obs.livetime(duration*deadc)
        obs.deadc(deadc)
        
        obs.id(str(number))
        obs.name('events_'+'M31_'+modelname+'_KSP'+'0'+str(number))

        Obs_list.append(obs)    
        number=number+1

    Obs_list.models(model)
    Obs_list.save(outfile)

    #Move observations to destination folder.Observations are simulated and stored in the present folder, but I prefer to store them somewhere else, so I move them there.

    for src_dir,dirs,files in os.walk(PATH_HERE):
        dst_dir = PATH_OBS
        if not os.path.exists(dst_dir):
            os.makedirs(dst_dir)
        for file in files:
            src_file = os.path.join(PATH_HERE,file)
            dst_file = os.path.join(PATH_OBS,file)
            if file.endswith('.fits'):
                if os.path.exists(dst_file):
                    os.remove(dst_file)
                shutil.move(src_file,PATH_OBS)


    
    # BIN THE DATA
    binn = ctools.ctbin()
    binn["inobs"]=outfile
    binn["outcube"]=cntcube
    binn["coordsys"] = "CEL"
    binn["proj"] = "CAR"
    binn["ebinalg"] = "LOG"
    binn["xref"]=centerx
    binn["yref"]=centery
    binn["nxpix"] = nxpix
    binn["nypix"] = nypix
    binn["binsz"] = binsz
    binn["enumbins"] = enumbins
    binn["emin"]=emin
    binn["emax"]=emax
    binn["debug"]=True
    binn.execute()

    #PRODUCE MODELCUBE
    mod = ctools.ctmodel()
    mod["inobs"]= outfile
    mod["inmodel"]= model
    mod["outcube"] = modcube
    mod["incube"] = cntcube
    mod["expcube"] = "NONE" 
    mod["psfcube"] = "NONE"
    mod["bkgcube"] = "NONE"
    mod["caldb"]=caldb_
    mod["irf"]=irf_
    mod["debug"]=True
    mod.execute()
   
