#include "/home/queenmab/GitHub/Math/matrixes.h"
#include "DMLimitsLib.h"
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

void Run()
{
  int Ndif = 6;
  int Nps = 10;
  int Nbar = Ndif+Nps;
  
  Init(20,20,20,Ndif,Nps);
  
  TString extended[Ndif] = {"Irf",
                            "Leptonic",
                            "Hadronic",
                            "3FHL_J0500.9-6945e",
                            "3FHL_J0530.0-6900e",
                            "3FHL_J0531.8-6639e"};
  
  TString point[Nps] = {"J0537-691",
                        "J0524.5-6937",
                        "J0534.1-6732",
                        "J0525.2-6614",
                        "J0535.3-6559",
                        "J0454.6-6825",
                        "J0537.0-7113",
                        "J0535-691",
                        "J0525-696",
                        "J0509.9-6418"};
  TString suf = "_KSPpointing_v2_";
  TString suf_DM = "_jfactorNFW";
  
  FillContainer_Bkg(extended,point,suf);
  FillContainer_DM(1,"W",suf_DM);
  FillContainer_Obs("Irf+CR+DiffuseSources+PS",true,suf);
  //Ntotal = DataSim(Obs_data);

  Number steps[Nbar+1]={10,0.001,1,1,1,1,0.5,0.001,0.5,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001};
  
  V Kpars; init(Kpars,Nbar+1);
  
  cout << "Maximizing Likelihood..." << endl;
  cout << endl;
  cout << calc_MaxlogL(Kpars,steps,false) << endl;

  cout << "Maximum Likelihood parameters: " << endl;
  for (int ii=0; ii<Nbar+1; ii++) cout << Kpars[ii] << "  ";
  cout << endl;
  cout << endl;
  cout << "Calculating Upper Limit on DM normalization..." << endl;
  Number UpperLimit = Upper_Minimizer(Kpars,false);
  cout << "Upper Limit: " << UpperLimit << endl;
  cout << endl;
  Number intervals[Nbar+1] = {UpperLimit,0.00025,0.25,0.5,0.02,0.06,0.015,0.1,0.15,4.5,1,4,2,1,0.1,0.1,200};
  V Cfactors;
  cout << "Calculating Correlation Factors..." << endl;
  calc_CorrFactors(Kpars,intervals,Cfactors);
}


void Run2()
{
  Init(20,20,20,1,0);
  
  TString extended[N_ext] = {"Irf"};
  TString point[N_ps];
  TString suf = "_KSPpointing_v2_";
  TString suf_DM = "_jfactorgamma1.5";
  
  FillContainer_Bkg(extended,point,suf);
  FillContainer_DM(1,"W",suf_DM);
  FillContainer_Obs("Irf",true,suf);
  //Ntotal = DataSim(Obs_data);
  Number steps[Nbar+1]={10,0.001};
  
  V Kpars; init(Kpars,Nbar+1);
  
  cout << "Maximizing Likelihood..." << endl;
  cout << endl;
  calc_MaxlogL(Kpars,steps,false);

  cout << "Maximum Likelihood parameters: " << endl;
  for (int ii=0; ii<Nbar+1; ii++) cout << Kpars[ii] << "  ";
  cout << endl;
  cout << endl;
  cout << "Calculating Upper Limit on DM normalization..." << endl;
  Number UpperLimit = Upper_Minimizer(Kpars,false);
  cout << "Upper Limit: " << UpperLimit << endl;
  cout << endl;
  Number intervals[Nbar+1] = {UpperLimit,0.00025};
  V Cfactors;
  cout << "Calculating Correlation Factors..." << endl;
  calc_CorrFactors(Kpars,intervals,Cfactors);
}

void CheckJFACTOR()
{
  Init(20,20,20,2,0);
  
  TString extended[N_ext] = {"Irf","M31disk"};
  TString point[N_ps];
  TString suf = "_";
  TString suf_DM = "_jfactorNFW";

  TString DMtypes[1] = {"_jfactorNFW"};
  Number DMmasses[8] = {0.100,0.200,0.500,1,5,10,50,100};
  TString DMparticles[4] = {"W"};

  Number steps[Nbar+1]={10,0.001,0.1};
  
  FillContainer_Bkg(extended,point,suf);
  FillContainer_Obs("Irf+M31point",true,suf);

  for (int npart=0; npart<1; npart++)
    {
      for (int ntype=0; ntype<1; ntype++)
	{
	  TString filenamelim = "/home/queenmab/GitHub/M31pipe/results/Limits"+DMparticles[npart]+DMtypes[ntype]+"_Disk";
	  TString filenamecfact = "/home/queenmab/GitHub/M31pipe/results/Cfactors"+DMparticles[npart]+DMtypes[ntype]+"_Disk";
	  ofstream outfilelim;
	  ofstream outfilecfact;
	  outfilelim.open(filenamelim,ios::app);
	  outfilecfact.open(filenamecfact,ios::app);
	  
	  for (int nmass=0; nmass<8; nmass++)
	    {
	      FillContainer_DM(DMmasses[nmass],DMparticles[npart],DMtypes[ntype]);
	      V Kpars; init(Kpars,Nbar+1);
	      
	      cout << "Maximizing Likelihood..." << endl;
	      cout << endl;
	      cout << calc_MaxlogL(Kpars,steps,false) << endl;
	      
	      cout << "Maximum Likelihood parameters: " << endl;
	      for (int ii=0; ii<Nbar+1; ii++) cout << Kpars[ii] << "  ";
	      cout << endl;
	      cout << endl;
	      cout << "Calculating Upper Limit on DM normalization..." << endl;
	      Number UpperLimit = Upper_Minimizer(Kpars,false);
	      cout << "Upper Limit: " << UpperLimit << endl;
	      cout << endl;
	      outfilelim << DMmasses[nmass] << "  " << UpperLimit << endl;
	      Number intervals[Nbar+1] = {UpperLimit,0.00025,3};
	      V Cfactors;
	      cout << "Calculating Correlation Factors..." << endl;
	      calc_CorrFactors(Kpars,intervals,Cfactors);
	      outfilecfact << DMmasses[nmass] << "  ";
	      for (int ii=0; ii<Nbar; ii++) outfilecfact << Cfactors[ii] << "  ";
	      outfilecfact << endl;
	    }
	}
    }
}
