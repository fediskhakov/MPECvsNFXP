% This is the main script to create all the data for the
% Iskhakov,Rust,Schjerning comment for Judd and Su 2012 Econometrica

%add path to the functions
currdir=pwd;
addpath([currdir filesep 'matlab_include']);

%global parameter structure
global param

%paths for knitro and ampl
param.ampl_command='ampl';
%Need to have correct env variable for knitro to be able to run under ampl!!!
setenv('DYLD_LIBRARY_PATH','/usr/local/knitro900/lib') 

%original parameters for Judd and Su code
param.beta=0.975;
param.nT = 120;
param.nBus = 50;
param.N = 175;
param.M = 5;
param.RC = 11.7257;
param.thetaCost = 2.4569;
param.thetaProbs = [ 0.0937
       0.4475
       0.4459
       0.0127
       0.0002 ];

%system and MC parameters
param.logfile=sprintf('RustBusMLETableX_MC_Multistart_beta%3.0f.out',1000*param.beta);
param.MC=2;	
param.multistarts=2;
param.figure=0;
param.runNFXP=0;
param.runMPECktrlink=0;

%call Judd and Su code
cd juddsu/
delete(param.logfile);
RustBusMLETableX_MC
cd ..


