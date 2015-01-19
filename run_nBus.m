% This script runs compares MPEC-AMPL to NFXP for various sample sizes
clear all

%add path to the functions
currdir=pwd;
addpath([currdir filesep 'matlab_include']);
addpath([currdir filesep 'nfxp_Rust87']);
addpath([currdir filesep 'nBus']);
addpath('..');
cd nBus

% paths for knitro and ampl
param.ampl_command='ampl';
%Need to have correct env variable for knitro to be able to run under ampl!!!
%setenv('DYLD_LIBRARY_PATH','/usr/local/knitro900/lib') %Fedor
setenv('DYLD_LIBRARY_PATH','/usr/local/knitro/lib') %Bertel

% sys switches
param.figure=0;					%should figures be drawn

% parameters
param.MC=1; %number of MC simulated samples
param.nT = 120;
param.N = 175;
param.beta=0.9999;
param.RC=11.7257;
param.thetaCost=2.45569; %cost corresponding to John's grid size
param.thetaProbs=[0.0937 0.4475 0.4459 0.0127 0.0002]';


%MAIN LOOP over dimensions of the grid
param.out_suff='_nBus';
param.dontAMPL=0; % 1 to skip AMPL, only compute X0
param.dontNFXP=0; % 1 to skip AMPL, only compute X0

param.beta=0.9999;

nBus=(100:100:10000);
mvec=zeros(size(nBus));
for iBus=1:numel(nBus)
	param.nBus = nBus(iBus);

	%number of non-zeros in rows of transition probability matrix
	param.M=numel(param.thetaProbs);
	mvec(iBus)=5;

	%solve the model
	myEV=trueEV(param);	

	%simulate data (all MC samples)
	[data.MC_dt,data.MC_xt,data.MC_dx]=simdata(param,myEV);

	%run AMPL solver (and create X0)
	ampl{iBus}=MPEC_AMPL(param,data);

	%NFXP
	if ~param.dontNFXP
		nfxp{iBus}=run_nfxp_nBus(param,data,ampl{iBus}.X0);
	end

	% Result table (appended on each iteration)
	% RTable(param,[dimvec;mvec],ampl,nfxp,1);
	RTable(param,nBus,ampl,nfxp,0); %no plot
	save('Results_nBus','ampl','nBus', 'mvec', 'param', 'nfxp');

end % end of loop over sample size

cd ..

for iBus=1:numel(nBus); disp([ampl{iBus}.runtime ampl{iBus}.runtime_matlab ]); end
plot_cpu_vs_samplesize



