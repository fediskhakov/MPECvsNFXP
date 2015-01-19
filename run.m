% This is the main script to create all the data for the
% Iskhakov,Rust,Schjerning comment for Su and Judd 2012 Econometrica
clear all

%add path to the functions
currdir=pwd;
addpath([currdir filesep 'matlab_include']);
addpath([currdir filesep 'nfxp_Rust87']);

% global parameter structure
global param

% paths for knitro and ampl
param.ampl_command='ampl';
%Need to have correct env variable for knitro to be able to run under ampl!!!
%setenv('DYLD_LIBRARY_PATH','/usr/local/knitro900/lib') %Fedor
setenv('DYLD_LIBRARY_PATH','/usr/local/knitro/lib') %Bertel

if 1 %RUN MC or read from saved file

% original model parameters for Su and Judd code (Rust 87)
% param.nT = 120;
% param.nBus = 50;
% param.N = 175;
% param.M = 5;
% param.beta=0.975;
% param.RC = 11.7257;
% param.thetaCost = 2.4569;
% param.thetaProbs = [ 0.0937
%        0.4475
%        0.4459
%        0.0127
%        0.0002 ];

% alternative set of true parameters
param.nT = 120;
param.nBus = 500;	    %larger sample (60,000) for Table 2
param.nBus = 50;		%normal sample size (6000) for Table 1
param.N = 175;
param.M = 5;
param.beta=0.995;
param.RC=11.7257;
param.thetaCost=2.45569;
param.thetaProbs=[0.0937 0.4475 0.4459 0.0127 0.0002]';

%system and MC parameters
param.MC=250;					  %number of MC iterations
param.multistarts=5;		%number of tries for each estimation
param.runInitEV=0;			%run AMPL to compute EV of true parameter values
param.figure=0;					%should figures be drawn
param.runMPECktrlink=1;	%run MPEC/ktrlink in JuddSu code
param.runNFXP=0;				%run NFXP-SA in JuddSu code?
param.FreqX0=1;					%use frequencies for starting values of prob params
param.logfile=sprintf('RustBusMLETableX_MC_Multistart_beta%3.0f.out',100000*param.beta);
param.runNKktrlink=0;		%run NFXP-NK with ktrlink?


%Loop over betas
betavec=[0.975 0.985 0.995 0.999 0.9995 0.9999];
for ibeta=1:numel(betavec)

	%update true parameter value
	param.beta=betavec(ibeta);

	%compute initial EV with true parameter values using NFXP if asked for
	%this is needed because AMPL fails for solve model with high beta
	if ~param.runInitEV 
		EV=trueEV(param);
	end

	%call MPEC in Judd-Su directory
	cd 'juddsu/'
	delete(param.logfile);
	RustBusMLETableX_MC
	cd '..'
	result_js{ibeta}=result;
	%clean up workspace to ensure no accidentelly matching var names
	clearbut('param','result_js','result_jr87','result_nk', ...
					 'ibeta','betavec','MC_dt','MC_xt','MC_dx','X0');

	%load simulated data created in Judd-Su code
	if ~exist('MC_xt') || ~exist('MC_dt') || ~exist('MC_dx')
		simdatafile=['RustBusTableXSimDataMC' num2str(param.MC) '_beta' num2str(100000*param.beta) '.mat'];
		load(['juddsu' filesep simdatafile]);
		disp(['Simulation data loaded from ' simdatafile]);
		clearbut('param','result_js','result_jr87','result_nk','ibeta',...
						 'MC_dt','MC_xt','MC_dx','X0');
	end
	%call NFXP
	result_jr87{ibeta}=run_nfxp(param,MC_dt,MC_xt,MC_dx,X0);


	%create intermediate tables (for each iteration of beta)
	ResultTables(param,betavec,result_js,result_jr87,false);

end %beta loop

%save all results
save('GrandResultsMC','MC_dt','MC_dx','MC_xt','X0','betavec','param','result_jr87','result_js');

else %RUN MC or read from saved file
load 'GrandResultsMC.mat'	
end %RUN MC or read from saved file

%create all tables and graphs
ResultTables(param,betavec,result_js,result_jr87,param.figure);
%ibeta=1;
%comp_lik{ibeta}=compare_lik(param,MC_dt,MC_xt,MC_dx,X0,result_jr87{ibeta});



