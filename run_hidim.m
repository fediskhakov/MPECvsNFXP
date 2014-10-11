% This script runs a high dimensional example
clear all

%add path to the functions
currdir=pwd;
addpath([currdir filesep 'matlab_include']);
addpath([currdir filesep 'nfxp_Rust87']);
addpath([currdir filesep 'high_dim']);
addpath('..');
cd high_dim

% paths for knitro and ampl
param.ampl_command='ampl';
%Need to have correct env variable for knitro to be able to run under ampl!!!
setenv('DYLD_LIBRARY_PATH','/usr/local/knitro900/lib') %Fedor
% setenv('DYLD_LIBRARY_PATH','/usr/local/ampl/knitro900/lib') %Bertel

% sys switches
param.figure=0;					%should figures be drawn

% parameters
param.MC=10; %number of MC simulated samples
param.nT = 120;
param.nBus = 50;
param.beta=0.975;
param.RC=11.7257;
param.thetaCost0=2.45569; %cost corresponding to John's grid size

% Parametric approximation by censored normal - discretized distribution
param.jr.alpha=3.2889;           % Mean of monthly mileage in Rust data set
param.jr.sigma_dx=1.4686;        % Std of monthly mileage in Rust data set
param.jr.upperbnd=450;


%MAIN LOOP over dimensions of the grid
dimvec=[200:100:2000];
%dimvec=[2100:200:5000];
param.out_suff='_1';
param.dontAMPL=0; % 1 to skip AMPL, only compute X0
param.dontNFXP=0; % 1 to skip AMPL, only compute X0

mvec=zeros(size(dimvec));
for idim=1:numel(dimvec)
	%number of grid points over milage
	param.N=dimvec(idim);
	
	%compute thetaProbs from truncated Normal fitted on John's data
	[param.thetaProbs, y0]=discretized_normal((0:param.N-1)'*param.jr.upperbnd/param.N,param.jr.alpha,param.jr.sigma_dx);
	param.thetaProbs(param.thetaProbs<1e-12)=[];
	% param.thetaProbs=[0.0937 0.4475 0.4459 0.0127 0.0002]'; %this is for testing
	%number of non-zeros in rows of transition probability matrix
	param.M=numel(param.thetaProbs);
	mvec(idim)=param.M;

	%adjust cost
	param.thetaCost=param.thetaCost0*(175/param.N);

	%plot of transition probabilities
	if 0
		h=figure('Color',[1 1 1],'NextPlot','new');
		ax=axes('Parent',h,'XTickLabel', ...
			{'0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25'}, ...
			'XTick',[1:numel(param.thetaProbs)]);
		box(ax,'on');
		hold(ax,'all');
		bar(param.thetaProbs,'Parent',ax);
	end

	%solve the model
	myEV=trueEV(param);	

	%simulate data (all MC samples)
	[data.MC_dt,data.MC_xt,data.MC_dx]=simdata(param,myEV);

	%run AMPL solver (and create X0)
	ampl{idim}=MPEC_AMPL(param,data);

	%NFXP
	if ~param.dontNFXP
	nfxp{idim}=run_nfxp_hidim(param,data,ampl{idim}.X0);
	end

	% Result table (appended on each iteration)
	% RTable(param,[dimvec;mvec],ampl,nfxp,1);
	RTable(param,[dimvec;mvec],ampl,nfxp,0); %no plot

end % end of loop over dimensions

cd ..


