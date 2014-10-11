function myEV=trueEV(param)	
%This funtion runs NFXP-NK solver to produce expected values for the data simulation
%Optionally, performs check against EV computed in JuddSu code
%INPUT: structure of common parameters
%OUTPUT: Judd-Su style EV

%Convert parameters to Bertel's notation
X.n=param.N;
X.grid= (0:X.n-1)';
S.N=param.nBus;
S.T=param.nT;

%parameters of algorithms
fxp_opt.printfxp=0;		% 0: No output,1: Fixed point iteration output 1   
ml_opt.title='Structural parameters in Rust(1987) Engine replacement model';
ml_opt.output=0;     % 0: No output, 1: Only final output,2: Detailed iteration output
ml_opt.method='Full MLE';  % 'Full MLE' or 'Partial MLE'

%intiial solver for expected values to be used for data simulations
mp0.p=param.thetaProbs(1:end-1);
mp0.RC=param.RC;
mp0.c=param.thetaCost;
mp0.beta=param.beta;
P0 = nfxp.statetransition(mp0.p,X.n);

%costs
cost0=0.001*mp0.c*X.grid;  

%run solver
[myEV, pr2keep]=nfxp.solve(0, P0, cost0, mp0, fxp_opt);

if param.figure
	%check against Judd-Su calculations of EV

	try
		load(['juddsu/truethetaEV_beta' num2str(1000*param.beta)]);
		x = (1:param.N)';
		P0 = 1./ (1 + exp( 0.001*thetaCost*x - param.beta.*EV - RC - 0.001*thetaCost*x(1)+ param.beta*EV(1))); %JS choice probs

		fig=figure('Color',[1 1 1],'NextPlot','new');
		ax1=subplot(1,2,1,'Parent',fig);
		ax2=subplot(1,2,2,'Parent',fig);
		%ev
		plot(ax1,EV);
		hold(ax1,'on');
		plot(ax1,myEV);
		%prob
		plot(ax2,P0);
		hold(ax2,'on');
		plot(ax2,pr2keep);

		title(ax1,'Value functions');
		title(ax2,'Choice probabilities');
	catch err
		warning 'Judd-Su computed true EV file not found, no graphs will be plotted'
	end
		
end

