function out=compare_lik(param,MC_dt,MC_xt,MC_dx,X0, options)	
%This funtion runs NFXP-NK estimator for the series of dataset
%The purpose of this function is to compare the optimized likelihood functions
%between NFXP and MPEC (and NFXP likelihood at MPEC estimated parameters)
%
%INPUT: structure of common parameters
%				MC_* simulated data
%				X0 matrix of starting values (first 7 rows are used, columns=multistarts)
%OUTPUT: self-explaining structures of results for Rust87 implementation (golden standard) 
%				 and Judd-Su design replicated exactly

homedir=pwd;
cd('juddsu')
load(['MC' num2str(param.MC) '_beta' num2str(100000*param.beta) '_result'])
cd(homedir)

% Use same optios as in run_nfxp.m
ml_opt=options.ml_opt;
fxp_opt=options.fxp_opt;


%Convert parameters to our notation
nMC=param.MC;
X.n=param.N;
X.grid= (0:X.n-1)';
S.N=param.nBus;
S.T=param.nT;

%Monte Carlo loop
for i_mc=1:param.MC

	%convert simulated data to Bertel's format
	datasim.d=reshape(MC_dt(2:S.T,:,i_mc)+1,(S.T-1)*S.N,1);
	datasim.x=reshape(MC_xt(2:S.T,:,i_mc),(S.T-1)*S.N,1);
	datasim.dx1=reshape(MC_dx(1:S.T-1,:,i_mc)-1,(S.T-1)*S.N,1);


	for reps=1:param.multistarts
		%multistarts only differ in starting values

		% Frequency estimator for transition probs
		tab=tabulate(datasim.dx1);
		tab=tab(tab(:,3)>0,:);
		p=tab(1:end-1,3)/100;
		P = nfxp.statetransition(p, X.n);
		% Initial and fixed values
		mp_start.beta=param.beta;
		mp_start.RC=X0(7,reps); %use JS initial values for the variables on the first stage
		mp_start.c=X0(1,reps);
		mp_start.ev=0;
		mp_start.p=p; %freq estimator

	  % Likelihood function wrappers
	    ll_p=@(mp) ll(datasim, X, P, mp, fxp_opt, 'pmle');
		ll_f=@(mp) ll(datasim, X, P, mp, fxp_opt, 'fullmle');

		% Partial likelihood for RC and c (Step 1)
		estimationtimer=tic;
   	    pnames={'RC','c'};
		[mphat, mpse, cov, g, llval_p, iterinfo1]=nfxp.maxlik(ll_p, mp_start, pnames, ml_opt);

		% Full likelihood (Step 2)
		pnames={'RC','c','p'};
		[mphat, mpse, cov, g, llval_f, iterinfo2]=nfxp.maxlik(ll_f, mphat, pnames, ml_opt);

		%save results
		j=(i_mc-1)*param.multistarts+reps; %index of the row in the results structure
		out.param=param;

		;

	    mp_js.beta=param.beta;
		mp_js.RC=thetaAMPLsol(7,i_mc); %use JS initial values for the variables on the first stage
		mp_js.c=thetaAMPLsol(1,i_mc);
		mp_js.ev=0;
		mp_js.p=thetaAMPLsol(2:6,i_mc); %freq estimator

		ll_nfxp(i_mc, reps)=sum(llval_f);
		out.ll_nfxp_js(i_mc)=sum(ll_f(mp_js));
	end
end
out.ll_nfxp=max(ll_nfxp,[],2)';
out.ll_MPEC_js=ObjValAMPL';

diary 'ComareLik.txt'
fprintf('Liklihood values\n');
fprintf('----------------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf('%18s %18s %18s\n','NFXP','NFXP(theta_MPEC)','MPEC(theta_MPEC)');
fprintf('----------------------------------------------------------------------------------------------------------------------------------------------------------\n');
for i_mc=1:param.MC;
	fprintf('%18.4f %18.4f %18.4f\n',out.ll_nfxp(i_mc),out.ll_nfxp_js(i_mc),out.ll_MPEC_js(i_mc));
end
fprintf('----------------------------------------------------------------------------------------------------------------------------------------------------------\n')
diary off




