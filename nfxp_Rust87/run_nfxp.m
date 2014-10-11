function out=run_nfxp(param,MC_dt,MC_xt,MC_dx,X0)	
%This funtion runs NFXP-NK estimator for the series of dataset
%
%INPUT: structure of common parameters
%				MC_* simulated data
%				X0 matrix of starting values (first 7 rows are used, columns=multistarts)
%OUTPUT: self-explaining structures of results for Rust87 implementation (golden standard) 
%				 and Judd-Su design replicated exactly

%Convert parameters to our notation
nMC=param.MC;
X.n=param.N;
X.grid= (0:X.n-1)';
S.N=param.nBus;
S.T=param.nT;

%parameters of algorithms
fxp_opt.printfxp=0;		% 0: No output,1: Fixed point iteration output 1   
ml_opt.title='Structural parameters in Rust(1987) Engine replacement model';
ml_opt.output=0;     % 0: No output, 1: Only final output,2: Detailed iteration output
ml_opt.method='Full MLE';  % 'Full MLE' or 'Partial MLE'

out.Runs=0;
out.TotalSuccess=0;

%Monte Carlo loop
for i_mc=1:param.MC

	%convert simulated data to Bertel's format
	datasim.d=reshape(MC_dt(:,:,i_mc)+1,S.T*S.N,1);
	datasim.x=reshape(MC_xt(:,:,i_mc),S.T*S.N,1);
	datasim.dx1=reshape(MC_dx(:,:,i_mc)-1,S.T*S.N,1);

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
		[mphat, mpse, cov, g, llval, iterinfo1]=nfxp.maxlik(ll_p, mp_start, pnames, ml_opt);
		% Full likelihood (Step 2)
		pnames={'RC','c','p'};
		[mphat, mpse, cov, g, llval, iterinfo2]=nfxp.maxlik(ll_f, mphat, pnames, ml_opt);
		timetoestimate=toc(estimationtimer);
		%save results
		j=(i_mc-1)*param.multistarts+reps; %index of the row in the results structure
		out.param=param;
		out.mc_iter(j)=i_mc;
		out.multistart_iter(j)=reps;
		out.Runs=out.Runs+1;
		out.TotalSuccess=out.TotalSuccess+iterinfo2.convflag;
		out.converged(j)=iterinfo2.convflag;
		%estimated parameters
		out.RC.est(j)=mphat.RC;
		out.RC.se(j)=mpse.RC;
		out.c.est(j)=mphat.c;
		out.c.se(j)=mpse.c;
		out.p.est(j,1:numel(mphat.p)+1)=[reshape(mphat.p,1,[]) 1-sum(mphat.p)];
		% out.p.se(j,1:numel(mphat.p)+1)=[reshape(mpse.p,1,[]) 0];
		%performance indicators
		out.runtime(j)=timetoestimate;
		out.MajorIter(j)=iterinfo1.MajorIter+iterinfo2.MajorIter;
		out.FuncEval(j)=iterinfo1.ll+iterinfo2.ll;
		out.BellmanIter(j)=iterinfo1.BellmanIter+iterinfo2.BellmanIter;
		out.NKIter(j)=iterinfo1.NKIter+iterinfo2.NKIter;

	end
end

if 1
	%text output
	c_tval	= (out.RC.est-param.RC) ./out.RC.se;
	RC_tval	= (out.c.est-param.thetaCost)./out.c.se;
	fprintf('\n');
	fprintf('---------------------------------------------------------------------------------------------------------\n');
	fprintf('Beta = %10.5f \n',param.beta);
	fprintf('RC = %10.5f \n',param.RC);
	fprintf('c = %10.5f \n',param.thetaCost);
	fprintf('Number of grid points = %g \n',param.N);
	fprintf('---------------------------------------------------------------------------------------------------------\n');
	fprintf('Rust 87 algorithm (2 stage likelihood, frequency based starting values\n');
	fprintf('---------------------------------------------------------------------------------------------------------\n');
	fprintf('Mean CPU time =%10.5f\n', mean(out.runtime)); 
	fprintf('NFXP converged %d times out of %d runs\n', out.TotalSuccess,out.Runs) ;
	fprintf('Average number of BHHH Hill Climbing iterations    %-g\n',   mean(out.MajorIter)) ;
	fprintf('Average number of likelihood evaluations           %-g\n',   mean(out.FuncEval)) ;
	fprintf('Average number of contraction iterations           %-g\n',   mean(out.BellmanIter)) ;
	fprintf('Average number of Newton Kantorowich iterations    %-g\n\n', mean(out.NKIter)) ;
	fprintf('        Bias       MCSD       Mean S.E    t(.025)     t(.25)    t(.50)     t(.75)    t(.975)\n'); 
	fprintf('RC   %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f  \n',  ...
			mean(out.RC.est-param.RC), std(out.RC.est-param.RC), mean(out.RC.se), quantile(RC_tval,[.025 .25 .50 .75 .975])) ;
	fprintf('c    %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f  \n',  ...
			mean(out.c.est-param.thetaCost ), std(out.c.est-param.thetaCost ), mean(out.c.se), quantile(c_tval,[.025 .25 .50 .75 .975])) ;
	fprintf('---------------------------------------------------------------------------------------------------------\n');
end


end %function



