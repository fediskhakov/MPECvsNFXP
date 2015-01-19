% NFXP: procedure to compute contraction operator
%
%       NFXP.bellman :
%           Procedure to compute bellman equation
%
% See also: nfxp.bellman, nfxp.dbellman, nfxp.solve, nfxp.maxlik

classdef nfxp
    methods (Static)
        function [data] = simdata(X, mp, S, P, pk)
            N=S.N;
            T=S.T;
            
            if nargin ==3  % solve model if solution is not provided
                fxp_opt.printfxp=0;
                P = nfxp.statetransition(mp.p, X.n);
                cost=0.001*mp.c*X.grid;
                [ev, pk]=nfxp.solve(0, P, cost, mp, fxp_opt);
            end
            
            id=repmat((1:N),T,1);
            t=repmat((1:T)',1,N);
            u_dx=rand(T,N);
            u_d=rand(T,N);
            
            csum_p=cumsum(mp.p);
            dx1=0;
            for i=1:numel(csum_p);
                dx1=dx1+ (u_dx>(csum_p(i)));
            end;
            
            x=nan(T, N);
            x1=nan(T, N);
            x(1,:)=ones(1, N); % Intial conditions
            
            for it=1:T;
                d(it,:)=(u_d(it,:)<(1-pk(x(it,:)')'))+1;  % Keep =1, Replace =2
                x1(it,:)=min(x(it,:).*(d(it,:)==1) +(d(it,:)==2) + dx1(it,:), X.n);
                if it<T;
                    x(it+1,:)=x1(it,:);
                end
            end
            
            data.id=id;
            data.t=t;
            data.d=d;
            data.x=x;
            data.x1=x1;
            data.dx1=dx1;
            
            pfields=fieldnames(data);
            for i=1:numel(pfields);
                data.(pfields{i})=reshape(data.(pfields{i}), T*N,1);
            end
            
            
            
        end
        
        function data = readdata(X, bustypes)
            %Read data
            dta=load('busdata.txt','-ascii');
            disp('got here');
            % Cols in busdata1234:
            %     1: busid,
            %     2: bustype
            %     3: year
            %     4: month
            %     5: replacement_dummy
            %     6: accumulated mileage, x
            %     7: lagged accumulated milage, x1
            
            % Select busses
            if nargin>1
                dta=dta(dta(:,2) <=bustypes,:);
            end
            
            % Populate data struct with selected busdata
            data.id=dta(:,1);
            data.d=dta(:,5)+1; % d=1: keep engine; d=2: replace engine
            
            data.x=dta(:,6)/1000;     % Discretized accumulated mileage, t
            data.x1=dta(:,7)/1000;    % Discretized accumulated mileage, t-1
            data.dx1=data.x1-(data.d==1).*data.x;   % Discretized monthly mileage
        end
        
        function [data, grid] = discretize(data, X)
            % Discrete dx and compute x from discretized x
            discretize = @(x) ceil(x*X.n/X.max);
            data.x=discretize(data.x);       % Discretized accumulated mileage, t
            data.x1=discretize(data.x1);     % Discretized accumulated mileage, t-1
            data.dx1=(data.x1-(data.d==1).*data.x);   % Discretized monthly mileage
        end
        
        function P = statetransition(p, n)
            p=[p; (1-sum(p))];
            P=0;
            for i=0:numel(p)-1;
                P=P+sparse(1:n-i,1+i:n,ones(1,n-i)*p(i+1), n,n);
                P(n-i,n)=1-sum(p(1:i));
            end
            P=sparse(P);
        end
        
        
        function [ev]=E(v, p)
            ev=p(1)*v;
            for i=1:numel(p)-1
                ev=ev+ p(i+1)*[v(1+i:end); repmat(v(end), i, 1)];
                
            end
        end
        
        
        function [ev1, pk]=bellman(ev, P, c, mp)
            % NFXP.BELMANN:     Procedure to compute bellman equation
            %
            % Syntax :          ev=nfpx.bellman(ev, P, c, mp);
            %
            % Inputs
            %  ev0              n x 1 matrix of expected values given initial
            %                   guess on value function
            % Output:
            %  ev1              n x 1 matrix of expected values given initial
            %                   guess of ev
            % See also:
            %   nfxp.maxlik, nfxp.solve, nfxp.dbellman
            VK=-c+mp.beta*ev;
            VR=-mp.RC-c(1)+mp.beta*ev(1);
            maxV=max(VK, VR);
            ev1=P*(maxV + log(exp(VK-maxV)  +  exp(VR-maxV)));
            if nargout>1
                %also compute choice probability from ev (initial input)
                pk=1./(1+exp((VR-VK)));
            end
        end % end of NFXP.bellman
                
        
        function dev=dbellman(pk,  P, mp)
            % NFXP.DBELMANN:     Procedure to compute Frechet derivative of Bellman operator
            % Syntax :          ev=nfpx.dbellman(pk,  P, mp);
            n=numel(pk);
            tmp=P(:,2:n).*repmat(pk(2:n,1)',n,1);
            dev=(mp.beta*[1-(sum(tmp,2)) tmp]);
        end
        
        
        function [ev, pk, F]=solve(ev, P, cost, mp, options)
            %NFXP.SOLVE: Contraction mapping fixed point poly-algorithm.
            %
            %  syntax:	[ Pk,EV,dEV]=nfxp.solve(ev, P, cost, mp, options):
            %
            %  INPUT:
            %     EV :      m x 1 matrix or scalar zero. Initial guess of choice specific expected value function, EV.
            %     P:        State transition matrix
            %     cost:     State transition matrix
            %     mp:       Structure containing model parameters: mp.beta, mp.RC, mp.c, mp.ev
            %     options:  Optional setting for Fixed Point Algorithm
            %
            %  OUTPUT:
            %     Pk:     m x 1 matrix of conditional choice probabilities, P(d=keep|x)
            %     EV:     m x 1 matrix of expected value functions, EV(x)
            %     F:      m x m matrix of derivatives of Identity matrix I minus 
            %             contraction mapping operator, I-T' where T' refers to derivative of the expected  value function
            %
            %  FOR OPTIONS: see nfxp.m
            %
            % See also:
            %   nfxp.maxlik, nfxp.bellman
            
            %-----------------------------------------------------------------------------------------------------------------------------
            % Default settings for Fixed point algorithm
            
            global BellmanIter NKIter;
            opt.max_fxpiter= 3;             % Maximum number of times to switch between Newton-Kantorovich iterations and contraction iterations.
            opt.min_cstp   = 4;             % Minimum number of contraction steps
            opt.max_cstp   = 20;            % Maximum number of contraction steps
            opt.ctol       = 0.02;          % Tolerance before switching to N-K algorithm
            opt.rtol       = 0.02;          % Relative tolerance before switching to N-K algorithm
            opt.nstep      = 20;            % Maximum number of Newton-Kantorovich steps
            opt.ltol0      = 1.0e-10;       % Final exit tolerance in fixed point algorithm, measured in units of numerical precision
            opt.printfxp   = 1;             % print iteration info for fixed point algorithm if > 0
            opt.rtolnk     = .5;             % Tolerance for discarding N-K iteration and move to SA (tolm1/tolm > 1+opt.rtolnk) 
            if nargin>4
                pfields=fieldnames(options);
                for i=1:numel(pfields);
                    opt.(pfields{i})=options.(pfields{i});
                end
            end
            n=numel(cost);
            
            %Initialize counters
            NKIter=0;
            BellmanIter=0;
            converged=false;%initialize convergence indicator
            tolm=1;

            for k=1:opt.max_fxpiter; %poli-algorithm loop (switching between SA and N-K and back)

                % SECTION A: CONTRACTION ITERATIONS
                if opt.printfxp>0
                    fprintf('\n');
                    fprintf('Begin contraction iterations (for the %d. time)\n',k);
                    fprintf('   j           tol        tol(j)/tol(j-1) \n');
                end;
                
                %SA contraction steps
                for j=1:opt.max_cstp;
                    ev1=nfxp.bellman(ev, P, cost, mp);

                    BellmanIter=BellmanIter+1;
                    
                    tolm1=max(abs(ev-ev1));
                    rtolm=tolm1/tolm;
                    
                    if opt.printfxp>0
                        fprintf(' %3.0f   %16.8f %16.8f\n',j, tolm1,rtolm);
                    end

                    %prepare for next iteration
                    ev=ev1;
                    tolm=tolm1;

                    %stopping criteria
                    if (j>=opt.min_cstp) && (tolm1 < opt.ctol)
                        %go to NK iterations due to absolute tolerance
                        break;
                    end;
                    if (j>=opt.min_cstp) && (abs(mp.beta-rtolm) < opt.rtol)
                        %go to NK iterations due to relative tolerance
                        break
                    end
                end
                %ev is produced after contraction steps
                %tolm will also be used below
            
                % SECTION 2: NEWTON-KANTOROVICH ITERATIONS
                if opt.printfxp>0
                    fprintf('\n');
                    fprintf('Begin Newton-Kantorovich iterations (for the %d. time)\n',k);
                    fprintf('  nwt          tol   \n');
                    
                end
                if opt.printfxp>1
                    %plots
                    hold on
                    subplot(2,1,2), plot(ev -ev(1), '--k');
                end

                %do initial contraction iteration which is part of first N-K iteration
                [ev1, pk]=nfxp.bellman(ev, P, cost, mp); %also return choice probs=function of ev

                for nwt=1:opt.nstep; %do at most nstep of N-K steps

                    NKIter=NKIter+1;

                    %Do N-K step
                    F=speye(n) - nfxp.dbellman(pk, P,  mp);%using pk from last call to nfxp.bellman
                    ev=ev-F\(ev-ev1); %resuing ev here
                    
                    %Contraction step for the next N-K iteration
                    [ev2, pk]=nfxp.bellman(ev, P, cost, mp); %also return choice probs=function of ev

                    %Measure the tolerance
                    tolm1=max(max(abs(ev-ev2)));

                    if opt.printfxp>0
                        fprintf('   %d     %16.8e  \n',nwt,tolm1);
                    end
                    if opt.printfxp>1
                        %plot ev
                        hold on
                        subplot(2,1,2), plot(ev -ev(1), '-r');
                    end

                    %Discard the N-K step if tolm1 got worse AND we can still switch back to SA
                    if (tolm1/tolm > 1+opt.rtolnk) && (k<opt.max_fxpiter);
                        if opt.printfxp>0
                            %Discrading the N-K step
                            fprintf('Discarding N-K step\n');
                        end
                        ev=ev1; %new contraction step should start from ev1
                        break;
                    else
                        ev1=ev2; %accept N-K step and prepare for new iteration
                    end;

                    %adjusting the N-K tolerance to the magnitude of ev
                    adj=ceil(log(abs(max(max(ev1)))));
                    ltol=opt.ltol0*10^adj;  % Adjust final tolerance
                    ltol=opt.ltol0;

                    if (tolm1 < ltol);
                        %N-K converged 
                        converged=true;
                        break
                    end

                end %Next N-K iteration

                if converged
                    if opt.printfxp>0
                        fprintf('Convergence achieved!\n\n');
                    end
                    break; %out of poly-algorithm loop
                else
                    if nwt>=opt.nstep
                        warning('No convergence! Maximum number of iterations exceeded without convergence!');
                        break; %out of poly-algorithm loop with no convergence
                    end
                end
            end
        end % end of solve        
        
        function [mphat, mpse, cov, g, ll, iterinfo]=maxlik(llfun, p0, pnames, options);
            %NFXP.MAXLIK:  Routine to maximize likelihood functions based on
            %               - Analytical gradients.
            %               - BHHH Algorithm
            %               - Line search (step-halfing)
            %
            %   Syntax:   [mphat, mpse, cov, g, ll, convflag]=maxlik(llfun, p0, pnames, options);
            %
            %
            %   OUTPUTS:
            %       mphat:      Estimated parameters (structure with as many fields as p0)
            %       mpse:       Standard errors of estimated parameters (structure with number of fields equal to number elements in pnames, nP)
            %       cov:        Variance Covariance Matrix (nP x nP matrix)
            %       g:          Score vectors (N x nP)
            %       ll:         Likelihood contributions (N x 1 vector)
            %       convflag    Bolean (true/false) idicator of convergence
            %
            %   INPUTS:
            %       llfun:      Function handle for likelihood function
            %
            %       p0:         Structure holding starting values of parameters
            %
            %       pnames:     k dimensional cell array, that holds names of the field in the parameter structure
            %                   to be estimated. Then likelihood is only maximized with respect to the parameters in pnames
            %                   and it is thus possible to fix parameters at the values set in p0
            %
            %                   Example:
            %                       Suppose parameter structure contains both p0.delta0 and p0.delta1
            %                       If pnames={'delta0','delta1'}, then both delta0 and delta1 are estimated
            %                       If pnames={'delta0'}, then only delta0 is estimated and delta1 is fixed at starting value (set in p0)
            %
            %
            %       options:    Structure holding voluntary options for ML. If unspecified default options are used
            %                   options structure has the following fields:
            %
            %           options.title (String):
            %                   A title for the model.
            %
            %           options.output (Scalar) :
            %                   0 if no output, 1 for final output, 2 for detailed iteration output (default is 2)
            %
            %           options.maxit (Scalar):
            %                   Maximum number of iterations (default is 100)
            %
            %           options.tolgdirec (Scalar):
            %                   Tolerance for gradient'direc (default is 1e-3)
            %
            %           options.lambda0:
            %                   The basic step length (default is 1).
            % See also:
            %   nfxp.solve, nfxp.bellman
            
            tt=tic;
            par.title={''};
            par.maxit=100;
            par.tolgdirec=1e-3;
            par.tol_bhhh=0.1;
            par.output=1;
            par.lambda0=1;
            par.hess=0;

            p1=p0;

            global BellmanIter NKIter
            iterinfo.MajorIter=0;
            iterinfo.ll=0;
            iterinfo.NKIter=0;
            iterinfo.BellmanIter=0;
                        
            if nargin==4
                pfields=fieldnames(options);
                for i=1:numel(pfields);
                    par.(pfields{i})=options.(pfields{i});
                end
            end
            
            lambda=par.lambda0;   % initialize stepsize
            k = length(pnames);  % Number of parameters to estimate
            
            % -------------------------
            %  ** BHHH Optimization **
            % -------------------------
            
            if (par.hess==1) 
                strhess='Analytical Hessian';
                [ll0,ev,s0, h]=llfun(p0);
                h=-h;
            elseif (par.hess==0);  
                strhess='Outer product of scores';
                [ll0,ev,s0]=llfun(p0);
                h=s0'*s0;
            end;

            iterinfo.ll=iterinfo.ll+1;
            iterinfo.NKIter=iterinfo.NKIter+NKIter;
            iterinfo.BellmanIter=iterinfo.BellmanIter+BellmanIter;


            g=sum(s0)';
            direc=h\g;
            tolm=g'*direc;
            if tolm<0; % convex area of likelihood, Newton-Raphson moves away from max
                h=s0'*s0;
                direc=h\g;
                tolm=g'*direc;
                if par.output>=1 
                    disp('Convex area of likelihood: Switch to BHHH')
                end
            end

            iterinfo.MajorIter=0;
            for it=0:par.maxit;
                for i=1:k;
                    np=numel(p1.(char(pnames(i))));
                    p1.(char(pnames(i)))=p0.(char(pnames(i)))+ lambda*direc(i:i+np-1);                    
                end
                p1.ev=ev;
                
                if lambda==par.lambda0;   % If previous step was accepted
                     if ((par.hess==0) | (abs(tolm) > par.tol_bhhh)); 
                        strhess='Outer product of the scores';
                        [ll,ev,s]=llfun(p1);
                        h=s'*s;
                    else
                        [ll,ev,s, h]=llfun(p1);
                        strhess='Analytical Hessian';
                        h=-h;
                    end;

                    iterinfo.ll=iterinfo.ll+1;
                    iterinfo.NKIter=iterinfo.NKIter+NKIter;
                    iterinfo.BellmanIter=iterinfo.BellmanIter+BellmanIter;


                else
                    [ll, ev, s]=llfun(p1); % Do we need to recalculate s, when line searching
                    iterinfo.ll=iterinfo.ll+1;
                    iterinfo.NKIter=iterinfo.NKIter+NKIter;
                    iterinfo.BellmanIter=iterinfo.BellmanIter+BellmanIter;


                end
                if sum(ll)<sum(ll0);
                    lambda=lambda/2;
                    if par.output>1;
                        fprintf('\n');
                        fprintf('.............LINESEARCHING\n');
                        fprintf('\n');
                    end
                    if lambda<=.01;
                        if par.output>1;
                            fprintf('.............WARNING: Failed to increase - downhill step accepted\n');
                        end
                        p0=p1;
                        ll0=ll;                 % Accept step
                        lambda=par.lambda0;
                    end
                else       % IF INCREASE IN LIKELIHOOD, CALCULATE NEW DIRECTION VECTOR
                    p0=p1;
                    ll0=ll;                 % Accept step
                    
                    
                    % ** Plot iteration info ***
                    if par.output>1;
                        fprintf('Iteration: %d   \n', iterinfo.MajorIter);
                        fprintf('grad*direc       %10.5f \n', g'*direc);
                        fprintf('||rel grad||     %10.5f \n', max(abs(g/sum(ll0))));
                        fprintf('Log-likelihood   %10.5f \n', sum(ll));
                        fprintf('Hessian update   %10s \n', strhess);
                        fprintf('Step length      %10.5f \n', lambda);
                        fprintf('---------------------------------------------------------------------\n');
                        fprintf('    Param.       Estimates      Direction         Gradient\n');
                        i=0;
                        for iP=1:numel(pnames);
                            for j=1:numel(p0.(char(pnames(iP))));
                                i=i+1;
                                fprintf('    %-10s %10.4f    %10.4f        %10.4f\n', char(pnames(iP)), p0.(char(pnames(iP)))(j), direc(i), g(i));
                            end;
                        end;
                        fprintf('---------------------------------------------------------------------\n');

                        fprintf('\n');
                    end
                    g=sum(s)';
                    iterinfo.MajorIter=iterinfo.MajorIter+1;
                    direc=h\g;
                    tolm=g'*direc;
                    tolgrad=max(abs(g/sum(ll0)));

                    if tolm<0; % convex area of likelihood, Newton-Raphson moves away from max
                        h=s'*s;
                        direc=h\g;
                        if par.output>=0
                            disp('Convex area of likelihood: Switch to BHHH')
                        end
                    end
                    lambda=par.lambda0; % and reset stepsize
                end
                
                if tolm < par.tolgdirec;  % Stopping rule
                    break;
                end
                
            end
            
            cov=inv(h); 
            se=(sqrt(diag(cov)));
            mphat=p0;
            i=0;
            for iP=1:numel(pnames);
                for j=1:numel(p0.(char(pnames(iP))))
                    i=i+1;
                    mpse.(char(pnames(iP)))(j)=se(i);
                end
            end
            
            % ---------------------------
            %  ** Plot final output  ***
            % ---------------------------
            if it<par.maxit;
                if par.output >= 1;
                    cr=corrcov(cov);
                    fprintf('---------------------------------------------------------------------\n');
                    fprintf('***                   Convergence Achieved                        ***\n');
                    fprintf('---------------------------------------------------------------------\n');
                    disp('                         _ ');
                    disp('                         \`\ ');
                    disp('                         |= | ');
                    disp('                        /-  ;.---. ');
                    disp('                  _ __.''     (____) ');
                    disp('                   `         (_____) ');
                    disp('                   _''  ._ .'' (____) ');
                    disp('                    `        (___) ');
                    disp('                   --`''------''` ');
                    fprintf('%s \n', char(par.title));
                    fprintf('Number of iterations: %d\n', iterinfo.MajorIter);
                    fprintf('grad*direc       %10.5f \n', g'*direc);
                    fprintf('Log-likelihood   %10.5f \n', sum(ll));
                    fprintf('Step length      %10.5f \n', lambda);
                    fprintf('\n');
                    fprintf('    %-13s %13s %13s %13s\n','Param.','Estimates','s.e.','t-stat');
                    fprintf('---------------------------------------------------------------------\n');
                     i=0;
                     for iP=1:numel(pnames);
                            for j=1:numel(p0.(char(pnames(iP))))
                                i=i+1;
                                myStr = char(pnames(iP));
                                if numel(p0.(char(pnames(iP))))>1; 
                                    myStr = sprintf('%-1s(%i)',char(pnames(iP)),j);
                                end;
                                fprintf('    %-13s %13.4f %13.4f %13.4f\n', myStr, p0.(char(pnames(iP)))(j), se(i), p0.(char(pnames(iP)))(j)/se(i));
                                mpse.(char(pnames(iP)))(j)=se(i);
                            end
                    end
                    fprintf('---------------------------------------------------------------------\n');

                    fprintf('\n');

                    fprintf('Correlation matrix of parameters\n');
                    fprintf('%18s','');
                    i=0;
                    for iP=1:k;
                        for j=1:numel(p0.(char(pnames(iP))))
                            i=i+1;
                            myStr = char(pnames(iP));
                            if numel(p0.(char(pnames(iP))))>1; 
                                myStr = sprintf('%-s(%i)',char(pnames(iP)),j);
                            end;

                            fprintf('%18s', myStr);
                        end;
                    end
                    fprintf('\n');

                    i=0;
                    for iP=1:k;
                        for j=1:numel(p0.(char(pnames(iP))))
                            i=i+1;
                            myStr = char(pnames(iP));
                            if numel(p0.(char(pnames(iP))))>1; 
                                myStr = sprintf('%-1s(%i)',char(pnames(iP)),j);
                            end;

                            fprintf('%-18s', myStr);
                            for j2=1:size(cr,1);
                                fprintf('%18.4f',cr(i,j2));
                            end
                            fprintf('\n');
                        end;
                    end
                    if ismac & par.output >= 3 
                        !say 'Success! Convergence Achieved'
                    end
                    tt_end=toc(tt);
                    fprintf('\n');
                    fprintf('Time to convergence is %3.0f min and %2.2f seconds\n', floor((tt_end)/60),tt_end - 60*floor((tt_end)/60));

                    fprintf('\n');
                    
                end
                % output a flag that we did converge
                iterinfo.convflag = true;
            else
                fprintf('----------------------------------------------------------------------------\n');
                fprintf('***   BHHH failed to converge: Maximum number of iterations is reached *** \n');
                fprintf('----------------------------------------------------------------------------\n');
                disp('                   _,....._ ')
                disp('                  (___     `''-.__ ')
                disp('                 (____ ')
                disp('                 (____ ')
                disp('                 (____         ___ ')
                disp('                      `)   .-''` ')
                disp('                      /  .'' ')
                disp('                     | =| ')
                disp('                      \_\ ')
                
                if ismac
                    !say 'B triple H failed to converge: Maximum number of iterations is reached without convergence. Better luck next time'
                end
                
                % no luck this time
                iterinfo.convflag = false;
                
            end

        end % end of estim.maxlik
        
    end % end of methods
    
    
end % end of estim class