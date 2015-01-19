function out=MPEC_AMPL(param,data)
%This function calls AMPL once 
%to estimate the model once for each MC sample
%Estimation: PARTIAL LIKELIHOOD (2 parameters--c and RC-- with thetaProbs fixed)

%unpack parameters
MC=param.MC;
beta = param.beta; 
nT = param.nT;
nBus = param.nBus;
N = param.N;
M = param.M;
RC = param.RC;
thetaCost = param.thetaCost;
thetaProbs = param.thetaProbs;
%unpack data
MC_dt=data.MC_dt;
MC_xt=data.MC_xt;

KnitroExitAMPL = -10000*ones(MC,1);
thetaCostAMPL = zeros(1,MC);
RCAMPL = zeros(1, MC);
EVAMPL = zeros(N, MC);
thetaProbsAMPL = zeros(M,MC);
SolveTimeAMPLMATLAB = zeros(MC,1);
SolveTimeAMPL = zeros(MC,1);
ObjValAMPL = zeros(MC,1);
IterAMPL = zeros(MC,1);
FunEvalAMPL = zeros(MC,1);
SuccessAMPL = zeros(MC,1);

% Starting Points AND fixed probability freq estimators for PML
% Frequecy estimator for transition probs
out.fixedProbs=zeros(M,MC);
for kk=1:MC
    tab=tabulate(reshape(data.MC_dx(:,:,kk),[],1)-1);
    tab=tab(tab(:,3)>0,:);
    out.fixedProbs(1+tab(:,1),kk)=tab(:,3)/100;
    % out.fixedProbs(:,kk)=param.thetaProbs;
    %starting values for parameters in AMPL
    %at true values
    out.X0(1,kk)=0*param.thetaCost; %c
    out.X0(2,kk)=0*param.RC; %RC
    out.X0(3,kk)=0; %EV
end
%true parameter value
thetatrue=[param.thetaCost;param.RC];

%possibility for early exit: starting values calculated, but MPEC is not run
if param.dontAMPL
    return;
end

% %Initial estimates for NFXP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 2) 
%%% Estimate the model using the constrained optimization approach with
%%% AMPL/KNITRO implementation;
%%%
%%% IMPLEMENTATION 1: MPEC/AMPL
%%% We implement the MPEC approach using AMPL modeling language
%%% with KNITRO as the solver.
%%% AMPL supplies first-order and second-order analytical derivatives
%%% and the sparsity pattern of the constraint Jacobian and the Hessian
%%% to the solver. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for kk = 1:MC
   
    fprintf('This is Monte Carlo run #%d out of %d\n \n', kk, MC);
    
    xt = MC_xt(:,:,kk);
    dt = MC_dt(:,:,kk);
        
       fprintf('Starting MPEC/AMPL Partial ML in Monte Carlo run #%d out of %d ... \n', kk, MC);
        
        fid = fopen('RustBusMLETableX.dat', 'w');
        fprintf (fid, '# Data generated: %s\n\n', datestr(now));  
        fprintf (fid, 'data; \n\n'); 
        fprintAmplParamCLSU(fid, 'N', N, 1);
        fprintAmplParamCLSU(fid, 'M', M, 1);
        fprintAmplParamCLSU(fid, 'nBus', nBus, 1);
        fprintAmplParamCLSU(fid, 'nT', nT, 1);
        fprintAmplParamCLSU(fid, 'beta', beta, 1);
        fprintAmplParamCLSU(fid, 'thetaProbs', out.fixedProbs(:,kk), 1);
        fprintAmplParamCLSU(fid, 'inithetaCost', out.X0(1,kk), 1);
        fprintAmplParamCLSU(fid, 'iniRC', out.X0(2,kk), 1);
        fprintAmplParamCLSU(fid, 'iniEV', out.X0(3,kk), 1);
        fprintAmplParamCLSU(fid, 'xt', xt, 1); 
        fprintAmplParamCLSU(fid, 'dt', dt, 1);
        %fprintAmplParamCLSU(fid, 'inithetaProbs', X0(2:1+numel(thetaProbs),:), 1);
        fclose(fid);
        % The ampl executable in my computer is in the directory:
        % "/Applications/AMPL/ampl64"
        % Change the path below to the location of your AMPL executable.
        strAmplCommand = param.ampl_command; 
        outname = ['output/MC' num2str(kk) '.out']; 
        strAmplSystemCall = sprintf('%s RustBusMLETableX_PML.run > %s', strAmplCommand, outname);
        startmapl=tic;
            [status,result] = system(strAmplSystemCall);
        matlabtimer=toc(startmapl);

        %check for success in AMPL run
        if exist('output/KnitroExit.sol')
            KnitroExitAMPL(kk) = csvread('output/KnitroExit.sol'); 
            delete('output/KnitroExit.sol'); %just in case so that next iteration is not confused

            if KnitroExitAMPL(kk) == 0
                SolveTimeAMPLMATLAB(kk)=matlabtimer;
                SuccessAMPL(kk) = 1;
                SolveTimeAMPL(kk) = csvread('output/solvetime.sol');
                ObjValAMPL(kk) = csvread('output/objval.sol');
                thetaCostAMPL(kk) = csvread('output/thetaCost.sol');
                RCAMPL(kk) = csvread('output/RC.sol');
                EVAMPL(:,kk) = csvread('output/EV.sol');

                fid = fopen('KnitroMessage.sol');
                line1 = fgetl(fid);
                line2 = fgetl(fid);
                line3 = fgetl(fid);
                IterAMPL(kk) = sscanf(line3, '%d',1);
                FunEvalAMPL(kk) = sscanf(line3, '%*s %*s %d',1);
                fclose(fid); 
                delete('KnitroMessage.sol'); %to preven the file being used again by mistake - if ample fails to run
                
                fprintf('\n');
                disp([line1])
                disp([line2])
                disp([line3])
                fprintf('\n');

            else
                SuccessAMPL(kk) = 0;
            end
        else
            error 'AMPL didn''t run';
        end

end

thetaAMPLsol = [thetaCostAMPL;RCAMPL];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 5
%%% Calculate summary statistics of the Monte Carlo experiments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out.TotalSuccess = sum(SuccessAMPL);
fprintf('MPEC/AMPL finishes successfully in #%d out of %d runs \n \n', out.TotalSuccess, MC);
out.mean_estimates = mean(thetaAMPLsol(:, KnitroExitAMPL==0),2);
% out.stdevthetaAMPL = std(thetaAMPLsol(:, KnitroExitAMPL==0),1,2);
out.bias_estimates = out.mean_estimates - thetatrue;
out.rmse_estimates = sqrt(mean((thetaAMPLsol(:, KnitroExitAMPL==0)- repmat(thetatrue,1,sum(KnitroExitAMPL==0))).^2,2));
out.ave_feval = mean(FunEvalAMPL);
out.run_time = mean(SolveTimeAMPL);
out.num_iter = mean(IterAMPL);
out.run_time = mean(SolveTimeAMPL);

out.run_time_conv = mean(SolveTimeAMPL(KnitroExitAMPL==0));
out.run_time_conv_matlab = mean(SolveTimeAMPLMATLAB(KnitroExitAMPL==0));
out.num_iter_conv = mean(IterAMPL(KnitroExitAMPL==0));
out.ave_feval_conv = mean(FunEvalAMPL(KnitroExitAMPL==0));

out.runtime = SolveTimeAMPL; %for each run
out.converged = SuccessAMPL; %for each run
out.iter = IterAMPL; %for each run
out.lleval = FunEvalAMPL; %for each run

out.runtime_matlab = SolveTimeAMPLMATLAB; %for each run


