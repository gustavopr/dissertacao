% *************************************************************************
% STRNKPC
% *************************************************************************
% Last update: July 16, 2008
% By: Waldyr Areosa; Marcelo C. Medeiros
% *************************************************************************
% Program for
% INSTRUMENTAL VARIABLE ESTIMATION OF SMOOTH TRANSITION REGRESSION MODELS
% WITH ENDOGENOUS VARIABLES
% By
% Waldyr Areosa; Michael McAleer; Marcelo C. Medeiros
% *************************************************************************
% Descripition
% Estimates a Smooth Transition Regression (STR) New-Keynesian Phillips
% Curve (NKPC) for the Brazilian Economy
% *************************************************************************
% Inputs: 
% -------
% Function Data.m:
%   y      : Regressand
%   xex    : Exogenous Regressors
%   xen    : Endogenous Regressors
%   z      : Logistic Function (Transition Variable)
%   wx     : Instruments for x
%   wz     : Instruments for z
%   ym     : Periods
%
% Outputs:
% --------
% Function LinTest.m (called in data.m)
%   Linearity test (F, f, isH0, pvalue)
%
% Linear model
% OLS:
%   alphaO        : parameter estimates
%   sd_hatO       : standard values
%
% 2SLS   
%   alphaL        : parameter estimates
%   sd_hatL       : standard values
% 
% Function Estimate.m
% STR model
% Standard Instruments
%   alphaf        : linear estimates (= par(1:nX,1))
%   sd_hat(1:nX,1): standard values
%
%   lambdaf       : nonlinear estimates (=par(nX+1:2*nX,1))
%   sd_hat(nX+1:2*nX,1): standard values
%
%   par(1:nX,1)+par(nX+1:2*nX,1): total effect
%
%   [gammaf cf]    : logistic parameters (=par(2*nX+1:2*nX+2,1))
%   sd_hat(2*nX+1:2*nX+2,1): standard values
%
% Optimal Instruments
%   alphaf        : linear estimates (= par2(1:nX,1))
%   sd_hat2(1:nX,1): standard values
%
%   lambdaf       : nonlinear estimates (=par2(nX+1:2*nX,1))
%   sd_hat2(nX+1:2*nX,1): standard values
%
%   par2(1:nX,1)+par2(nX+1:2*nX,1): total effect
%
%   [gammaf cf]    : logistic parameters (=par2(2*nX+1:2*nX+2,1))
%   sd_hat2(2*nX+1:2*nX+2,1): standard values
%
% Functions required:
% -------------------
% 1. Data.m          : Data input
% 2. LinTest.m       : Linearity test
% 3. Estimate.m      : Estimate
%    3.1. IVCost.m   : Minimization criteria
%
%
% *************************************************************************
% 1. PARAMETERS SETUP
% *************************************************************************
flag         = 0;                      % =1 do not drop intercept
                                       % <>1 drop intercept
                                       
sig          = 0.95;                   % significance level

instr_Amemya = 0;                      % =0 uses Standart Instruments
                                       % <>0 uses Optimal Instruments (Amemya)


% *************************************************************************
% 2. READ DATA
% *************************************************************************
% Call function Data.m
% ---------------------
[y,xex,xen,z,wx,wz,ym] = Data;


% *************************************************************************
% 3. LINEARITY TEST
% *************************************************************************
% Call function LinTest.m
% ------------------------
[F, f, isH0, pvalue] = LinTest(y,xex,xen,z,wx,wz(:,1),flag,sig,instr_Amemya)


% *************************************************************************
% 4. ESTIMATION
% *************************************************************************
% Prepare data
% ------------
[T,neq] = size(y);                     % sample size & # equations             
    
Y       = y;                           % regressand(s)
X       = [ones(T,1) xex xen];         % regressor(s)
Z       = [ones(T,1) wx wz];           % instrumets
q       = z;                           % transition variable            
qZ      = wz(:,1);                     % transition variable instrument

nxex    = size(xex,2)+1;               % # exogenous regressors
nxen    = size(xen,2);                 % # endogenous regressors
nX      = size(X,2);                   % # regressors
nZ      = size(Z,2);                   % # instruments

% -------------------------------------------------------------------------
% 4.1. LINEAR ESTIMATES (OLS & 2SLS)
% -------------------------------------------------------------------------
W   = [Z];                         % instrument set 1 (IS1) or instrument
                                    % set 3 (IS3)
%W   = [Z Z(:,2:end).^2];            % instrument set 2 (IS2) or IS4
nW  = size(W,2);                    % size of expanded set of instruments
PW  = W*pinv(W'*W)*W';              % instruments' projection matrix

% 4.1.1. OLS
% ----------
alphaO   = inv(X'*X)*X'*Y;          % linear OLS estimates
YhatO    = X*alphaO;                % linear OLS projection
ehatO    = y - YhatO;               % linear OLS error
v_ehatO  = var(ehatO);              
sd_hatO  = sqrt(diag(v_ehatO*inv(X'*X)));
%t_statO  = alphaO./sd_hatO;
%pvalueO  = 1-tcdf(t_statO,T-nX);

% 4.1.2. 2SLS
% -----------
alphaL   = inv(X'*PW*X)*X'*PW*Y;    % linear 2SLS estimates
YhatL    = X*alphaL;                % linear 2SLS projection
ehatL    = y - YhatL;               % linear 2SLS error
v_ehatL  = var(ehatL);              
sd_hatL  = sqrt(diag(v_ehatL*inv(X'*PW*X)));
%t_statL  = alphaL./sd_hatL;
%pvalueL  = 1-tcdf(t_statL,T-nX);

% 4.1.3. TEST J
% ---------------
m = W'*ehatL/T;
M = -W'*X/T;
S = v_ehatL*(W'*W)/T;
%PM = M*pinv(M'*M)*M';
%cov_m = [eye(size(PM)) - PM]*S*[eye(size(PM)) - PM]'/T;

% Model under the null
% --------------------
J0   = chi2inv(sig,size(W'*ehatL,1)- size(X,2));

% Model under the alternative
% ---------------------------
JA = T*m'*pinv(S)*m;
%JA = ehatL'*PW*ehatL/v_ehatL;

% Calculate and Compare Statistics
% --------------------------------
J0 = J0';
isH0J = JA <= J0;
pvalueJ=1-chi2cdf(JA,size(W'*ehatL,1)- size(X,2));
disp('#m #par J0 JA isH0J pvalueJ')
[size(W'*ehatL,1) size(X,2) J0 JA isH0J pvalueJ]


% -------------------------------------------------------------------------
% 4.2. NONLINEAR ESTIMATES (STANDART & OPTIMAL INSTRUMENTS)
% -------------------------------------------------------------------------
%
% 4.2.1. STANDART INSTRUMENTS (NLIV or MNLIV)
% -------------------------------------------
% Call function Estimate.m
% -------------------------
model = eye(2*nX+2);
ord = model*(1:2*nX+2)';
PV = ones(2*nX+2,1);

[par1, sd_hat1, VM1] = Estimate(y,xex,xen,z,W,PV,ord,model,0);
%[J0,JA,isH0J,pvalueJ,F1,f1,isH01,pvalue1,F2,f2,isH02,pvalue2] = diagnostics(par1,X,y,z,W,VM1);

% 4.2.2. OPTIMAL INSTRUMENTS (NLIV or MNLIV)
% ------------------------------------------
% Optimal Instruments
% -------------------
%W2   = [Z Z(:,2:end).^2 Z(:,2:end).^3 Z(:,2:end).^4 Z(:,2:end).^5];
W2 = Z;
for i=1:size(Z,2)
    for j=i:size(Z,2)
    W1 = Z(:,i).*Z(:,j);
    W2 = [W2 W1];
    end
end
nW2  = size(W2,2);                    
PW2  = W2*pinv(W2'*W2)*W2';              
%Ghat1 = PW2*Ghat;

% % Call function Estimate.m
% -------------------------
[par2, sd_hat2, VM2] = Estimate(y,xex,xen,z,W2,PV,ord,model,par1);
%[J0,JA,isH0J,pvalueJ,F1,f1,isH01,pvalue1,F2,f2,isH02,pvalue2] =
%diagnostics(par2,X,y,z,W2,VM2);

% *************************************************************************
% 5. SHOW RESULTS
% *************************************************************************
disp('Lineariy Test [F f isH0 pvalue]:')
[F  f  isH0  pvalue]

disp('OLS Linear model [Parameters SD Par/SD t-stat pvalue]:')
[alphaO sd_hatO alphaO./sd_hatO]

disp('2SLS Linear model [Parameters SD Par/SD t-stat pvalue]:')
[alphaL sd_hatL alphaL./sd_hatL]

disp('STR model - linear estimates [Parameters SD Par/SD]:')
[par1(1:nX,1) sd_hat1(1:nX,1) par1(1:nX,1)./sd_hat1(1:nX,1)]

disp('STR model - nonlinear estimates [Parameters SD Par/SD]:')
[par1(nX+1:2*nX,1) sd_hat1(nX+1:2*nX,1) par1(nX+1:2*nX,1)./sd_hat1(nX+1:2*nX,1)]

disp('STR model - total effect (vs Linear Model):')
[par1(1:nX,1)+par1(nX+1:2*nX,1) alphaL]

disp('Logistic [Parameters SD Par/SD]:')
[par1(2*nX+1:2*nX+2,1) sd_hat1(2*nX+1:2*nX+2,1) par1(2*nX+1:2*nX+2,1)./sd_hat1(2*nX+1:2*nX+2,1)]

disp('Optimal instruments - linear estimates [Parameters SD Par/SD]:')
[par2(1:nX,1) sd_hat2(1:nX,1) par2(1:nX,1)./sd_hat2(1:nX,1)]

disp('Optimal instruments - nonlinear estimates [Parameters SD Par/SD]:')
[par2(nX+1:2*nX,1) sd_hat2(nX+1:2*nX,1) par2(nX+1:2*nX,1)./sd_hat2(nX+1:2*nX,1)]

disp('Optimal instruments - total effect (vs Linear Model):')
[par2(1:nX,1)+par2(nX+1:2*nX,1) alphaL]

disp('Optimal instruments - Logistic [Parameters SD Par/SD]:')
[par2(2*nX+1:2*nX+2,1) sd_hat2(2*nX+1:2*nX+2,1) par2(2*nX+1:2*nX+2,1)./sd_hat2(2*nX+1:2*nX+2,1)]


%**************************************************************************
% 6. PLOT RESULTS
%**************************************************************************
%figure('Name','Residual Autocorrelation')
%autocorr(ehat)
%print -deps 'f2_IS1.eps'       % MNLIV (Modified NLIV) + IS1 
%print -deps 'f2_IS2.eps'       % MNLIV + IS2
%print -deps 'f2_IS3.eps'       % MNLIV + IS3
%print -deps 'f2_IS4.eps'       % MNLIV + IS4
%print -deps 'f2_IS1s.eps'      % NLIV (Standard Nonlinear IV) + IS1 
%print -deps 'f2_IS2s.eps'      % NLIV + IS2
%print -deps 'f2_IS3s.eps'      % NLIV + IS3
%print -deps 'f2_IS4s.eps'      % NLIV + IS4
%autocorr(ehat2)
%print -deps 'f2_IS1o.eps'      % MNLIV* (Optimal MNLIV) + IS1 
%print -deps 'f2_IS2o.eps'      % MNLIV* + IS2 
%print -deps 'f2_IS3o.eps'      % MNLIV* + IS3 
%print -deps 'f2_IS4o.eps'      % MNLIV* + IS4 
%print -deps 'f2_IS1so.eps'      % NLIV* (Optimal NLIV) + IS1 
%print -deps 'f2_IS2so.eps'      % NLIV* + IS2 
%print -deps 'f2_IS3so.eps'      % NLIV* + IS3 
%print -deps 'f2_IS4so.eps'      % NLIV* + IS4

%figure('Name','Logistic Function')
%subplot(121), plot(q,fX,'.'), grid, title('Impact of transition variable on logistic Function'),
%subplot(122), plot(ym,fX,'.'), grid, title('Logistic Function through time')
%print -deps 'f3_IS1.eps'       % MNLIV + IS1
%print -deps 'f3_IS2.eps'       % MNLIV + IS2
%print -deps 'f3_IS3.eps'       % MNLIV + IS3
%print -deps 'f3_IS4.eps'       % MNLIV + IS4
%print -deps 'f3_IS1s.eps'      % NLIV + IS1
%print -deps 'f3_IS2s.eps'      % NLIV + IS2
%print -deps 'f3_IS3s.eps'      % NLIV + IS3
%print -deps 'f3_IS4s.eps'      % NLIV + IS4
%subplot(121), plot(q,fX2,'.'), grid, title('Impact of transition variable on logistic Function')
%subplot(122), plot(ym,fX2,'.'), grid, title('Logistic Function through time')
%print -deps 'f3_IS1o.eps'      % MNLIV* (Optimal MNLIV) + IS1  
%print -deps 'f3_IS2o.eps'      % MNLIV* + IS2  
%print -deps 'f3_IS3o.eps'      % MNLIV* + IS3  
%print -deps 'f3_IS4o.eps'      % MNLIV* + IS4
%print -deps 'f3_IS1so.eps'      % NLIV* (Optimal NLIV) + IS1 
%print -deps 'f3_IS2so.eps'      % NLIV* + IS2 
%print -deps 'f3_IS3so.eps'      % NLIV* + IS3 
%print -deps 'f3_IS4so.eps'      % NLIV* + IS4


%**************************************************************************
% 7. SAVE RESULTS
%**************************************************************************
MQO     = [alphaO' zeros(1,size(alphaO,1)+2);  sd_hatO' zeros(1,size(alphaO,1)+2)];
TSLS    = [alphaL' zeros(1,size(alphaL,1)+2);  sd_hatL' zeros(1,size(alphaL,1)+2)];
NKPC = [MQO; TSLS; par1'; sd_hat1'; par2'; sd_hat2'];
csvwrite('STRNKPC',NKPC);