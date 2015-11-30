function [PSI, sd_hat, VM, ehat, aic, fX1, fX2] = Estimate(y,xex,xen,z,W,PV,ord,model,par0)
% *************************************************************************
% FUNCTION Estimate
% *************************************************************************
% Last update: July 22, 2008
% By: Waldyr Areosa; Marcelo C. Medeiros
% *************************************************************************
% Function for
% INSTRUMENTAL VARIABLE ESTIMATION OF SMOOTH TRANSITION REGRESSION MODELS
% WITH ENDOGENOUS VARIABLES
% By
% Waldyr Areosa; Michael McAleer; Marcelo C. Medeiros
% *************************************************************************
% Descripition
% Estimates STR Model
% *************************************************************************
% Inputs:
% -------
% y      : Regressand
% xex    : Exogenous Regressors
% xen    : Endogenous Regressors
% z      : Logistic Function (Transition Variable)
% W      : Instruments
% PV     : True parameter values
% ord    : List the order of the parameters to be estimated
% model  : Matriz that reduce the dimension of the original model
%
% Outputs:
% --------
% par_hat : Parameters estimates
% sd_hat  : Standart Erros
%
% Functions required:
% -------------------
% 1. parestlm.m : Starting values
% 2. IVCost     : Minimization Criteria
%
%
%**************************************************************************
% 1. READ AND PREPARE DATA
%**************************************************************************
T = size(y,1);
x = [xex xen];
X   = [ones(T,1) x];

PW = W*pinv(W'*W)*W';
nX = size(X,2);                     % # regressors
nW = size(W,2);                     % # instruments


%**************************************************************************
% 2. NONLINEAR ESTIMATES (NLIV or MNLIV)
%**************************************************************************
% STEP 1: Start Values
% --------------------
%[alpha0,lambda0,gamma0,c0] = parestlm(y,X,z,nX,[],[]);
gamma10  = 20;
gamma20  = 20;
z1 = z(:,1); % hiato do produto
z2 = z(:,2); % desvios da inflacao
%c0      = median(z);  
%fX      = (1./(1+exp(-gamma0.*(z - c0)))) ;
fX1      = (1./(1+exp(-gamma10.*(z1 - 0))));
fX2     = (1./(1+exp(gamma20.*(z2 - 0))));
%fW      = (1./(1+exp(-gamma0.*(W(:,end) - c0))));
fW1      = (1./(1+exp(-gamma10.*(W(:,end) - 0))));
fW2     = (1./(1+exp(gamma20.*(W(:,end) - 0))));
fX      = fX1 + fX2;
fW      = fW1 + fW2;
X2      = [X repmat(fX,1,nX).*X]; 
W2      = [W repmat(fW,1,nW).*W];
PW2      = W2*pinv(W2'*W2)*W2';
theta0  = pinv(X2'*PW2*X2)*X2'*PW2*y;
alpha0  = theta0(1:nX);
lambda0 = reshape(theta0(nX+1:end),nX,1);

% STEP 2: Estimation
% ------------------
if par0==0
    %PSI     = [alpha0;lambda0;gamma0;c0];
    PSI     = [alpha0;lambda0;gamma10;gamma20];
    PSI     = model*PSI;            % Reduced set of parameters to be estimated
else
    PSI     = model*par0;
end

func    = 'IVCost';
%options = optimset('Display','off','GradObj','on','MaxFunEvals',1e10,'LargeScale','off','MaxIter',4000,'TolFun',1e-8,'DerivativeCheck','off','LevenbergMarquardt','on','LineSearchType','cubicpoly','TolX',1e-8);
options = optimset('Display','off','Jacobian','on','MaxFunEvals',1e10,'LargeScale','off','MaxIter',6000,'TolFun',1e-8,'DerivativeCheck','off','LevenbergMarquardt','on','LineSearchType','cubicpoly','TolX',1e-8);

LB0     = -inf*ones(2*nX,1);
UB0     = inf*ones(2*nX,1);
%LB2     = [LB0;1;prctile(z2,15)];
LB2     = [LB0;1;-30];
%UB2     = [UB0;60;prctile(z2,85)];
UB2     = [UB0;60;+30];

for i = 1:length(ord)
    LB(i) = LB2(ord(i));
    UB(i) = UB2(ord(i));    
end

PSI     = fmincon(func,PSI,[],[],[],[],LB,UB,[],options,X,y,z1,z2,PW,PV,ord,model);

par_hat = PV;
for i = 1:length(ord)
    par_hat(ord(i))=PSI(i);
end

alphaf  = par_hat(1:nX);
lambdaf = par_hat(nX+1:2*nX);
gammaf1  = par_hat(2*nX+1);
gammaf2  = par_hat(2*nX+2);
%cf      = par_hat(2*nX+2);


%**************************************************************************
% 3. CALCULATE STANDART DEVIATIONS (NLIV OR MNLIV)
%**************************************************************************
%fX      = (1./(1+exp(-gammaf.*(z - cf))));
fX      = (1./(1+exp(-gammaf1.*(z1 - 0)))) + (1./(1+exp(gammaf2.*(z2 - 0))));
Yhat    = X*alphaf + (X*lambdaf).*fX;
ehat    = y - Yhat;
Ghat    = [X X.*repmat(fX,1,nX) (X*lambdaf).*fX.*(1-fX).*(z1 - 0) (X*lambdaf).*fX.*(1-fX).*(z2 - 0) - gammaf1*(X*lambdaf).*fX.*(1-fX) + gammaf2*(X*lambdaf).*fX.*(1-fX)];
Ghat    = Ghat*model';
V       = [xen z]  - PW*[xen z];
XZ      = [ehat V];
SIGMA   = cov(XZ);
S11     = SIGMA(1,1);
S22     = SIGMA(2:end,2:end);
S12     = SIGMA(1,2:end);
S21     = SIGMA(2:end,1);
S11_star= S11 - S12*pinv(S22)*S21;
PV      = V*pinv(V'*V)*V';
G       = Ghat'*(eye(size(PV))-PV)*Ghat;
H       = Ghat'*PW*Ghat;
VM      = pinv(G)*(S11_star*G+(S11-S11_star)*H)*pinv(G);
sd_hat  = sqrt(diag(VM)); 
aic     = 2*(size(par_hat(1:nX,1),1)+2)+T*log(sumsqr(ehat)/T);