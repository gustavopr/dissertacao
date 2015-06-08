function [f] = IVCost(PSI,X,y,z,PW,PV,ord,model)
%function [f,g] = IVCost(PSI,X,y,z,PW,PV,ord,model)
% *************************************************************************
% FUNCTION IVCost
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
% Minimization criteria for STR model
% *************************************************************************
% Inputs:
% -------
% PSI    : Parameter Values
% X      : Regressors 
% y:     : Regressand
% z      : Logistic Function (Transition Variable)
% PW     : Instruments' Projection Matrix
% PV     : True parameter values
% ord    : List the order of the parameters to be estimated
% model  : Matrix that reduce the dimension of the original model

% Outputs:
% --------
% f      : Weighted Squared Sum of Errors
% PSI    : Parameter Values
%
% Functions required: None
% -------------------
%
%
%**************************************************************************
% 1. READ AND PREPARE DATA
%**************************************************************************
nX = size(X,2);

param = PV;
for i=1:length(ord)
    param(ord(i))=PSI(i);
end

alpha  = param(1:nX);
lambda = param(nX+1:2*nX);
gamma  = param(2*nX+1);
c      = param(2*nX+2);


%**************************************************************************
% 2. CALCULATE f (NLIV or MNLIV)
%**************************************************************************
fX   = (1./(1+exp(-gamma.*(z - c))));
Yhat = X*alpha +(X*lambda).*fX;
e = y - Yhat;

% NLIV
% -----
%PZ = PW;

% MNLIV
% ------
M = eye(size(PW)) - PW;
PZ = eye(size(PW)) - M*X*pinv(X'*M*X)*X'*M;

f = e'*PZ*e;


%**************************************************************************
% 3. CALCULATE g (NLIV or MNLIV)
%**************************************************************************
%G = -[X X.*repmat(fX,1,nX) (X*lambda).*fX.*(1-fX).*(z - c) -gamma*(X*lambda).*fX.*(1-fX)];
%G = G*model';

%g = 2*G'*PZ*e;