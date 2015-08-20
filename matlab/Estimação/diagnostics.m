function [J0,JA,isH0J,pvalueJ,F1,f1,isH01,pvalue1,F2,f2,isH02,pvalue2] = diagnostics(par,X,y,q,W,VM)
% diagnostics.m
% Model Evaluation
% Last update: December 01, 2008
% By: Marcelo Medeiros & Waldyr Areosa
%
% inputs:
% ------
% par    : Parameter Values
% X      : Regressors 
% y:     : Regressand
% q      : Logistic Function (Transition Variable)
% PW     : Instruments' Projection Matrix
% VM     : Covariance Matrix
%
% outputs:
% -------
% F: F statistic
% f: theoretical 
% isH0: = F <= f
% pvalue


%**************************************************************************
% 1. READ AND PREPARE DATA
%**************************************************************************
sig = 0.95;
[T, nX] = size(X);

alpha  = par(1:nX);
lambda = par(nX+1:2*nX);
gamma  = par(2*nX+1);
c      = par(2*nX+2);

fX   = (1./(1+exp(-gamma.*(q - c))));
Yhat = X*alpha +(X*lambda).*fX;
u = y - Yhat;
Z = W;
PW  = W*pinv(W'*W)*W';

% NLIV
% -------
PZ = PW;

% MNLIV
% -------
%M = eye(size(PW)) - PW;
%PZ = eye(size(PW)) - M*X*pinv(X'*M*X)*X'*M;

f = u'*PZ*u;

Gd = -[X X.*repmat(fX,1,nX) (X*lambda).*fX.*(1-fX).*(q - c) -gamma*(X*lambda).*fX.*(1-fX)];

%**************************************************************************
% 2. TEST J
%**************************************************************************
%m = PZ*u/T;
m = Z'*u/T;
%M = PZ'*Gd/T;
M = Z'*Gd/T;
S = m*m'/(T - size(par,1));

%cov_m = [eye(size(M,1)) - M*pinv(M'*M)*M']*S*[eye(size(S)) - M*pinv(M'*M)*M']'/T;

V = pinv(M'*pinv(S)*M);

cov_m = (S - M*V*M')/T;

% Model under the null
% --------------------
J0 = chi2inv(1-sig,size(m,1)- size(par,1));

% Model under the alternative
% ---------------------------
JA = T*m'*pinv(cov_m)*m;


% Calculate and Compare Statistics
% --------------------------------
J0 = J0';
isH0J = JA <= J0;
pvalueJ=1-chi2cdf(JA,size(m,1)- size(par,1));    


%**************************************************************************
% 3. TEST 1
%**************************************************************************
% Model under the null
% --------------------
SSE10 = u'*u;

% Model under the alternative
% ---------------------------
X1A = PZ*Gd;
b1A = inv(X1A'*X1A)*X1A'*u;
e1A = u - X1A*b1A;
SSE1A = e1A'*e1A;

% Calculate and Compare Statistics
% --------------------------------
k1 = size(X1A,2);
kk1 = k1;

F1 = ((SSE10 - SSE1A)/k1)/(SSE1A/(T-kk1));
f1 = finv(1-sig,k1,T-kk1);           
f1 = f1';
isH01 = F1 <= f1;
pvalue1=1-fcdf(F1,k1,T-kk1);    


%**************************************************************************
% 4. TEST 2
%**************************************************************************
% Model under the null
% --------------------
X20 = X1A;
b20 = b1A;
e20 = e1A;
SSE20 = SSE1A;

% Model under the alternative
% ---------------------------
X2A = PZ*[Gd X.*repmat(q,1,nX) X.*repmat(q,1,nX).^2 X.*repmat(q,1,nX).^3];
b2A = pinv(X2A'*X2A)*X2A'*e1A;
e2A = e1A - X2A*b2A;
SSE2A = e2A'*e2A;

% Calculate and Compare Statistics
% --------------------------------
k2 = size(X2A,2) - size(X20,2);
kk2 = size(X2A,2);
    
F2 = ((SSE20 - SSE2A)/k2)/(SSE2A/(T-kk2));
f2 = finv(1-sig,k2,T-kk2);           
f2 = f2';
isH02 = F2 <= f2;
pvalue2=1-fcdf(F2,k2,T-kk2);    
