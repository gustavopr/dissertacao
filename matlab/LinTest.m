function [F,f,isH0,pvalue,e0,eA] = LinTest(y,xex,xen,z,wx,wz,flag,sig,instr_Amemya)
% *************************************************************************
% FUNCTION LinTest
% *************************************************************************
% Last update: July 16, 2008
% By: Waldyr Areosa; Marcelo C. Medeiros
% *************************************************************************
% Function for
% INSTRUMENTAL VARIABLE ESTIMATION OF SMOOTH TRANSITION REGRESSION MODELS
% WITH ENDOGENOUS VARIABLES
% By
% Waldyr Areosa; Michael McAleer; Marcelo C. Medeiros
% *************************************************************************
% Descripition
% Linearity test for STR model
% 
%    y = X*B1 + (X*B2).*f + u;
%    f = (1./(1+exp(-GAMMA*(z - C))));   
% 
% *************************************************************************
% Inputs:
% -------
% y            : Regressand
% xex          : Exogenous Regressors
% xed          : Endogenous Regressors
% z            : Logistic Function (Transition Variable)
% wx           : Instruments for x
% wz           : Instruments for z
% flag         : =1 if z belongs to x (do not drop intercept)
%                =0 otherwise (drop intercept)
% sig          : Significance level
% instr_Amemya : =1 uses Optimal Instruments (Amemya)
%                <>1 uses Standart Instruments
% Outputs:
% --------
% F            : Test Statistic
% f            : Asymptotic Critical Value
% isH0         : 1 if F <= f
%                0 cc 
% pvalue       : pvalue of F
% e0           : Linear Error
% e2           : Nonlinear Error
%
% Functions required: None
% -------------------


%**************************************************************************
% 1. READ AND PREPARE DATA
%**************************************************************************
x = [xex xen];
[T,N] = size(x);
Nw = size(wx,2);

X   = [ones(T,1) x];

if flag == 1
    XA = [ones(T,1) x x.*repmat(z,1,N) x.*repmat(z,1,N).^2 x.*repmat(z,1,N).^3];
    W  = [ones(T,1) wx wx.*repmat(wz,1,Nw) wx.*(repmat(wz,1,Nw).^2) wx.*(repmat(wz,1,Nw).^3)];
else
    XA = [X X.*repmat(z,1,N+1) X.*repmat(z,1,N+1).^2 X.*repmat(z,1,N+1).^3];
    wxaux = [ones(T,1) wx];
    W  = [wxaux wxaux.*repmat(wz,1,Nw+1) wxaux.*(repmat(wz,1,Nw+1).^2) wxaux.*(repmat(wz,1,Nw+1).^3)];
end
Pw = W*inv(W'*W)*W';


%**************************************************************************
% 2. ESTIMATE HYPOTHESIS
%**************************************************************************
if instr_Amemya == 0 
% -------------------------------------------------------------------------
% 2.1. Standart Instruments
% -------------------------------------------------------------------------

    % -------------------------------------------------------------------------
    % 2.1.1. Model under the null
    % -------------------------------------------------------------------------
    b0 = inv(X'*Pw*X)*X'*Pw*y;
    e0 = y - X*b0;
    SSEr = e0'*Pw*e0;

    % -------------------------------------------------------------------------
    % 2.1.2. Model under the alternative
    % -------------------------------------------------------------------------
    bA = inv(XA'*Pw*XA)*XA'*Pw*y;
    eA = y - XA*bA;
    SSEu = eA'*Pw*eA;

else
% -------------------------------------------------------------------------
% 2.2. Amemya Instruments
% -------------------------------------------------------------------------
    % -------------------------------------------------------------------------
    % 2.2.1. Model under the null
    % -------------------------------------------------------------------------
    M   = eye(size(Pw)) - Pw;
    Pz  = eye(size(Pw)) - M*XA*pinv(XA'*M*XA)*XA'*M;
    b0  = inv(X'*Pz*X)*X'*Pz*y;
    e0  = y - X*b0;
    SSEr  = e0'*Pz*e0;

    % -------------------------------------------------------------------------
    % 2.2.2. Model under the alternative
    % -------------------------------------------------------------------------
    bA = inv(XA'*Pz*XA)*XA'*Pz*y;
    eA = y - XA*bA;
    SSEu = eA'*Pz*eA;

end


%**************************************************************************
% 3. CALCULATE AND COMPARE STATISTICS
%**************************************************************************
k = size(XA,2) - size(X,2);
kk = size(XA,2);
    
F = ((SSEr - SSEu)/(k))/(sumsqr(eA)/(T-kk));
F  = ones(length(sig),1)*F;

f = finv(sig,k,T-kk);           
f = f';

isH0 = F <= f;
pvalue = 1-fcdf(F,k,T-kk);