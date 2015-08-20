% *************************************************************************
% main
% *************************************************************************

% Inputs: 
% -------
% Function Data.m:
% y      : regressando
% xex    : regressores exogenos
% xen    : regressores endogenos
% z      : variavel de transicao (funcao logistica)
% wx     : instrumentos para x
% wz     : instrumentos para z
% ym     : periodos
%
%
% Dependencias:
% -------------------
% 1. Data.m          : dados
% 2. LinTest.m       : teste de linearidade
% 3. Estimate.m      : estimacao
%    3.1. IVCost.m   : criterio de minimizacao
%
%
% *************************************************************************
% 1. PARAMETROS
% *************************************************************************
flag         = 1;                      % com intercepto 
                                       
sig          = 0.95;                   % nivel de significancia

instr_Amemya = 0;                      % =0 usa instrumentos padrao
                                       % <>0 instrumentos otimos (Amemya)


% *************************************************************************
% 2. DADOS
% *************************************************************************

[y,xex,xen,z,wx,wz,ym] = Data;

% teste ADF, selic com tendência linear
[h,pValue,stat,cValue,reg] = adftest(y,'model','TS','lags',1);
[h,pValue,stat,cValue,reg] = adftest(xen(:,1),'lags',1);
[h,pValue,stat,cValue,reg] = adftest(xen(:,1),'lags',1);


% *************************************************************************
% 3. TESTE DE NÃO LINEARIDADE
% *************************************************************************
% funcao LinTest.m
% ------------------------
% primeiro teste: hiato do produto
[F, f, isH0, pvalue] = LinTest(y,xex,xen,z(:,1),wx,wz(:,1),flag,sig,instr_Amemya)

% segundo teste: desvios da inflação
[F, f, isH0, pvalue] = LinTest(y,xex,xen,z(:,2),wx,wz(:,1),flag,sig,instr_Amemya)


% *************************************************************************
% 4. ESTIMACAO
% *************************************************************************
% dados
% ------------
[T,neq] = size(y);                                
    
Y       = y;                          
X       = [ones(T,1) xex xen];         
Z       = [ones(T,1) wx wz];          
q       = z;                                     
qZ      = wz(:,1);                     

nxex    = size(xex,2)+1;              
nxen    = size(xen,2);                 
nX      = size(X,2);                  
nZ      = size(Z,2);                 

% -------------------------------------------------------------------------
% 4.1. LINEAR (OLS & 2SLS)
% -------------------------------------------------------------------------
W   = Z;                       
                                
nW  = size(W,2);                   
PW  = W*pinv(W'*W)*W';              

% 4.1.1. OLS
% ----------
alphaO   = inv(X'*X)*X'*Y;         
YhatO    = X*alphaO;               
ehatO    = y - YhatO;             
v_ehatO  = var(ehatO);              
sd_hatO  = sqrt(diag(v_ehatO*inv(X'*X)));

% 4.1.2. 2SLS
% -----------
alphaL   = inv(X'*PW*X)*X'*PW*Y;    
YhatL    = X*alphaL;                
ehatL    = y - YhatL;            
v_ehatL  = var(ehatL);              
sd_hatL  = sqrt(diag(v_ehatL*inv(X'*PW*X)));

% 4.1.3. TESTE J
% ---------------
m = W'*ehatL/T;
M = -W'*X/T;
S = v_ehatL*(W'*W)/T;


% 
% --------------------
J0   = chi2inv(sig,size(W'*ehatL,1)- size(X,2));

% modelo sobre Ha
% ---------------------------
JA = T*m'*pinv(S)*m;

% estatisticas
% --------------------------------
J0 = J0';
isH0J = JA <= J0;
pvalueJ=1-chi2cdf(JA,size(W'*ehatL,1)- size(X,2));
disp('#m #par J0 JA isH0J pvalueJ')
[size(W'*ehatL,1) size(X,2) J0 JA isH0J pvalueJ]


% -------------------------------------------------------------------------
% 4.2. ESTIMATIVAS NAO LINEARES
% -------------------------------------------------------------------------
%