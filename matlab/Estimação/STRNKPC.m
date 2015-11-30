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

% y	m	i	inf_esperada_12m	x	meta	dif_inflacao
base = [
2003	03	26.31631424	9.3	-0.16	4	5.3
2003	04	26.31899996	8.52	-0.39	4	4.52
2003	05	26.30809397	8.18	-0.58	4	4.18
2003	06	26.08627211	7.13	-0.7	4	3.13
2003	07	25.35926648	6.67	-0.75	4	2.67
2003	08	23.49672891	6.33	-0.68	4	2.33
2003	09	21.01967367	6.36	-0.51	4	2.36
2003	10	19.53652471	6.08	-0.31	4	2.08
2003	11	18.31230749	5.91	-0.14	4	1.91
2003	12	16.9094661	5.93	0	4	1.93
2004	01	16.31714166	5.85	0.09	5.5	0.35
2004	02	16.29944199	5.46	0.18	5.5	-0.04
2004	03	16.19170007	5.48	0.29	5.5	-0.02
2004	04	15.9559313	5.54	0.43	5.5	0.04
2004	05	15.77047387	5.82	0.6	5.5	0.32
2004	06	15.79904585	6.19	0.77	5.5	0.69
2004	07	15.77363446	6.1	0.93	5.5	0.6
2004	08	15.86136092	6.26	1.04	5.5	0.76
2004	09	16.09371241	6.18	1.07	5.5	0.68
2004	10	16.40974796	6.24	1.01	5.5	0.74
2004	11	16.95824411	6.05	0.91	5.5	0.55
2004	12	17.50060027	5.83	0.79	5.5	0.33
2005	01	17.93070333	5.61	0.67	4.5	1.11
2005	02	18.47307719	5.45	0.57	4.5	0.95
2005	03	18.97246952	5.46	0.48	4.5	0.96
2005	04	19.31644278	5.59	0.37	4.5	1.09
2005	05	19.60518219	5.21	0.23	4.5	0.71
2005	06	19.7454544	5	0.03	4.5	0.5
2005	07	19.72238084	4.99	-0.19	4.5	0.49
2005	08	19.7482603	4.92	-0.41	4.5	0.42
2005	09	19.60613588	4.69	-0.6	4.5	0.19
2005	10	19.25174474	4.64	-0.7	4.5	0.14
2005	11	18.86529897	4.63	-0.72	4.5	0.13
2005	12	18.24160574	4.62	-0.69	4.5	0.12
2006	01	17.65215817	4.57	-0.67	4.5	0.07
2006	02	17.28055416	4.43	-0.71	4.5	-0.07
2006	03	16.73518857	4.25	-0.82	4.5	-0.25
2006	04	16.19444312	4.13	-0.95	4.5	-0.37
2006	05	15.70135636	4.21	-1.05	4.5	-0.29
2006	06	15.17619029	4.34	-1.08	4.5	-0.16
2006	07	14.97924594	4.43	-1.05	4.5	-0.07
2006	08	14.65821386	4.53	-0.99	4.5	0.03
2006	09	14.17049994	4.1	-0.89	4.5	-0.4
2006	10	13.95069436	4.06	-0.75	4.5	-0.44
2006	11	13.64894844	4.18	-0.55	4.5	-0.32
2006	12	13.18599989	4.06	-0.33	4.5	-0.44
2007	01	13.12949356	4	-0.1	4.5	-0.5
2007	02	12.93222215	3.79	0.1	4.5	-0.71
2007	03	12.73813381	3.7	0.27	4.5	-0.8
2007	04	12.57943374	3.5	0.42	4.5	-1
2007	05	12.43	3.4	0.56	4.5	-1.1
2007	06	12.02982179	3.52	0.69	4.5	-0.98
2007	07	11.72518394	3.65	0.81	4.5	-0.85
2007	08	11.4295652	3.74	0.93	4.5	-0.76
2007	09	11.21943636	3.84	1.03	4.5	-0.66
2007	10	11.17954544	3.83	1.1	4.5	-0.67
2007	11	11.18	3.96	1.12	4.5	-0.54
2007	12	11.18	3.96	1.14	4.5	-0.54
2008	01	11.17818175	4.3	1.19	4.5	-0.2
2008	02	11.17947366	4.26	1.32	4.5	-0.24
2008	03	11.17899996	4.35	1.51	4.5	-0.15
2008	04	11.36977279	4.41	1.75	4.5	-0.09
2008	05	11.62699905	4.79	1.94	4.5	0.29
2008	06	12.08937277	5.29	2.01	4.5	0.79
2008	07	12.36430236	5.37	1.85	4.5	0.87
2008	08	12.91904758	5.26	1.41	4.5	0.76
2008	09	13.38989019	5.16	0.68	4.5	0.66
2008	10	13.66478246	5.32	-0.26	4.5	0.82
2008	11	13.63999978	5.44	-1.29	4.5	0.94
2008	12	13.65999996	4.96	-2.19	4.5	0.46
2009	01	13.32425139	4.69	-2.8	4.5	0.19
2009	02	12.66	4.58	-3.07	4.5	0.08
2009	03	11.70312969	4.21	-3.09	4.5	-0.29
2009	04	11.10928529	4.14	-2.92	4.5	-0.36
2009	05	10.15949998	4.02	-2.62	4.5	-0.48
2009	06	9.539403967	4.09	-2.26	4.5	-0.41
2009	07	9.005840253	4.08	-1.85	4.5	-0.42
2009	08	8.65	4.07	-1.42	4.5	-0.43
2009	09	8.65	4.25	-0.98	4.5	-0.25
2009	10	8.65	4.29	-0.59	4.5	-0.21
2009	11	8.65	4.4	-0.28	4.5	-0.1
2009	12	8.649545435	4.4	-0.04	4.5	-0.1
2010	01	8.65	4.51	0.14	4.5	0.01
2010	02	8.65	4.49	0.29	4.5	-0.01
2010	03	8.65	4.62	0.43	4.5	0.12
2010	04	8.722269664	4.72	0.55	4.5	0.22
2010	05	9.398571373	4.76	0.63	4.5	0.26
2010	06	9.93661689	4.85	0.66	4.5	0.35
2010	07	10.31884537	4.97	0.67	4.5	0.47
2010	08	10.66	5.02	0.68	4.5	0.52
2010	09	10.66	5.14	0.7	4.5	0.64
2010	10	10.66	5.16	0.73	4.5	0.66
2010	11	10.66	5.36	0.77	4.5	0.86
2010	12	10.6613043	5.35	0.8	4.5	0.85
2011	01	10.85306598	5.51	0.81	4.5	1.01
2011	02	11.16949998	5.35	0.82	4.5	0.85
2011	03	11.62228422	5.34	0.82	4.5	0.84
2011	04	11.73573528	5.34	0.81	4.5	0.84
2011	05	11.92	5.04	0.76	4.5	0.54
2011	06	12.0985145	5.18	0.67	4.5	0.68
2011	07	12.2532715	5.4	0.52	4.5	0.9
2011	08	12.41782601	5.45	0.32	4.5	0.95
2011	09	11.90523798	5.71	0.09	4.5	1.21
2011	10	11.69923214	5.62	-0.11	4.5	1.12
2011	11	11.4	5.54	-0.26	4.5	1.04
2011	12	10.90045453	5.48	-0.34	4.5	0.98
2012	01	10.69518142	5.31	-0.38	4.5	0.81
2012	02	10.4	5.3	-0.38	4.5	0.8
2012	03	9.820005896	5.4	-0.37	4.5	0.9
2012	04	9.349382149	5.54	-0.37	4.5	1.04
2012	05	8.870403781	5.5	-0.37	4.5	1
2012	06	8.390499978	5.5	-0.35	4.5	1
2012	07	8.071550755	5.56	-0.31	4.5	1.06
2012	08	7.84642948	5.65	-0.26	4.5	1.15
2012	09	7.39	5.52	-0.21	4.5	1.02
2012	10	7.230841682	5.37	-0.16	4.5	0.87
2012	11	7.14	5.39	-0.14	4.5	0.89
2012	12	7.157485408	5.52	-0.13	4.5	1.02
2013	01	7.11	5.47	-0.11	4.5	0.97
2013	02	7.12388432	5.62	-0.07	4.5	1.12
2013	03	7.147498886	5.43	0.04	4.5	0.93
2013	04	7.259024779	5.58	0.2	4.5	1.08
2013	05	7.4237569	5.67	0.36	4.5	1.17
2013	06	7.9	5.65	0.5	4.5	1.15
2013	07	8.225824709	5.89	0.59	4.5	1.39
2013	08	8.445359524	6.12	0.62	4.5	1.62
2013	09	8.9	6.22	0.61	4.5	1.72
2013	10	9.247583549	6.19	0.59	4.5	1.69
2013	11	9.449897463	6.09	0.57	4.5	1.59
2013	12	9.899523789	6.09	0.55	4.5	1.59
2014	01	10.17244589	6	0.54	4.5	1.5
2014	02	10.42497456	6.12	0.5	4.5	1.62
2014	03	10.65	6.1	0.43	4.5	1.6
2014	04	10.8749746	5.95	0.33	4.5	1.45
2014	05	10.9	6.01	0.21	4.5	1.51
2014	06	10.9	5.9	0.08	4.5	1.4
2014	07	10.9	6.03	0	4.5	1.53
2014	08	10.90142852	6.24	0	4.5	1.74
2014	09	10.9	6.33	0.01	4.5	1.83
2014	10	10.92171679	6.38	-0.02	4.5	1.88
2014	11	11.15	6.57	-0.15	4.5	2.07
2014	12	11.58168596	6.61	-0.41	4.5	2.11
2015	01	11.81641849	6.61	-0.83	4.5	2.11
2015	02	12.15	6.54	-1.37	4.5	2.04
2015	03	12.58168714	6.17	-1.95	4.5	1.67
];

selic = [
2003	03	26.32
2003	04	26.32
2003	05	26.31
2003	06	26.09
2003	07	25.36
2003	08	23.50
2003	09	21.02
2003	10	19.54
2003	11	18.31
2003	12	16.91
2004	01	16.32
2004	02	16.30
2004	03	16.19
2004	04	15.96
2004	05	15.77
2004	06	15.80
2004	07	15.77
2004	08	15.86
2004	09	16.09
2004	10	16.41
2004	11	16.96
2004	12	17.50
2005	01	17.93
2005	02	18.47
2005	03	18.97
2005	04	19.32
2005	05	19.61
2005	06	19.75
2005	07	19.72
2005	08	19.75
2005	09	19.61
2005	10	19.25
2005	11	18.87
2005	12	18.24
2006	01	17.65
2006	02	17.28
2006	03	16.74
2006	04	16.19
2006	05	15.70
2006	06	15.18
2006	07	14.98
2006	08	14.66
2006	09	14.17
2006	10	13.95
2006	11	13.65
2006	12	13.19
2007	01	13.13
2007	02	12.93
2007	03	12.74
2007	04	12.58
2007	05	12.43
2007	06	12.03
2007	07	11.73
2007	08	11.43
2007	09	11.22
2007	10	11.18
2007	11	11.18
2007	12	11.18
2008	01	11.18
2008	02	11.18
2008	03	11.18
2008	04	11.37
2008	05	11.63
2008	06	12.09
2008	07	12.36
2008	08	12.92
2008	09	13.39
2008	10	13.66
2008	11	13.64
2008	12	13.66
2009	01	13.32
2009	02	12.66
2009	03	11.70
2009	04	11.11
2009	05	10.16
2009	06	9.54
2009	07	9.01
2009	08	8.65
2009	09	8.65
2009	10	8.65
2009	11	8.65
2009	12	8.65
2010	01	8.65
2010	02	8.65
2010	03	8.65
2010	04	8.72
2010	05	9.40
2010	06	9.94
2010	07	10.32
2010	08	10.66
2010	09	10.66
2010	10	10.66
2010	11	10.66
2010	12	10.66
2011	01	10.85
2011	02	11.17
2011	03	11.62
2011	04	11.74
2011	05	11.92
2011	06	12.10
2011	07	12.25
2011	08	12.42
2011	09	11.91
2011	10	11.70
2011	11	11.40
2011	12	10.90
2012	01	10.70
2012	02	10.40
2012	03	9.82
2012	04	9.35
2012	05	8.87
2012	06	8.39
2012	07	8.07
2012	08	7.85
2012	09	7.39
2012	10	7.23
2012	11	7.14
2012	12	7.16
2013	01	7.11
2013	02	7.12
2013	03	7.15
2013	04	7.26
2013	05	7.42
2013	06	7.90
2013	07	8.23
2013	08	8.45
2013	09	8.90
2013	10	9.25
2013	11	9.45
2013	12	9.90
2014	01	10.17
2014	02	10.42
2014	03	10.65
2014	04	10.87
2014	05	10.90
2014	06	10.90
2014	07	10.90
2014	08	10.90
2014	09	10.90
2014	10	10.92
2014	11	11.15
2014	12	11.58
2015	01	11.82
2015	02	12.15
2015	03	12.58
];


%**************************************************************************
% PREPARACAO DOS DADOS
%**************************************************************************
[T,s]  = size(base);
%li    = log(1+(selic(:,3)/100)); % log-anualizado

% dados
% y      : regressando
% xex    : regressores exogenos
% xen    : regressores endogenos
% z      : variavel de transicao (funcao logistica)
% wx     : instrumentos para x
% wz     : instrumentos para z
% ym     : periodos
% --------------
% contemporaneo
y   = selic(5:T,3);
xex = [selic(4:T-1,3) base(4:T-1,5) base(4:T-1,7)];
xen = [base(5:T,5) base(5:T,7)]; % hiato do produto e diferencial de inflacao
x   = [xex xen];
z   = xen;
wx  = [base(4:T-1,5) base(3:T-2,5) base(2:T-3,5) base(1:T-4,5) base(4:T-1,7) base(3:T-2,7) base(2:T-3,7) base(1:T-4,7) selic(4:T-1,3) selic(3:T-2,3) selic(2:T-3,3) selic(1:T-4,3)]; % uso como instrumento proprias variaveis defasadas
wz  = [base(4:T-1,5) base(3:T-2,5) base(2:T-3,5) base(1:T-4,5) base(4:T-1,7) base(3:T-2,7) base(2:T-3,7) base(1:T-4,7)]; % uso como instrumento proprias variaveis defasadas
ym  = base(5:T,1)+(base(5:T,2)/12);


% teste ADF, selic com tend�ncia linear
[h,pValue,stat,cValue,reg] = adftest(y,'model','TS','lags',1);
[h,pValue,stat,cValue,reg] = adftest(xen(:,1),'lags',1);
[h,pValue,stat,cValue,reg] = adftest(xen(:,1),'lags',1);


% *************************************************************************
% 3. TESTES DE N�O LINEARIDADE
% *************************************************************************
% funcao LinTest.m
% ------------------------

% PRIMEIRO TESTE: VARI�VEIS SEM LAG - TESTE CONTEMPORANEO
% primeiro teste: hiato do produto
[F, f, isH0, pvalue] = LinTest(y,xex,xen,z(:,1),wx,wz(:,1),flag,sig,instr_Amemya);

% segundo teste: desvios da infla��o
[F, f, isH0, pvalue] = LinTest(y,xex,xen,z(:,2),wx,wz(:,2),flag,sig,instr_Amemya);



% SEGUNDO TESTE: LAG = 1
%xen = [base(4:T-1,5) base(4:T-1,7)]; % hiato do produto e diferencial de inflacao
%x   = [xex xen];
%z   = xen;

% primeiro teste: hiato do produto
%[F, f, isH0, pvalue] = LinTest(y,xex,xen,z(:,1),wx,wz(:,1),flag,sig,instr_Amemya);

% segundo teste: desvios da infla��o
%[F, f, isH0, pvalue] = LinTest(y,xex,xen,z(:,2),wx,wz(:,2),flag,sig,instr_Amemya);



% TERCEIRO TESTE: LAG = 2
%xen = [base(3:T-2,5) base(3:T-2,7)]; % hiato do produto e diferencial de inflacao
%x   = [xex xen];
%z   = xen;

% primeiro teste: hiato do produto
%[F, f, isH0, pvalue] = LinTest(y,xex,xen,z(:,1),wx,wz(:,1),flag,sig,instr_Amemya);

% segundo teste: desvios da infla��o
%[F, f, isH0, pvalue] = LinTest(y,xex,xen,z(:,2),wx,wz(:,2),flag,sig,instr_Amemya);



% QUARTO TESTE: LAG = 3
%xen = [base(2:T-3,5) base(2:T-3,7)]; % hiato do produto e diferencial de inflacao
%x   = [xex xen];
%z   = xen;

% primeiro teste: hiato do produto
%[F, f, isH0, pvalue] = LinTest(y,xex,xen,z(:,1),wx,wz(:,1),flag,sig,instr_Amemya);

% segundo teste: desvios da infla��o
%[F, f, isH0, pvalue] = LinTest(y,xex,xen,z(:,2),wx,wz(:,2),flag,sig,instr_Amemya);


% QUINTO TESTE: LAG = 4
%xen = [base(1:T-4,5) base(1:T-4,7)]; % hiato do produto e diferencial de inflacao
%x   = [xex xen];
%z   = xen;

% primeiro teste: hiato do produto
%[F, f, isH0, pvalue] = LinTest(y,xex,xen,z(:,1),wx,wz(:,1),flag,sig,instr_Amemya);

% segundo teste: desvios da infla��o
%[F, f, isH0, pvalue] = LinTest(y,xex,xen,z(:,2),wx,wz(:,2),flag,sig,instr_Amemya);



% *************************************************************************
% 4. ESTIMACAO
% *************************************************************************
% dados
% ------------
y   = selic(5:T,3);
xex = [selic(4:T-1,3) base(4:T-1,5) base(4:T-1,7)];
xen = [base(5:T,5) base(5:T,7)]; % hiato do produto e diferencial de inflacao
x   = [xex xen];
z   = xen;
wx  = [base(4:T-1,5) base(3:T-2,5) base(2:T-3,5) base(1:T-4,5) base(4:T-1,7) base(3:T-2,7) base(2:T-3,7) base(1:T-4,7) selic(4:T-1,3) selic(3:T-2,3) selic(2:T-3,3) selic(1:T-4,3)]; % uso como instrumento proprias variaveis defasadas
wz  = [base(4:T-1,5) base(3:T-2,5) base(2:T-3,5) base(1:T-4,5) base(4:T-1,7) base(3:T-2,7) base(2:T-3,7) base(1:T-4,7)]; % uso como instrumento proprias variaveis defasadas
ym  = base(5:T,1)+(base(5:T,2)/12);

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
disp('#m #par J0 JA isH0J pvalueJ');
[size(W'*ehatL,1) size(X,2) J0 JA isH0J pvalueJ];


% -------------------------------------------------------------------------
% 4.2. ESTIMATIVAS NAO LINEARES
% -------------------------------------------------------------------------
%

model = eye(2*nX+2);
ord = model*(1:2*nX+2)';
PV = ones(2*nX+2,1);

% estima��o
[par1, sd_hat1, VM1, ehat1, aic1, fX1, fX2] = Estimate(y,xex,xen,z,W,PV,ord,model,0);
%[J0,JA,isH0J,pvalueJ,F1,f1,isH01,pvalue1,F2,f2,isH02,pvalue2] = diagnostics(par1,X,y,z(:,2),W,VM1);

% primeiro teste: hiato do produto
[F, f, isH0, pvalue] = LinTest(y,xex,xen,z(:,1),wx,wz(:,1),flag,sig,instr_Amemya);

disp('Lineariy Test - hiato do produto [F f isH0 pvalue]:')
[F  f  isH0  pvalue]

% segundo teste: desvios da infla��o
[F, f, isH0, pvalue] = LinTest(y,xex,xen,z(:,2),wx,wz(:,2),flag,sig,instr_Amemya);

disp('Lineariy Test - desvios da infla��o [F f isH0 pvalue]:')
[F  f  isH0  pvalue]

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

disp('Sum squared errors')
sumsqr(ehat1)

disp('Akaike Information Criteria')
sumsqr(aic1)

figure('Name','Autocorrela��o dos res�duos')
autocorr(ehat1)

figure('Name','Fun��o log�stica')
subplot(221), plot(q(:,1),fX1,'.k'), grid, title('Impacto do hiato do produto na fun��o log�stica'),
subplot(222), plot(ym,fX1,'.k'), grid, title('Fun��o log�stica referente ao hiato do produto no tempo')
subplot(223), plot(q(:,2),fX2,'.k'), grid, title('Impacto dos desvios da infla��o na fun��o log�stica'),
subplot(224), plot(ym,fX2,'.k'), grid, title('Fun��o log�stica referente aos desvios da infla��o no tempo')
