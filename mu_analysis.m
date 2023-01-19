%% Nominal plant and controller
global rp w_rp
close all
load('dataset');
s = tf('s');
sys = log_vars.sys;
S = log_vars.S;
T = log_vars.T;
K = log_vars.K;
J = get_linearization();
A_i = J.A_i;
B_i = J.B_i;
C = J.C;
D = J.D;
% autovalori = sigma(T);
% sigma(T)
omega = logspace(-1,6,302);
%% Weighting filter for uncertainty modelling
%wi = 10^2*1/(s)^5/(1+10^5*s)*(1+s*10^2)^6;
Gp = C*(s*eye(5)-A_i)^(-1)*B_i; %G incerta
G_inv = inv(sys'*sys)*sys';
sigma(G_inv*(Gp-sys));
%[freq,response] = ginput(5);
freq = [1.40494020600125e-14
6.69677303868978e-07
0.249740006959577
2368.13236993714
67387776523264.6];

response = [849.322493224932
504.607046070460
10.2981029810296
-93.7669376693768
-76.4227642276426];
system = frd(10.^(response/20),freq');
wi = fitmagfrd(system,4,0);
wi = tf(wi);
rp_tau = w_rp/(rp);
%wi = rp_tau*rp*s/(1+rp*s);
Wi = blkdiag(wi,wi);
WP = log_vars.WP;
WU = log_vars.WU;

%% Generalized plant P with Wi, Wu and Wp

% systemnames = 'sys WP';
% inputvar = '[w(4) ; u(2)]';
% outputvar = '[WP; -sys-w]';
% input_to_sys= '[u]';
% input_to_WP = '[sys+w]'; 
% sysoutname = 'P'; cleanupsysic = 'yes';
% sysic;
% 
% nmeas=4;
% nu=2;
% [K,CL,gamma] = hinfsyn(P,nmeas,nu);


systemnames = 'sys WP WU Wi';
inputvar = '[udel{2}; w{4}; u{2}]';
outputvar = '[Wi ;WU; WP ; -w-sys]';
input_to_sys = '[u+udel]';
input_to_WP = '[sys]';
input_to_WU = '[u]';
input_to_Wi = '[u]';
sysoutname = 'P';
cleanupsysic = 'yes';
sysic;
P = minreal(ss(P));
% 
% nmeas=4;
% nu=2;
% [K,CL,gamma] = hinfsyn(P,nmeas,nu);

%% MDelta system
N = lft(P,K);
Nf = frd(N,omega);

% Matrix N
Delta1 = ultidyn('Delta1',[1 1]);
Delta2 = ultidyn('Delta2',[1 1]);
Delta = blkdiag(Delta1,Delta2);
Gpp = sys*(eye(2)+Wi*Delta);
sigma(Gpp, 'r'); hold on; sigma(sys, 'b');
% G_pinv = inv(sys'*sys)*sys';
% bodemag(usample(G_pinv*(Gp-Gp.NominalValue),100))
% sigma(G_pinv*(Gp-sys));
% hold on;
% sigma(Wi);

%bodemag(sys_inv*tf(Gp-sys),'r'); hold on; bodemag(Wi,'b');
%bodemag(usample(sys_inv*(Gp-sys),50),'r'); hold on; bodemag(Wi,'b');
M = lft(Delta,N);
Mf = frd(M,omega);

%% RS with mussv, rea

%osservo autovalori di N per capire se è NS
eig(N)
% M = N(1,1), per la robusta stabilità la norma infinito di M deve essere
% minore di 1
Nrs = Nf(1:2,1:2);
[mubnds,muinfo] = mussv(Nrs,[1 1; 1 1],'a');
muRS = mubnds(:,1);
[muRSinf,muRSw] = norm(muRS,inf);

%per la Performance Nominale devo fare un controllo sulla N22
% DELTAP è una matrice complessa piena delle stesse dimensioni
%di F' dove F è lft(N,Delta)=M;
%Quindi essendo M una matrice 6x4 M' è 4x6
Nnp=Nf(3:8,3:6); % Picking out wP*Si
[mubnds,muinfo]=mussv(Nnp,[4 6],'a');
muNP = mubnds(:,1);
[muNPinf,muNSw]=norm(muNP,inf);

%Performance robusta
[mubnds,muinfo]=mussv(Nf,[1 1;1 1;4 6],'a');
muRP = mubnds(:,1);
[muRPinf,muNSw]=norm(muRP,inf);

%% plots

figure(1);
sigma(mubnds,'r-'); hold on; sigma(1/Wi,'g'); 
bodemag(mubnds(:,2),'r-'); hold on; bodemag(1/Wi(2,2),'g'); hold on; bodemag(frd(autovalori(1,:),omega),'b')
% figure(3);
%bodemag(mubnds,'r-'); hold on; bodemag(1/Wi,'g');
%legend('mu(T)', '1/|wp|', '$$\bar{\sigma}(T)$$', 'Interpreter','latex')

