%% Nominal plant and controller
close all
load('dataset');
s = tf('s');
sys = log_vars.sys;
S = log_vars.S;
T = log_vars.T;
K = log_vars.K;
J = get_linearization_lqg();
A_i = J.A_i;
B_i = J.B_i;
C = J.C;
D = J.D;
% autovalori = sigma(T);
% sigma(T)
omega = logspace(-1,6,302);
%% Weighting filter for uncertainty modelling
Wi = log_vars.Wi;
WP = log_vars.WP;
WU = log_vars.WU;

%% Generalized plant P
systemnames = 'sys WP WU Wi';
inputvar = '[udel{2}; ref{4}; u{2}]';
outputvar = '[Wi ; WP ; WU ; ref-sys]';
input_to_sys = '[u+udel]';
input_to_WP = '[sys]';
input_to_WU = '[u]';
input_to_Wi = '[u]';
sysoutname = 'P';
cleanupsysic = 'yes';
sysic;


%% MDelta system
N = lft(P,K);
Nf = frd(N,omega);

% Matrix N
Delta1 = ultidyn('Delta1',[1 1]);
Delta2 = ultidyn('Delta2',[1 1]);
Delta = blkdiag(Delta1,Delta2);
Gp = C*(s*eye(5)-A_i)^(-1)*B_i;
G_pinv = inv(sys'*sys)*sys';
bodemag(usample(G_pinv*(Gp-Gp.NominalValue),100))
sigma(G_pinv*(Gp-sys));
hold on;
sigma(Wi);

%bodemag(sys_inv*tf(Gp-sys),'r'); hold on; bodemag(Wi,'b');
%bodemag(usample(sys_inv*(Gp-sys),50),'r'); hold on; bodemag(Wi,'b');
M = lft(Delta,N);
Mf = frd(M,omega);

% RS with mussv, rea
% M = N(1,1), per la robusta stabilità la norma infinito di M deve essere
% minore di 1
% N = lft(P,k) di dimensione 6x8
Nrs = Nf(1:4,1:2);
[mubnds,muinfo] = mussv(Nrs,[-1 0],'a');
muRS = mubnds(:,1);
[muRSinf,muRSw] = norm(muRS,inf);
% per la Performance Nominale devo fare un controllo sulla N22
Nnp=Nf(5:8,3:6); % Picking out wP*Si
[mubnds,muinfo]=mussv(Nnp,[2 2],'c');
muNP = mubnds(:,1);
[muNPinf,muNSw]=norm(muNP,inf);
% Mrp = Mf(5:6,5:6);

%% plots


%bodemag(mubnds(:,1),'r-'); hold on; bodemag(1/Wi(1,1),'g'); hold on; bodemag(frd(autovalori(1,:),omega),'b')
% figure(1);
% bodemag(mubnds(:,2),'r-'); hold on; bodemag(1/Wi(2,2),'g'); hold on; bodemag(frd(autovalori(1,:),omega),'b')
% figure(3);
%bodemag(mubnds,'r-'); hold on; bodemag(1/Wi,'g');
%legend('mu(T)', '1/|wp|', '$$\bar{\sigma}(T)$$', 'Interpreter','latex')

