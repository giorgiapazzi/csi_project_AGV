%% Nominal plant
close all
s = tf('s');
Jo = get_linearization_lqg();    % matrices of the linearized system
A = Jo.A;
B = Jo.B;
C = Jo.C;
D = Jo.D;
SYS = ss(A,B,C,D);    % transfer function
G = minreal(Gnom);


%% Controller
load('dataset');
K = log_vars.K;



[A B C D] = ssdata(K);
% Ks = C*(s*eye(6)-A)^(-1)*B;
% Ks = minreal(Ks);
[num1,den1] = ss2tf(A,B,C,D,1);
[num2,den2] = ss2tf(A,B,C,D,2);
K1 = tf(num1(1,:),den1);
K2 = tf(num1(2,:),den1);
K3 = tf(num2(1,:),den2);
K4 = tf(num2(2,:),den2);
Kt = [K1 K3;K2 K4];
Kt = minreal(Kt);

S = inv(eye(2)+G*K);
eig(S)
T = S*G*K;
autovalori = sigma(T);
sigma(T)
omega = logspace(-1,6,302);
%% Weighting filter for uncertainty modelling
Wi = log_vars.Wi

wP1 = makeweight(30,0.09,0.1);
wP2 = makeweight(30,0.75,0.5);
wP3 = makeweight(30,0.6,0.5);
WP =  [wP1 0;0 wP3];

Wu = eye(2);

%% Generalized plant P
systemnames = 'G WP Wu Wi';
inputvar = '[udel{2}; ref{2}; w{2}; u{2}]';
outputvar = '[Wi ; WP ; Wu ; ref-G-w]';
input_to_G = '[u+udel]';
input_to_WP = '[G+w]';
input_to_Wu = '[u]';
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
M = lft(Delta,N);
Mf = frd(M,omega);

% RS with mussv, rea
Mrs = Nf(1:2,1:2);
[mubnds,muinfo] = mussv(Mrs,[1 1 ; 1 1],'a');
muRS = mubnds(:,1);
[muRSinf,muRSw] = norm(muRS,inf);
Nnp=Nf(5:6,5:6); % Picking out wP*Si
[mubnds,muinfo]=mussv(Nnp,[1 1;1 1],'c');
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

