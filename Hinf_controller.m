close all

s = tf('s');
Jo = get_linearization_lqg();    % matrices of the linearized system
A = Jo.A;
B = Jo.B;
C = Jo.C;
D = Jo.D;
SYS = ss(A,B,C,D);    % transfer function
Gnom = minreal(tf(SYS));
[Anom Bnom Cnom Dnom] = ssdata(Gnom);
sys = minreal(ss(Anom,Bnom,Cnom,Dnom));

%Definizione dei parametri
M = 2;
AP = 2;
wBp = 0.01;
wBt = 1;
% Matrici di peso
Wu = s/(s+wBt);
wP = (s/(M)^1/2+wBp)^2/(s+wBp*(AP)^1/2)^2;
wT = s/(s+wBt);
% wP = makeweight(60,9e-2,0.1);
% wP2 = makeweight(250,7.5e-6,0.1);
% wP3 = makeweight(60,3e-5,0.1);
% wP4 = makeweight(250,3e-6,0.1);
% wT1 = makeweight(0.1,1.2,30);
% wT3 = makeweight(0.1,30,30);
% wT = blkdiag(wT1,wT3);
WP = blkdiag(wP,wP,wP,wP);
WT = blkdiag(wT,wT,wT,wT);
WU = blkdiag(Wu,Wu);
% bodemag(WP1); hold on; bodemag(wP3);
% wP2 = (s/M1+wB2)/(s+wB2*A1);
% wP3 = (s/M1+wB3)/(s+wB3*A1);
% wP4 = (s/M1+wB4)/(s+wB4*A1);
% WP =  blkdiag(wP1, wP2, wP3, wP4);

% Controllore H inf
[khinf,ghinf,gopt] = mixsyn(sys,WP,WU,WT);


% Valori singolari
K = khinf;
S = inv(eye(4)+sys*K);
T = sys*K*S;
% sigma(S)
log_vars.K = K;

%Plot
figure(1);
sigma(S); hold on; sigma(1/WP);
figure(2);
sigma(T); hold on; sigma(1/WT);
figure(3);
sigma(S,'b',K*S,'r',T,'g',gopt/WP,'b-.',ss(gopt/WU),'m-.',gopt/WT,'g-.',{1e-3,1e3})
legend('S','KS','T','GAM/W1','GAM/W2','GAM/W3','Location','SouthWest');
save('dataset','log_vars');
















