close all

s = tf('s');
Jo = get_linearization();    % matrices of the linearized system
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
%AP = 3;
AP = 10^-2;
%wBp = 1;
wBp = 10^-2;
wBt = 1;
wBu = 1;
% Matrici di peso
Wu = s/(s+wBu);
wP = (s/(M)^1/2+wBp)^2/(s+wBp*(AP)^1/2)^2;
%wP = (s/M+wBp)/(s+wBp*AP);
%wT = s/(s+wBt);
WP = blkdiag(wP,wP,wP,wP);
WT = blkdiag(wT,wT,wT,wT);
WU = blkdiag(Wu,Wu);

% Controllore H inf
[khinf,ghinf,gopt] = mixsyn(sys,WP,WU,WT);
[Ahinf Bhinf Chinf Dhinf] = ssdata(khinf);


% Valori singolari
K = khinf;
S = inv(eye(4)+sys*K);
val_sing_S = sigma(S);
val_sing_max_S = val_sing_S(1,:);
omega_S = logspace(-1,6,146);
T = sys*K*S; %GK(I+GK)^-1
val_sing_T = sigma(T);
val_sing_max_T = val_sing_T(1,:);
omega_T = logspace(-1,6,251);
% sigma(S)
log_vars.K = K;
log_vars.sys = sys;
log_vars.S = S;
log_vars.T = T;
log_vars.WP = WP;
log_vars.WU = WU;

%Plot
figure(1);
sigma(S); hold on; sigma(1/WP);
% figure(2);
% sigma(T); hold on; sigma(1/WT);
% figure(3);
% sigma(K*S); hold on; sigma(1/WU);
% figure(4);
% sigma(S,'b',K*S,'r',T,'g',gopt/WP,'b-.',ss(gopt/WU),'m-.',gopt/WT,'g-.',{1e-3,1e3})
% legend('S','KS','T','GAM/W1','GAM/W2','GAM/W3','Location','SouthWest');
save('dataset','log_vars');













