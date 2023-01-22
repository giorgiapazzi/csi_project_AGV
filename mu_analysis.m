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
Gp = C*(s*eye(5)-A_i)^(-1)*B_i; %G incerta
G_inv = inv(sys'*sys)*sys';
%sigma(G_inv*(Gp-sys));

%wi creata per punti con la ginput
%muRSinf = 0.0014; muNPinf = 0.1070, muRPinf = 38.8978
%non funziona con la dk
freq = [1.76236720195385e-14
3.70891120355911e-12
7.73255735566483e-11
9.88344501101653e-09
1.96923412162587e-07
3.12786351032379e-06
1.67375984003505e-05
0.000107371579399216
0.00178454496337455
0.0389306738745860
0.564763912507189
2.88818212427404
187.022601461826
8426.78501816061
784211.360930332
24586644.0157908
1159190427.60994
129323321380.526
6379991425502.34
452340045456298];

response = [899.186991869919
814.634146341463
758.265582655826
682.384823848238
636.856368563686
565.311653116531
489.430894308943
409.214092140921
274.796747967480
166.395663956640
49.3224932249320
14.6341463414633
-4.87804878048792
8.13008130081289
5.96205962059594
10.2981029810296
-0.542005420054466
3.79403794037921
31.9783197831976
36.3143631436312];
system = frd(10.^(response/20),freq);
%sigma(G_inv*(Gp-sys)); hold on; sigma(system);
% wi = fitmagfrd(system,2,2);
% wi = tf(wi);
% sigma(G_inv*(Gp-sys)); hold on; sigma(wi);  % forse dobbiamo stampare la wi per vedere se va bene invece che system

%%Funzione peso che sta sopra i valori singolari da 10^-5 in poi 
%muRSinf = 0.1322; muNPinf = 0.0729; muRPinf = 0.2841
%%Non funziona con la dk
% wi = 1/(1+s*10^10)^2*1/(1+s*10^6)^2*(1+s*10^3)^4*(1+s*10^17)^10*1/(1+s*10^15)^10;
% sigma(G_inv*(Gp-sys)); hold on; sigma(wi);



% %con questa funzione peso sto sopra ai valori singolari G_inv*(Gp-sys) da
% 10^-8 in poi, ottengo muRSinf=0.0138, muNPinf=0.0729, muRPinf = 0.1878;
% %Non funziona con la dk
% wi = 1/(1+s*10^8)^3*1/(1+s*10^6)^2*(1+s*10^3)^5*(1+s*10^17)^10*1/(1+s*10^15)^10;
% figure(1);
% sigma(G_inv*(Gp-sys)); hold on; sigma(wi);


% 
% wi il cui valore singolare sta sotto ai valori singolari di G_inv*(Gp-sys)
% muRSinf = 0.0028; muNPinf = 0.0729; muRPinf = 0.0782;
% Funziona con la dk
rp_tau = w_rp/(rp);
wi = rp_tau*rp*s/(1+rp*s);
sigma(G_inv*(Gp-sys)); hold on; sigma(wi);


%wi creata con ucover, funziona con la dk
% Garray = usample(Gp,50);
% orderWt = 2;
% Garrayg = frd(Garray,logspace(-3,3,60));
% [Usys,Info] = ucover(Garrayg,Gp.NominalValue,orderWt,'in');
% wi = tf(Info.W1); % Funziona con la dk
% sigma(G_inv*(Gp.NominalValue-Garray),'b--'); hold on; sigma(wi)

Wi = blkdiag(wi,wi);
log_vars.Wi = Wi;
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
outputvar = '[Wi ; WP ; WU; -w-sys]';
input_to_sys = '[u+udel]';
input_to_WP = '[sys]';
input_to_WU = '[u]';
input_to_Wi = '[u]';
sysoutname = 'P';
cleanupsysic = 'yes';
sysic;
%P = minreal(ss(P));
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
figure(2);
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
%eig(N)
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
save('dataset','log_vars');
%% plots

% figure(1);
% sigma(mubnds,'r-'); hold on; sigma(1/Wi,'g'); 
% bodemag(mubnds(:,2),'r-'); hold on; bodemag(1/Wi(2,2),'g'); hold on; bodemag(frd(autovalori(1,:),omega),'b')
% figure(3);
%bodemag(mubnds,'r-'); hold on; bodemag(1/Wi,'g');
%legend('mu(T)', '1/|wp|', '$$\bar{\sigma}(T)$$', 'Interpreter','latex')

