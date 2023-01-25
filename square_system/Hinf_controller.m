close all
global rp w_rp
s = tf('s');
Jo = get_linearization();    %richiama la funzione per la linearizzazione del sistema
A = Jo.A;
B = Jo.B;
C = Jo.C;
D = Jo.D;

SYS = ss(A,B,C,D);    % funzione di trasferiemnto del sistema nominale
Gnom = minreal(tf(SYS));
[Anom Bnom Cnom Dnom] = ssdata(Gnom);
sys = minreal(ss(Anom,Bnom,Cnom,Dnom));

% Dfinizione dei parametri di performance
M = 2; %picco massimo di S che da prassi garantisce buoni margini di guadagno sul sistema
AP = 1.2; %errore massimo a regime
wBp = 0.2; %frequenza minima di banda per la performance
%wP = (s/(M)^1/2+wBp)^2/(s+wBp*(AP)^1/2)^2; %wp per maggiore pendenza 
wP = (s/M+wBp)/(s+wBp*AP); %peso sulla performance
wBt = 1; %frequenza minima di banda per attenuazione di rumore di misura
% Matrici di peso
Wu = tf(1);  %peso sullo sforzo di controllo

%Definizione dei parametri
% M = 2;
% AP = 10^-2;
% wBp = 10^-2;
% wBt = 1;
% wBu = 1;
% % Matrici di peso
% Wu = s/(s+wBu);
% %wP = (s/M+wBp)/(s+wBp*AP);

wT = s/(s+wBt);%peso sul rumore di misura
WP = blkdiag(wP,wP);
WT = blkdiag(wT,wT);
WU = blkdiag(Wu,Wu);

%% Controllore H inf con MIXSYN
[khinf,ghinf,gopt] = mixsyn(sys,WP,WU,WT);
[Ahinf Bhinf Chinf Dhinf] = ssdata(khinf);


% Valori singolari
K = khinf;
S = inv(eye(2)+sys*K); %funzione di sensitività
omega = logspace(-1,6,146);
T = sys*K*S; %GK(I+GK)^-1 %funzione di sensitività complementare

%Plot
figure(1);
sigma(S, 'b'); hold on; sigma(1/WP,'b-.');
legend('S','gopt/W1');
figure(2);
sigma(T,'r'); hold on; sigma(1/WT,'r-.');
legend('T','gopt/WT','r-.');
figure(3);
sigma(K*S,'g'); hold on; sigma(1/WU,'g-.');
legend('KS','gopt/WU');
figure(4)
sigma(S,'b',K*S,'r',T,'g',gopt/WP,'b-.',gopt/WU,'m-.',gopt/WT,'g-.',{1e-3,1e3})
legend('S','KS','T','GAM/W1','GAM/W2','GAM/W3','Location','SouthWest');
save('dataset','log_vars');

%% Controllore con HINF (funziona peggio di mixsyn)

% nmeas=4;
% nu=2;
% [K,CL,gamma] = hinfsyn(P,nmeas,nu);


%% MU-ANALISI
%Matrici dinamiche incerte
A_i = Jo.A_i;
B_i = Jo.B_i;
omega = logspace(-1,6,302);
% Funzione di trasferimento incerta
Gp = C*(s*eye(5)-A_i)^(-1)*B_i; %G incerta
% Gp = usample(Gp,50);
% Gp = tf(Gp);
%sys = tf(sys);
G_inv = inv(sys);%pseudoinversa sinistra
%sigma((Gp-sys)*G_inv);

%wi creata per punti con la ginput
%muRSinf = 0.0014; muNPinf = 0.1070, muRPinf = 38.8978
%non funziona con la dk
% freq = [1.76236720195385e-14
% 3.70891120355911e-12
% 7.73255735566483e-11
% 9.88344501101653e-09
% 1.96923412162587e-07
% 3.12786351032379e-06
% 1.67375984003505e-05
% 0.000107371579399216
% 0.00178454496337455
% 0.0389306738745860
% 0.564763912507189
% 2.88818212427404
% 187.022601461826
% 8426.78501816061
% 784211.360930332
% 24586644.0157908
% 1159190427.60994
% 129323321380.526
% 6379991425502.34
% 452340045456298];
% 
% response = [899.186991869919
% 814.634146341463
% 758.265582655826
% 682.384823848238
% 636.856368563686
% 565.311653116531
% 489.430894308943
% 409.214092140921
% 274.796747967480
% 166.395663956640
% 49.3224932249320
% 14.6341463414633
% -4.87804878048792
% 8.13008130081289
% 5.96205962059594
% 10.2981029810296
% -0.542005420054466
% 3.79403794037921
% 31.9783197831976
% 36.3143631436312];
% system = frd(10.^(response/20),freq);
%sigma(G_inv*(Gp-sys)); hold on; sigma(system);
% wi = fitmagfrd(system,2,2);
% wi = tf(wi);
% sigma(G_inv*(Gp-sys)); hold on; sigma(wi);  % forse dobbiamo stampare la wi per vedere se va bene invece che system

%%Funzione peso che sta sopra i valori singolari da 10^-5 in poi 
%muRSinf = 0.1322; muNPinf = 0.0729; muRPinf = 0.2841
%%Non funziona con la dk
% wi = 1/(1+s*10^10)^2*1/(1+s*10^6)^2*(1+s*10^3)^4*(1+s*10^17)^10*1/(1+s*10^15)^10;
% sigma(G_inv*(Gp-sys)); hold on; sigma(wi);

%Funzione che sta tutta sopra ma che non dà valori buoni sulla mu sintesi

% %con questa funzione peso sto sopra ai valori singolari G_inv*(Gp-sys) da
% 10^-8 in poi, ottengo muRSinf=0.0138, muNPinf=0.0729, muRPinf = 0.1878;
% %Non funziona con la dk
% wi = 1/(1+s*10^8)^3*1/(1+s*10^6)^2*(1+s*10^3)^5*(1+s*10^17)^10*1/(1+s*10^15)^10;
% figure(1);
% sigma(G_inv*(Gp-sys)); hold on; sigma(wi);



%wi il cui valore singolare sta sotto ai valori singolari di G_inv*(Gp-sys)
%muRSinf = 0.0028; muNPinf = 0.0729; muRPinf = 0.0782;
%Funziona con la dk
% rp_tau = w_rp/(rp);
% wi = rp_tau*rp*s/(1+rp*s);
% sigma(G_inv*(Gp-sys)); hold on; sigma(wi);


% M = 2; %picco massimo di S che da prassi garantisce buoni margini di guadagno sul sistema
% AP = 0.01; %errore massimo a regime
% wBp = 1; %frequenza minima di banda per la performance
% wP = (s/(M)^1/2+wBp)^2/(s+wBp*(AP)^1/2)^2; %wp per maggiore pendenza
% sigma(G_inv*(Gp-sys)); hold on; sigma(wP);

% wi = makeweight(2.2,0.9,0.4);
% sigma(G_inv*(Gp-sys)); hold on; sigma(wi);

%wi creata con ucover, funziona con la dk
Garray = usample(Gp,50);
orderWt = 2;
Garrayg = frd(Garray,logspace(-3,3,60));
[Usys,Info] = ucover(Garrayg,Gp.NominalValue,orderWt,'in');
wi = tf(Info.W1); % Funziona con la dk
sigma(G_inv*(Gp.NominalValue-Garray),'b--'); hold on; sigma(wi,'r')

Wi = blkdiag(wi,wi);
%Realizzazione del sistema in forma di DOYLE
systemnames = 'sys WP WU Wi';
inputvar = '[udel{2}; w{2}; u{2}]';
outputvar = '[Wi ; WP ; WU; -w-sys]';
input_to_sys = '[u+udel]';
input_to_WP = '[sys+w]';
input_to_WU = '[u]';
input_to_Wi = '[u]';
sysoutname = 'P';
cleanupsysic = 'yes';
sysic;
P = minreal(ss(P));

% Struttura MDelta
N = lft(P,K); %funzione di trasferimento tra w e z
Nf = frd(N,omega);%risposta in frequenza di N

%Si definiscono le incertezze
Delta1 = ultidyn('Delta1',[1 1]);
Delta2 = ultidyn('Delta2',[1 1]);
%Incertezza strutturata 
Delta = blkdiag(Delta1,Delta2);

Gpp = sys*(eye(2)+Wi*Delta);%G incerta con incertezza moltiplicativa
figure(2);
sigma(Gpp, 'r'); hold on; sigma(sys, 'b');

M = lft(Delta,N);%Si ottiene la fdt tra udel e ydel
Mf = ufrd(M,omega);%risposta in frequenza

%% RS con mussv
%osservo autovalori di N per capire se è NS
eig(N);
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
Nnp=Nf(3:6,3:4); % Picking out wP*Si
[mubnds,muinfo]=mussv(Nnp,[2 4],'a');
muNP = mubnds(:,1);
[muNPinf,muNSw]=norm(muNP,inf);

%Performance robusta
[mubnds,muinfo]=mussv(Nf,[1 1;1 1;2 4],'a');
muRP = mubnds(:,1);
[muRPinf,muNSw]=norm(muRP,inf);

%% RS con robuststab
%altrimenti is può far applicare la robusta stabilità direttamente da
%matlab, 
% looptranfer = loopsens(Gp, K);
% Ti = looptranfer.Ti;
% Tif = ufrd(Ti, omega);
opt = robopt('Display','on');
% con il comando successivo mi dice sulla robusta stabilità
%in particolare stabmarg mi indica gli upper e lower bound, l'inverso del
%lower bound deve essere uguale al muRSinf ottenuto con il comando mussv
%destabunc mi indica esattamente l'incertezza massima tollerabile dal
%sistema oltre la quale non è robustamente stabile
[stabmarg, destabunc, report] = robuststab(Tif,opt) %incertezze massime che mi porterebbero all'instabilità

%for RS with new "robstab" passing the whole  Nf
[stabmarg, destabunc, info] = robstab(Mf, opt)

%RP analysis
[stabmarg, destabunc, info] = robustperf(Mf, opt)

%% plots

% figure(1);
% sigma(mubnds,'r-'); hold on; sigma(1/Wi,'g'); 
% bodemag(mubnds(:,2),'r-'); hold on; bodemag(1/Wi(2,2),'g'); hold on; bodemag(frd(autovalori(1,:),omega),'b')
% figure(3);
%bodemag(mubnds,'r-'); hold on; bodemag(1/Wi,'g');
%legend('mu(T)', '1/|wp|', '$$\bar{\sigma}(T)$$', 'Interpreter','latex')












