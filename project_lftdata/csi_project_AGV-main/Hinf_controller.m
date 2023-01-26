close all
global rp w_rp
s = tf('s');
Jo = get_linearization();    %richiama la funzione per la linearizzazione del sistema
A = Jo.A;
B = Jo.B;
C = Jo.C;
D = Jo.D;
A_i = Jo.A_i; %Matrice della dinamica A incerta
B_i = Jo.B_i; %Matrice B incerta

%Costruzione della funzione di trasferimento incerta
%Gp = C*(s*eye(5)-A_i)^(-1)*B_i;
Gp = ss(A_i,B_i,C,D);
G_i= minreal(tf(Gp));
[Ap Bp Cp Dp] = ssdata(G_i);
Gi = minreal(ss(Ap,Bp,Cp,Dp));
%Definizione dei pesi di performance

%Costruzione della Wp
M = 1; %picco massimo di S che da prassi garantisce buoni margini di guadagno sul sistema
AP = 10^-1; %errore massimo a regime
wBp = 1; %frequenza minima di banda per la performance
wP = 5*10^-2*(s/M+wBp)/(s+wBp*AP); %peso sulla performance
%wP = (s/(M)^1/2+wBp)^2/(s+wBp*(AP)^1/2)^2; %wp per maggiore pendenza 
WP = blkdiag(wP,wP,wP,wP);

%Costruzione della WU
wBu = 1;
WU = 1/10*tf(eye(2)); %peso sullo sforzo di controllo
%WU = blkdiag(Wu,Wu);

%Costruzione della WT
wBt = 1; %frequenza minima di banda per attenuazione di rumore di misura
Wt = makeweight(0.0001,50,500);
%Wt = s/(s+wBt);%peso sul rumore di misura
WT = blkdiag(Wt,Wt,Wt,Wt);


%Funzione di trasferimento del sistema nominale
SYS = ss(A,B,C,D);   
Gnom = minreal(tf(SYS));
[Anom Bnom Cnom Dnom] = ssdata(Gnom);
sys = minreal(ss(Anom,Bnom,Cnom,Dnom));



%Definizione dei parametri
% M = 2;
% AP = 10^-2;
% wBp = 10^-2;
% wBt = 1;
% wBu = 1;
% % Matrici di peso
% Wu = s/(s+wBu);
% %wP = (s/M+wBp)/(s+wBp*AP);

%% Controllore H inf con MIXSYN

[K,ghinf,gopt] = mixsyn(Gp,WP,WU,WT);

%Simulink
[K_mix,ghinf_mix,gopt_mix] = mixsyn(sys,WP,WU,WT);
[Amix Bmix Cmix Dmix] = ssdata(K_mix);


%% Interconnessione e Controllore con HINF 

systemnames = 'Gp WP WU WT';
inputvar = '[w{4}; u{2}]';
outputvar = '[WP ; WU; WT; -w-Gp]';
input_to_Gp = '[u]';
input_to_WP = '[-w-Gp]';
input_to_WU = '[u]';
input_to_WT = '[Gp]';
sysoutname = 'P_i';
cleanupsysic = 'yes';
sysic;
[P, Delta, blk] = lftdata(P_i);

nmeas=4; %Numero uscite dall'impianto
nu=2; %Numero ingressi all'impianto
[Khinf,CL,gamma] = hinfsyn(P_i,nmeas,nu);
[Khinf,CL,gamma] = hinfsyn(P_i,nmeas,nu);
[Ahinf,Bhinf,Chinf,Dhinf] = ssdata(Khinf);
looptrans = loopsens(Gp,Khinf);

figure(1);
sigma(looptrans.So, 'b'); hold on; sigma(1/WP,'r-.');
legend('S','gopt/W1');
figure(2);
sigma(looptrans.To,'r'); hold on; sigma(1/WT,'r-.');
legend('T','gopt/WT','r-.');
figure(3);
sigma(Khinf*looptrans.So,'g'); hold on; sigma(1/WU,'g-.');
legend('KS','gopt/WU');
% sigma(S,'b',K*S,'r',T,'g',gopt/WP,'b-.',gopt/WU,'m-.',gopt/WT,'g-.',{1e-3,1e3})
% legend('S','KS','T','GAM/W1','GAM/W2','GAM/W3','Location','SouthWest');
save('dataset','log_vars');

%% MU-ANALISI
%Matrici dinamiche incerte
omega = logspace(-1,6,100);

% Struttura MDelta
N = lft(P,K_mix); %funzione di trasferimento tra w e z
Nf = frd(N,omega);%risposta in frequenza di N

M = lft(Delta,N); %Si ottiene la fdt tra udel e ydel
Mf = ufrd(M,omega); %risposta in frequenza

%% RS con mussv
%osservo autovalori di N per capire se è NS
eig(N);
% M = N(1,1), per la robusta stabilità la norma infinito di M deve essere
% minore di 1
Nrs = Nf(1:9,1:9);
[mubnds,muinfo] = mussv(Nrs,blk,'a');
muRS = mubnds(:,1);
[muRSinf,muRSw] = norm(muRS,inf);
semilogx(mubnds(:,1), mubnds(:,2));

%per la Performance Nominale devo fare un controllo sulla N22
% DELTAP è una matrice complessa piena delle stesse dimensioni
%di F' dove F è lft(N,Delta)=M;
%Quindi essendo M una matrice 6x4 M' è 4x6
Nnp=Nf(10:end,10:end); % Picking out wP*Si
[mubnds,muinfo]=mussv(Nnp,[4 10],'a');
muNP = mubnds(:,1);
[muNPinf,muNSw]=norm(muNP,inf);

%Performance robusta

[mubnds,muinfo]=mussv(Nf,[9 0;4 10],'a');
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
%[stabmarg, destabunc, report] = robuststab(Tif,opt) %incertezze massime che mi porterebbero all'instabilità

%for RS with new "robstab" passing the whole  Nf
[stabmarg, destabunc, info] = robstab(Mf, opt)

%RP analysis with robgain
Delta_P = ultidyn('Delta_P', [4 10]);
delta = blkdiag(Delta, Delta_P);
N = lft(P, K_mix);
F = lft(delta, N);
Ff = ufrd(F, omega);
[perfmarg, destabunc, info] = robgain(Ff, 1, opt)













