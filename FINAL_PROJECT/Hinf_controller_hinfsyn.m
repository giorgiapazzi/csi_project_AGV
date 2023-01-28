%inizializzazione
init
close all

s = tf('s');
Jo = get_linearization();    %richiama la funzione per la linearizzazione del sistema
%Matrici della dinamica certe
A = Jo.A; 
B = Jo.B;
C = Jo.C;
D = Jo.D;

%Matrici dinamiche incerte
A_i = Jo.A_i;
B_i = Jo.B_i;
omega = logspace(-1,6,146);



%funzione di trasferiemnto del sistema nominale
SYS = ss(A,B,C,D);  
Gnom = minreal(tf(SYS));
[Anom Bnom Cnom Dnom] = ssdata(Gnom);
sys = minreal(ss(Anom,Bnom,Cnom,Dnom));

% Funzione di trasferimento incerta
Gp = C*(s*eye(5)-A_i)^(-1)*B_i; %G incerta
G_inv = inv(sys'*sys)*sys';%pseudoinversa sinistra

%Peso sulla performance
M = 2; %picco massimo di S che da prassi garantisce buoni margini di guadagno sul sistema
AP = 3; %errore massimo a regime
wBp = 0.1; %frequenza minima di banda per la performance
%wP = (s/(M)^1/2+wBp)^2/(s+wBp*(AP)^1/2)^2; %wp per maggiore pendenza 
wP = (s/M+wBp)/(s+wBp*AP); %peso sulla performance
WP = blkdiag(wP,wP,wP,wP);


% Peso sulla WU
Wu = tf(1);  %peso sullo sforzo di controllo
WU = blkdiag(Wu,Wu);

%Peso sulla WT
wBt = 1; %frequenza minima di banda per attenuazione di rumore di misura
wT = s/(s+wBt);%peso sul rumore di misura
WT = blkdiag(wT,wT,wT,wT);

%% Costruzione della Wi

%wi creata con ucover, funziona con la dk
Garray = usample(Gp,50);
orderWt = 2;
Garrayg = frd(Garray,logspace(-3,3,60));
[Usys,Info] = ucover(Garrayg,Gp.NominalValue,orderWt,'in');
wi = Info.W1; % Funziona con la dk
figure(1)
sigma(G_inv*(Gp.NominalValue-Garray),'r--'); hold on; sigma(wi)

%Creazione di una Wi piena
Wi = [wi wi;wi wi];

%% Forma di Doyle

%Realizzazione del sistema in forma di DOYLE
systemnames = 'sys Wi WP WU WT';
inputvar = '[udel{2}; w{4}; u{2}]';
outputvar = '[Wi ; WP ; WU; WT; -w-sys]';
input_to_sys = '[u+udel]';
input_to_WP = '[sys]';
input_to_WU = '[u]';
input_to_WT = '[sys]';
input_to_Wi = '[udel]';
sysoutname = 'P';
cleanupsysic = 'yes';
sysic;
P = minreal(ss(P));

%% Costruzione controllore con Hinfsyn
nmeas=4;
nu=2;
[K,CL,gamma] = hinfsyn(P,nmeas,nu);

%% Verifica delle funzioni peso
looptrans = loopsens(sys,K);
%Plot
figure(2);
sigma(looptrans.So, 'b'); hold on; sigma(1/WP,'b-.');
legend('S','gopt/W1','location','south'); hold off;
figure(3);
sigma(looptrans.To,'r'); hold on; sigma(1/WT,'r-.');
legend('T','gopt/WT'); hold off;
figure(4);
sigma(K*looptrans.So,'g'); hold on; sigma(1/WU,'g-.');
legend('KS','gopt/WU'); hold off;

%% Costruzione di N e M
N = lft(P,K); %funzione di trasferimento tra w e z
Nf = frd(N,omega);%risposta in frequenza di N

%Si definisce l'incertezza
Deltai = ultidyn('Deltai',[2 2]);
DeltaP = ultidyn('DeltaP',[4 10]);
Delta = blkdiag(Deltai,DeltaP);
%Plot dei valori singolari della G incerta e della G nominale
Gpp = sys*(eye(2)+Wi*Deltai);%G incerta con incertezza moltiplicativa
figure(5);
sigma(Gpp, 'r'); hold on; sigma(sys, 'b'); hold off;

M = lft(Delta,N); %Si ottiene la fdt tra udel e ydel
Mf = ufrd(M,omega); %risposta in frequenza

%% Mu analisi
%osservo autovalori di N per capire se è NS
eig(N);


% M = N(1,1), per la robusta stabilità la norma infinito di M deve essere
% minore di 1
Nrs = Nf(1:2,1:2);
[mubnds,muinfo] = mussv(Nrs,[2 2],'a');
muRS = mubnds(:,1);
[muRSinf,muRSw] = norm(muRS,inf);

%per la Performance Nominale devo fare un controllo sulla N22
% DELTAP è una matrice complessa piena delle stesse dimensioni
%di F' dove F è lft(N,Delta)=M;
Nnp=Nf(3:12,3:6); % Picking out wP*Si
[mubnds,muinfo]=mussv(Nnp,[4 10],'a');
muNP = mubnds(:,1);
[muNPinf,muNSw]=norm(muNP,inf);

%Performance robusta
[mubnds,muinfo]=mussv(Nf,[2 2;4 10],'a');
muRP = mubnds(:,1);
[muRPinf,muNSw]=norm(muRP,inf);

fprintf('\n*****************************************************************************************************\n')
fprintf('Stabilità robusta muRS: %f \n',muRSinf)
fprintf('Prestazione nominale muNPinf: %f \n',muNPinf)
fprintf('Robusta prestazione muRPinf: %f \n',muRPinf)
fprintf('*****************************************************************************************************\n\n')

%% RS con robuststab
opt = robopt('Display','on');
N = lft(P,K);
F = lft(Delta,N); %serve per performance robusta
Ff = ufrd(F,omega);
M = lft(Deltai,N); %serve per stabilità robusta
Mf = ufrd(M,omega);
[stabmarg, destabunc, report] = robstab(Mf,opt)% Robusta performance
[stabmarg, destabunc, info] = robustperf(Ff, opt)