%Inizializzazione
init
s = tf('s');
Jo = get_linearization();    %richiama la funzione per la linearizzazione del sistema
A = Jo.A;
B = Jo.B;
C = Jo.C;
D = Jo.D;
A_i = Jo.A_i; %Matrice della dinamica A incerta
B_i = Jo.B_i; %Matrice B incerta

%Costruzione della funzione di trasferimento incerta
Gp = ss(A_i,B_i,C,D);
% G_i= minreal(tf(Gp));
% [Ap Bp Cp Dp] = ssdata(G_i);
% Gi = minreal(ss(Ap,Bp,Cp,Dp));

% Definizione dei pesi di performance:
% Costruzione della WP, peso sulla performance
M = 1;  % picco massimo di S che da prassi garantisce buoni margini di guadagno sul sistema
AP = 10^-1; % errore massimo a regime
wBp = 1;    % frequenza minima di banda per la performance
wP = 5*10^-2*(s/M+wBp)/(s+wBp*AP);
%wP = (s/(M)^1/2+wBp)^2/(s+wBp*(AP)^1/2)^2; % wp per maggiore pendenza 
WP = blkdiag(wP,wP,wP,wP);

%Costruzione della WU, peso sullo sforzo di controllo
wBu = 1;
WU = 1/20*tf(eye(2));

%Costruzione della WT, peso sul rumore di misura
wBt = 1;    % frequenza minima di banda per attenuazione di rumore di misura
Wt = makeweight(10^-2,20,500);
WT = blkdiag(Wt,Wt,Wt,Wt);


%Funzione di trasferimento del sistema nominale
SYS = ss(A,B,C,D);   
Gnom = minreal(tf(SYS));
[Anom Bnom Cnom Dnom] = ssdata(Gnom);
sys = minreal(ss(Anom,Bnom,Cnom,Dnom));


%Interconnessione
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
[P, Delta, blk] = lftdata(P_i); % estrae la matrice certa P e la matrice delle incertezze Delta

%% Controllore Hinf con MIXSYN

[Kmyx,ghinf,gopt] = mixsyn(Gp,WP,WU,WT);
looptrans = loopsens(Gp,Kmyx); % Funzione di trasferimento a ciclo chiuso

%Simulink
[K_mix,ghinf_mix,gopt_mix] = mixsyn(sys,WP,WU,WT);
[Amix Bmix Cmix Dmix] = ssdata(K_mix);

% Plot delle funzioni peso
figure(1);
sigma(looptrans.So, 'b'); hold on; sigma(1/WP,'r-.'); 
legend('S','gopt/W1'); hold off;
figure(2);
sigma(looptrans.To,'r'); hold on; sigma(1/WT,'r-.');
legend('T','gopt/WT','location','south'); hold off;
figure(3);
sigma(Kmyx*looptrans.So,'g'); hold on; sigma(1/WU,'g-.');
legend('KS','gopt/WU'); hold off;
% figure(4)
% sigma(S,'b',K*S,'r',T,'g',gopt/WP,'b-.',gopt/WU,'m-.',gopt/WT,'g-.',{1e-3,1e3})
% legend('S','KS','T','GAM/W1','GAM/W2','GAM/W3','Location','SouthWest');



%% MU-ANALISI
omega = logspace(-1,6,100);

% Modificare questa variabile per scegliere su quale controllore fare
% mu-analisi: Khinf per controllore con hinfsyn o Kmyx per controllore con mixsyn
K = Kmyx;

% Matrici dinamiche incerte
% Struttura MDelta
N = lft(P,K);   % funzione di trasferimento tra w e z
Nf = frd(N,omega);  % risposta in frequenza di N

M = lft(Delta,N);   % Si ottiene la fdt tra udel e ydel
Mf = ufrd(M,omega); % risposta in frequenza

%%%%%%%%%
% MUSSV
%%%%%%%%%

%osservo autovalori di N per capire se è NS
eig(N); % per la stabilità nominale devono essere tutti negativi

% M = N(1,1), per la robusta stabilità la norma infinito di M deve essere minore di 1
Nrs = Nf(1:9,1:9);
[mubnds,muinfo] = mussv(Nrs,blk,'a');
muRS = mubnds(:,1);
[muRSinf,muRSw] = norm(muRS,inf);

% per la performance nominale devo fare un controllo sulla N(2,2)
% DELTAP è una matrice complessa piena delle stesse dimensioni di F' dove F è lft(N,Delta) = M
% Quindi essendo M una matrice 10x4 M' è 4x10
Nnp=Nf(10:end,10:end);  % Picking out wP*Si
[mubnds,muinfo]=mussv(Nnp,[4 10],'a');
muNP = mubnds(:,1);
[muNPinf,muNSw]=norm(muNP,inf);

%Performance robusta: si fa il controllo su tutta la matrice N
[mubnds,muinfo]=mussv(Nf,[9 0;4 10],'a');
muRP = mubnds(:,1);
[muRPinf,muNSw]=norm(muRP,inf);

% Per controllare tali proprietà di robustezza stampiamo i valori muNpinf,
% muRSinf, muRPinf (devono essere tutti minori di 1 per avere buono risultati):
fprintf('\n*****************************************************************************************************\n')
fprintf('Gli autovalori del sistema sono:')
eig(feedback(sys*K,eye(4)))
fprintf('Stabilità robusta muRS: %f \n',muRSinf)
fprintf('Prestazione nominale muNPinf: %f \n',muNPinf)
fprintf('Robusta prestazione muRPinf: %f \n',muRPinf)
fprintf('*****************************************************************************************************\n\n')


%%%%%%%%%%%%%%%%%%%%%
% Robusta stabilità:
%%%%%%%%%%%%%%%%%%%%%
% looptranfer = loopsens(Gp, K);
% Ti = looptranfer.Ti;
% Tif = ufrd(Ti, omega);
opt = robopt('Display','on');
% con il comando successivo mi dice sulla robusta stabilità
%in particolare stabmarg mi indica gli upper e lower bound, l'inverso del
%lower bound deve essere uguale al muRSinf ottenuto con il comando mussv
%destabunc mi indica esattamente l'incertezza massima tollerabile dal
%sistema oltre la quale non è robustamente stabile

%%%%%%%%%%%%%%%%%%%%%%%
% Robusta Prestazione:
%%%%%%%%%%%%%%%%%%%%%%%
% Con robgain
Delta_P = ultidyn('Delta_P', [4 10]); % Incertezze fittizie
delta = blkdiag(Delta, Delta_P);    % Matrice delle incertezze finale
F = lft(delta, N);
Ff = ufrd(F, omega);

%for RS with new "robstab" passing the whole  Nf
[stabmarg, destabunc, info] = robstab(Mf, opt)
% Con robustperf
[stabmarg, destabunc, info] = robustperf(Ff, opt)