
%inizializzazione
init
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


% Funzione di trasferimento incerta
Gp = C*(s*eye(5)-A_i)^(-1)*B_i; %G incerta

%funzione di trasferiemnto del sistema nominale
SYS = ss(A,B,C,D);  
Gnom = minreal(tf(SYS));
[Anom Bnom Cnom Dnom] = ssdata(Gnom);
sys = minreal(ss(Anom,Bnom,Cnom,Dnom));
%Calcolo della pseudoinversa sinistra
G_inv = inv(sys'*sys)*sys';

%Peso sulla performance
M = 2; %picco massimo di S che da prassi garantisce buoni margini di guadagno sul sistema
AP = 3; %errore massimo a regime
wBp = 1; %frequenza minima di banda per la performance
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



%% Controllore H inf con MIXSYN
[Kmix,ghinf,gopt] = mixsyn(sys,WP,WU,WT);
[Amix, Bmix, Cmix, Dmix] = ssdata(Kmix);
%Funzione di trasferimento a ciclo chiuso
looptrans = loopsens(sys,Kmix);
%Plot
figure(1)
sigma(looptrans.So, 'b'); hold on; sigma(1/WP,'b-.'); 
legend('S','gopt/W1','location','south'); hold off;
figure(2)
sigma(looptrans.To,'r'); hold on; sigma(1/WT,'r-.'); 
legend('T','gopt/WT'); hold off;
figure(3)
sigma(Kmix*looptrans.So,'g'); hold on; sigma(1/WU,'g-.'); 
legend('KS','gopt/WU'); hold off; 
% figure(4)
% sigma(S,'b',K*S,'r',T,'g',gopt/WP,'b-.',gopt/WU,'m-.',gopt/WT,'g-.',{1e-3,1e3})
% legend('S','KS','T','GAM/W1','GAM/W2','GAM/W3','Location','SouthWest');



%% MU-ANALISI


%Funzione peso che sta sopra i valori singolari da 10^-5 in poi 
% wi = 1/(1+s*10^10)^2*1/(1+s*10^6)^2*(1+s*10^3)^4*(1+s*10^17)^10*1/(1+s*10^15)^10;
% sigma(G_inv*(Gp-sys)); hold on; sigma(wi);



%con questa funzione peso sto sopra ai valori singolari G_inv*(Gp-sys) da
%10^-8 in poi
% wi = 1/(1+s*10^8)^3*1/(1+s*10^6)^2*(1+s*10^3)^5*(1+s*10^17)^10*1/(1+s*10^15)^10;
% figure(1);
% sigma(G_inv*(Gp-sys)); hold on; sigma(wi);


%Funziona con la dk
% rp_tau = w_rp/(rp);
% wi = rp_tau*rp*s/(1+rp*s);
% sigma(G_inv*(Gp-sys)); hold on; sigma(wi);


%wi creata con ucover, funziona con la dk
Garray = usample(Gp,50);
orderWt = 2;
Garrayg = frd(Garray,logspace(-3,3,60));
[Usys,Info] = ucover(Garrayg,Gp.NominalValue,orderWt,'in');
wi = Info.W1; % Funziona con la dk
w_i = tf(wi);
figure(4)
sigma(G_inv*(Gp.NominalValue-Garray),'r--'); hold on; sigma(wi)

%Creazione di una Wi piena
Wi = [wi wi;wi wi];


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

N = lft(P,Kmix); %funzione di trasferimento tra [udel w] e [ydel z]
Nf = frd(N,omega);%risposta in frequenza di N

%Si definisce l'incertezza
Delta = ultidyn('Delta',[2 2]);

%Plot dei valori singolari della G incerta e della G nominale
Gpp = sys*(eye(2)+Wi*Delta);%G incerta con incertezza moltiplicativa
figure(5);
sigma(Gpp, 'r'); hold on; sigma(sys, 'b');

M = lft(Delta,N); %Si ottiene la fdt tra w e z, è la Fu(N,delta)
Mf = ufrd(M,omega); %risposta in frequenza

%% RS con mussv
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
fprintf('Gli autovalori del sistema sono:')
eig(feedback(sys*Kmix,eye(4)))
fprintf('Stabilità robusta muNP: %f \n',muRSinf)
fprintf('Prestazione nominale muNPinf: %f \n',muNPinf)
fprintf('Robusta prestazione muRPinf: %f \n',muRPinf)
fprintf('*****************************************************************************************************\n\n')
%% RS con robuststab

opt = robopt('Display','on');

%for RS with new "robstab" passing the whole  Nf
% [stabmarg, destabunc, info] = robstab(Mf, opt)

%RP analysis
[stabmarg, destabunc, info] = robustperf(Mf, opt)













