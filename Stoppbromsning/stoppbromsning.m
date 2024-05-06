%% stoppbromsning.m
% Huvudprogram f�r simulering av stoppbromsning med l�st skivbroms.
%
% Anders S�derberg, KTH - Maskinkonstruktion, 2018-08-31

%% St�dar
clear all; close all

%% Systemparametrar
% Definera systemparametar som globala
global r L1 L2 H1 H2 m g f c rho A mu

% Data f�r gokart
m   = 120;         % Massa gokart inkl. f�rare [kg]
r   = 0.2794/2;       % Hjulradie [m]
L1  = 0.658;       % Masscentrum relativt framaxel [m]
L2  = 0.600;       % Masscentrum relativit bakaxel [m]
H1  = 0.40;        % Masscentrum relativit markniv� [m]
H2  = 0.10;        % Angeppspunkt f�r Fluft relativit masscentrum [m]
A   = 0.5;         % Projicerad area f�r luftmotst�nd [m^2]
x2=0.5*10^-3;


% Data f�r bromsok/-cylinder
rb = 85.34 *10^-3;     % Mittradien f�r bromsbel�gget p� skivan [m]
dKok= 25.4 * 10^-3; % Diameter f�r kolv i bromsok [m] 
dKmcL=9/8*25.40*10^-3;    % Diameter f�r kolv i Mastercylinder[m]
dKmcS=5/8*25.4*10^-3;
Aok=(pi*dKok^2)/4;    % Kolvarea i ok [m^2]
AmcS=(pi*dKmcS^2)/4;    % Kolvarea i Mastercylinder [m^2]
AmcL=(pi*dKmcL^2)/4;
alfa= 16;               %Vinkel bromsok

% Fysikaliska parametrar
mu  = 0.8;         % Friktionstal mellan d�ck och v�gbana [-]
f   = 0.012;       % Rullmotst�ndskoefficient [-]
g   = 9.81;        % Tyngdacceleration [m/s^2]
c   = 0.6;         % Luftmotst�ndskoefficient [-]
rho = 1.22;        % Denstitet hos luft [kg/m^3]  
Mub=0.35;          % Friktionskoefficient mellan skiva och bel�gg [-]
%% Simulering av stoppbromsing
% Definiera tidsintervall och begynnelsevillkor
tstart = 0;                                                 % Starttid f�r simulering [s]
tend   = 10;                                                % Sluttid f�r simulering [s] 
vstart = 50/3.6;                                            % Utg�ngshastighet [m/s]
sstart = 0;                                                 % Utg�ngsstr�cka [m]
tspan  = [tstart tend];                                     % Start- och sluttid f�r simulering [s]
IC     = [vstart sstart];                                   % S�tter initialvillkor s� att vi bromsar fr�n utg�ngshastigheten och utg�ngsposition 
% L�s problem med l�mplig ode-l�sare
opt    = odeset('RelTol',1e-9);                             % S�tter toleranser p� l�saren
[T,Y] = ode45('derivatabromsning',tspan,IC,opt);            % Anropar ode-l�sare
% Dela upp tillst�ndsmatris p� vektorer
V      = Y(:,1);                                            % Hastighet [m/s]
S      = Y(:,2);                                            % Str�cka [m]
% Ber�kna kontaktkrafter och erforderlig bromsmoment
Frull  = f*m*g*(V>0);                                       % Rullmotst�nd [N]
Fluft  = 0.5*c*rho*A*V.^2.*(V>0);                           % Luftmotst�nd [N]
N2     = (m*g*L1-Frull*H1+Fluft*H2)./(L1+L2+mu*H1*(V>0));   % Normalkraft mellan bakd�ck och v�gbana [N]
N1     =  m*g-N2;                                           % Normalkraft mellan framd�ck och v�gbana [N]
F      = -mu*N2.*(V>0);                                     % Bromsande friktionskraft mellan bakd�ck och v�gbana [N]
M      =  mu*N2*r.*(V>0);                                   % Erforderligt bromsmoment f�r att l�sa bakaxeln [Nm]

% Utv�rdering av inbromsingsf�rlopp
i = find(V<1e-4,1);     % Hittar index till still�stende
t   = T(i);             % Inbromsningstid [s]
s   = S(i);             % Bromsstr�cka [m]
Mb  = max(M);           % Maximalt erfoderligt bromsmoment [Nm]

Ff= (Mb/rb)/2;          % Friktionskraften p� skivan fr�n ett bromsbel�gg [N]
Fk=Ff/Mub;              % Kl�mkraften fr�n oket [N]

%% Presentation av resultat
% Grafer som visar hastighet och str�cka
figure(1)
subplot(2,1,1)
plot(T,V*3.6,'LineWidth',2)
grid on
xlabel('t [s]')
ylabel('v [km/h]')
subplot(2,1,2)
plot(S,V*3.6,'LineWidth',2)
grid on
xlabel('s [m]')
ylabel('v [km/h]')
% Grafer som visar kontaktkrafter mellan d�ck och v�gbana samt erforderligt
% bromsmoment
figure(2)
subplot(2,1,1)
plot(T,[N1 N2],'LineWidth',2)
grid on
xlabel('t [s]')
ylabel('N_1, N_2 [N]')
legend('N_1','N_2')
subplot(2,1,2)
plot(T,M,'LineWidth',2)
grid on
xlabel('t [s]')
ylabel('M [Nm]')
% Utskrift i command window

psys=Fk/Aok;
FMcS=psys*AmcS;
FMcL=psys*AmcL;
xmcs=(2*Aok*x2)/AmcS;
xmcL=(2*Aok*x2)/AmcL;


disp(['Stopptid:     ' num2str(t) ' s'])
disp(['Stoppstr�cka: ' num2str(s) ' m'])
disp(['Bromsmoment:  ' num2str(Mb) ' Nm'])
disp(['Kl�mkraft:    ' num2str(Fk) ' N per sida'])
disp(['Bromstryck:   ' num2str(psys) 'Pa' ])

%% Kraftanalys bakaxel bromsning
bhl=108.1*10^-3;
blb=38.7*10^-3;
bld=202.38*10^-3;
b=879*10^-3;

Nb=max(N2); Ffr=-min(F);
H1z = Nb/2;
H2z = Nb/2;
H1y = -Ffr/2;
H2y= -Ffr/2;
H1x=0;
H2x=0;
Bz=(r/rb)*sind(alfa)*Ffr;
By=-(r/rb)*cosd(alfa)*Ffr;
Dy=0; 
Mh1y=r*(-H1x);
Mh2y=r*(-H2x);

L1x=-(H2x+H1x);
L2z=(Mh1y+bhl*H1z-blb*Bz-(b-bhl)*H2z+Mh2y)/(b-2*bhl);
L1z=-H1z-Bz-L2z-H2z;
L2y=(bhl*H1y-blb*By-(b-2*bhl-bld)*Dy-(b-bhl)*H2y)/(b-2*bhl);
L1y=-H1y-By-Dy-L2y-H2y;

L1R=sqrt(L1y^2+L1z^2); L2R=sqrt(L2y^2+L2z^2);
L1A=L1x;
Bromsvec=[L1A,L1R,L2R];
% 
% disp(['Lager kraft 1x: ' num2str(L1x) ' N'])
% disp(['Lager kraft 1y: ' num2str(L1y) ' N'])
% disp(['Lager kraft 1z: ' num2str(L1z) ' N'])
% disp(['Lager kraft 2y: ' num2str(L2y) ' N'])
% disp(['Lager kraft 2z: ' num2str(L2z) ' N'])
%% Kraftanlys bakaxel sv�ngning h�ger
bhl=108.1*10^-3;
blb=38.7*10^-3;
bld=202.38*10^-3;
b=879*10^-3;
Ky=0.11; Ki=0.89;

Nb=max(N2); Ffr=-min(F);
H1z = Nb*Ki;
H2z = Nb*Ky;
H1y = 0;
H2y= 0;
H1x=mu*H1z;
H2x=mu*H2z;
Bz=0;
By=0;
Dy=0; 
Mh1y=r*(-H1x);
Mh2y=r*(-H2x);

L1x=-(H2x+H1x);
L2z=(Mh1y+bhl*H1z-blb*Bz-(b-bhl)*H2z+Mh2y)/(b-2*bhl);
L1z=-H1z-Bz-L2z-H2z;
L2y=(bhl*H1y-blb*By-(b-2*bhl-bld)*Dy-(b-bhl)*H2y)/(b-2*bhl);
L1y=-H1y-By-Dy-L2y-H2y;

L1R=sqrt(L1y^2+L1z^2); L2R=sqrt(L2y^2+L2z^2);
L1A=L1x;
hsvangvec=[L1A,L1R,L2R];
%% Kraftanlys bakaxel sv�ngning v�nster
bhl=108.1*10^-3;
blb=38.7*10^-3;
bld=202.38*10^-3;
b=879*10^-3;

Nb=max(N2); Ffr=-min(F);
H1z = Nb*Ky;
H2z = Nb*Ki;
H1y = 0;
H2y= 0;
H1x=-mu*H1z;
H2x=-mu*H2z;
Bz=0;
By=0;
Dy=0; 
Mh1y=r*(-H1x);
Mh2y=r*(-H2x);

L1x=-(H2x+H1x);
L2z=(Mh1y+bhl*H1z-blb*Bz-(b-bhl)*H2z+Mh2y)/(b-2*bhl);
L1z=-H1z-Bz-L2z-H2z;
L2y=(bhl*H1y-blb*By-(b-2*bhl-bld)*Dy-(b-bhl)*H2y)/(b-2*bhl);
L1y=-H1y-By-Dy-L2y-H2y;

L1R=sqrt(L1y^2+L1z^2); L2R=sqrt(L2y^2+L2z^2);
L1A=L1x;
vsvangvec=[L1A,L1R,L2R];

%% Konstant hastightet


%% Positionering
% d=610*10^-3; % Avst�nd mellan lagern[m]
% Ff=Fk*Mub; % Friktionskraft
% Rx= @(alfa) Ff*cosd(alfa); Ry= @(alfa) Ff*sind(alfa);
% alfavec=linspace(0,2*pi,1000);
% dvec=linspace(0.1,d/2,1000);
% 
% 
% FLvvec=zeros(1,1000);FLhvec=zeros(1,1000);
% 
% for i=1:1000
% i=1; % Variation av rotationsvinkeln har ingen inverkan p� lagerkrafter
% alfa=alfavec(i);
% for ii=1:1000
% 
% 
% Sv=dvec(ii); % Avst�nd mellan bromsskiva och v�nstra lagret
% Sh=d-Sv; 
% FLvx=(Sh/d)*Rx(alfa); FLvy=(Sh/d)*Ry(alfa); 
% FLhx=(Sv/d)*Rx(alfa); FLhy=(Sv/d)*Ry(alfa);
% 
% FLvvec(ii)=sqrt(FLvx^2+FLvy^2);  FLhvec(ii)=sqrt(FLhx^2+FLhy^2);
% end
% 
% [FLv,ind1]=min(FLvvec); [FLh,ind2]=min(FLhvec); 
% 
% FLvec=FLhvec+FLvvec; [FL,ind]=min(FLvec);
% 
% disp(['Lagerkraft:   ' num2str(FLvvec(ind)) ' N p� v�nsterlager'])
% disp(['Lagerkraft:   ' num2str(FLhvec(ind)) ' N p� h�gerlager'])
% 
% 
% 
% 
%% Isaks test p� bakaxelber�kningar (M�ste ha k�rt "simulering av stoppbromsning" eftersom vi beh�ver bromskraften F.)
%Val av Y-lager:
%YAT 206 (sid 458 i SKF)
f_0 = 14; %Taget fr�n tabell 8 sid 445 i SKF
C_0 = 11.2 * 10^3; %[N]
C = 19.5 * 10^3; %[N]

p = 3; %F�r ber�kning av nominell livsl�ngd

%Parametrar som inte beror p� lastfall
bHL = 108.1*10^-3; %[m] Str�cka fr�n hjul till lager
bLB = 38.7*10^-3; %[m] Str�cka fr�n lager till bromsok
bLD = 202.38*10^-3; %[m] Str�cka fr�n lager till drev
b = 879*10^-3; %[m] Hela bakaxelns str�cka
rb = 91.43 *10^-3;%[m] Mittradien f�r bromsbel�gget p� skivan [m]
alfa = 16; %[grader] Vinkel f�r bromsoket med z-axeln

L10Vektor1 = [];
L10Vektor2 = [];
varvtalsVektor = [];




% Inbromsning parametrar
Dy = 0; %[N] Kraft fr�n draget
Pb = 0.42; %[1] procent av bilens vikt p� bakhjulen
Pv = 0.5; %[1] procent av bilens vikt p� de v�nstra hjulen
Ph = 0.5; %[1] procent av bilens vikt p� de h�gra hjulen
Fb = -min(F); %[N] Bromsande kraften
mux = 0; %[1] friktionstal i centrifugalriktningen

[H1z,H2z,H1y,H2y,H1x,H2x,By,Bz,MH1y,MH2y] = bakaxelkrafter(Pb,Pv,Ph,Fb,mux,rb,alfa);
[L1x,L2z,L1z,L2y,L1y] = lagerkrafter(H1z,H2z,H1y,H2y,H1x,H2x,By,Bz,MH1y,MH2y, Dy, bHL, bLB, bLD, b);
Fr1 = sqrt(L1z^2+L1y^2);
Fa1 = abs(L1x);
Fr2 = sqrt(L2z^2+L2y^2);
Fa2 = 0;

if (Fa1==0)
    P1 = Fr1;
else
    P1 = ekvivalentDynamiskLagerbelastning(Fr1,Fa1,f_0,C_0);
end

L101 = nominellLivslangd(C,P1,p);
P2 = Fr2;
L102 = nominellLivslangd(C,P2,p);

L10Vektor1 = [L10Vektor1, L101];
L10Vektor2 = [L10Vektor2, L102];
varvtalsVektor = [varvtalsVektor, 0.15];




% V�nstersv�ng parametrar
Dy = 0; %[N] Kraft fr�n draget
Pb = 0.52; %[1] procent av bilens vikt p� bakhjulen. Taget fr�n systemanalys Resultat-dokumentet
Pv = 0.11; %[1] procent av bilens vikt p� de v�nstra hjulen
Ph = 0.89; %[1] procent av bilens vikt p� de h�gra hjulen
Fb = 0; %[N] Bromsande kraften
mux = -mu; %[1] friktionstal i centrifugalriktningen

[H1z,H2z,H1y,H2y,H1x,H2x,By,Bz,MH1y,MH2y] = bakaxelkrafter(Pb,Pv,Ph,Fb,mux,rb,alfa);
[L1x,L2z,L1z,L2y,L1y] = lagerkrafter(H1z,H2z,H1y,H2y,H1x,H2x,By,Bz,MH1y,MH2y, Dy, bHL, bLB, bLD, b);
Fr1 = sqrt(L1z^2+L1y^2);
Fa1 = abs(L1x);
Fr2 = sqrt(L2z^2+L2y^2);
Fa2 = 0;

if (Fa1==0)
    P1 = Fr1;
else
    P1 = ekvivalentDynamiskLagerbelastning(Fr1,Fa1,f_0,C_0);
end
L101 = nominellLivslangd(C,P1,p);
P2 = Fr2;
L102 = nominellLivslangd(C,P2,p);

L10Vektor1 = [L10Vektor1, L101];
L10Vektor2 = [L10Vektor2, L102];
varvtalsVektor = [varvtalsVektor, 0.25];

% H�gersv�ng parametrar
Dy = 0; %[N] Kraft fr�n draget
Pb = 0.52; %[1] procent av bilens vikt p� bakhjulen. Taget fr�n systemanalys Resultat-dokumentet
Pv = 0.89; %[1] procent av bilens vikt p� de v�nstra hjulen
Ph = 0.11; %[1] procent av bilens vikt p� de h�gra hjulen
Fb = 0; %[N] Bromsande kraften
mux = mu; %[1] friktionstal i centrifugalriktningen

[H1z,H2z,H1y,H2y,H1x,H2x,By,Bz,MH1y,MH2y] = bakaxelkrafter(Pb,Pv,Ph,Fb,mux,rb,alfa);
[L1x,L2z,L1z,L2y,L1y] = lagerkrafter(H1z,H2z,H1y,H2y,H1x,H2x,By,Bz,MH1y,MH2y, Dy, bHL, bLB, bLD, b);
Fr1 = sqrt(L1z^2+L1y^2);
Fa1 = abs(L1x);
Fr2 = sqrt(L2z^2+L2y^2);
Fa2 = 0;

if (Fa1==0)
    P1 = Fr1;
else
    P1 = ekvivalentDynamiskLagerbelastning(Fr1,Fa1,f_0,C_0);
end
L101 = nominellLivslangd(C,P1,p);
P2 = Fr2;
L102 = nominellLivslangd(C,P2,p);

L10Vektor1 = [L10Vektor1, L101];
L10Vektor2 = [L10Vektor2, L102];
varvtalsVektor = [varvtalsVektor, 0.25];

% Acceleration parametrar



% Konstant hastighet parametrar?

%Ber�kning av livsl�ngd med delskadeteori?
L101tot = livslangdVarierandeLagerbelastning(L10Vektor1, varvtalsVektor)
L102tot = livslangdVarierandeLagerbelastning(L10Vektor2, varvtalsVektor)

mil = 2*r*pi*min(L101tot,L102tot)*10^6/10000
%% Funktioner
function [H1z,H2z,H1y,H2y,H1x,H2x,By,Bz,MH1y,MH2y] = bakaxelkrafter(Pb,Pv,Ph,Fb,mux,rb,alfa)
global  r m g
rh = r; %[m] Hjulens radie. r �r en global variabel.
H1z = m*g*Pb*Pv; %[N] Normalkraft(z-riktning) p� v�nster bakhjul. m och g �r globala variabler
H2z = m*g*Pb*Ph; %[N] Normalkraft(z-riktning) p� h�ger bakhjul. m och g �r globala variabler
H1y = -Fb/2; %[N] Kraft p� v�nster bakhjul i k�rriktning(y-riktning). Fb �r bromskraften.
H2y = H1y; %[N] Kraft p� v�nster bakhjul i k�rriktning(y-riktning). Fb �r bromskraften.
H1x = mux * H1z; %[N] Kraft p� v�nster bakhjul i x-led.
H2x = mux * H2z; %[N] Kraft p� h�ger bakhjul i -led.
By = rh/rb * -Fb * cosd(alfa); %[N] Kraften fr�n bromsoket i y-led.
Bz = rh/rb * Fb * sind(alfa); %[N] Kraften fr�n bromsoket i z-led.
MH1y = rh * -H1x; %[Nm] Momentet p� v�nster hjul i y-led.
MH2y = rh * -H2x; %[Nm] Momentet p� h�ger hjul i y-led.
end

function [L1x,L2z,L1z,L2y,L1y] = lagerkrafter(H1z,H2z,H1y,H2y,H1x,H2x,By,Bz,MH1y,MH2y, Dy, bHL, bLB, bLD, b)
global  r m g
L1x = -(H1x+H2x); %[N] Det styrande(v�nster) lagrets kraft i x-riktning.
L2z = (MH1y+bHL*H1z-bLB*Bz-(b-bHL)*H2z+MH2y)/(b-2*bHL);  %[N] H�ger lagrets kraft i z-led (upp)
L1z = -H1z-Bz-L2z-H2z; %[N] V�nster lagrets kraft i z-led (upp)
L2y = (bHL*H1y-bLB*By-(b-2*bHL-bLD)*Dy-(b-bHL)*H2y)/(b-2*bHL);
L1y = -H1y-By-Dy-L2y-H2y;
end

function y = linjarInterpolering(x,x_vektor,y_vektor)
    i = 1;
    while x > x_vektor(i) 
        i = i+1;
    end
    x_2 = x_vektor(i);
    x_1 = x_vektor(i-1);
    
    y_2 = y_vektor(i);
    y_1 = y_vektor(i-1);

    y = y_1 + (y_2-y_1)/(x_2-x_1) * (x-x_1);
end

function P = ekvivalentDynamiskLagerbelastning(F_r,F_a,f_0,C_0)

    t = f_0 * F_a/C_0
    t_vektor = [0.172, 0.345, 0.689, 1.03, 1.38, 2.07, 3.45, 5.17, 6.89];
    
    %F�r sp�rkullager:
    %e_vektor = [0.19, 0.22, 0.26, 0.28, 0.30, 0.34, 0.38, 0.42, 0.44];
    %X_vektor = [0.56, 0.56, 0.56, 0.56, 0.56, 0.56, 0.56, 0.56, 0.56]; %samma f�r alla t!
    %Y_vektor = [2.30, 1.99, 1.71, 1.55, 1.45, 1.31, 1.15, 1.04, 1.00];
    
    %F�r Y-lager f�rutom 17262 och 17263:
    e_vektor = [0.29, 0.32, 0.36, 0.38, 0.4, 0.44, 0.49, 0.54, 0.54];
    X_vektor = [0.46, 0.46, 0.46, 0.46, 0.46, 0.46, 0.46, 0.46, 0.46]; %samma f�r alla t!
    Y_vektor = [1.88, 1.71, 1.52, 1.41, 1.34, 1.23, 1.1, 1.01, 1.00];
    
    e = linjarInterpolering(t,t_vektor,e_vektor);
    X = linjarInterpolering(t,t_vektor,X_vektor);
    Y = linjarInterpolering(t,t_vektor,Y_vektor);
    
    if (F_a/F_r <= e)
        P = F_r;
    else
        P = X*F_r + Y*F_a;
    end
end

function L10 = nominellLivslangd(C,P,p)
    % p = 3 f�r kullager och p = 10/3 f�r rullager
    L10 = (C/P)^p
end

function L10 = livslangdVarierandeLagerbelastning(L10Vektor, varvtalsVektor)
    %Tagen fr�n sid 81 i SKF-katalogen.
    N = sum(varvtalsVektor);
    L10 = 1/(sum(varvtalsVektor/N ./L10Vektor));
end