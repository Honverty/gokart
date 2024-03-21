%% stoppbromsning.m
% Huvudprogram för simulering av stoppbromsning med låst skivbroms.
%
% Anders Söderberg, KTH - Maskinkonstruktion, 2018-08-31

%% Städar
clear all; close all

%% Systemparametrar
% Definera systemparametar som globala
global r L1 L2 H1 H2 m g f c rho A mu

% Data för gokart
m   = 120;         % Massa gokart inkl. förare [kg]
r   = 0.135;       % Hjulradie [m]
L1  = 0.658;       % Masscentrum relativt framaxel [m]
L2  = 0.600;       % Masscentrum relativit bakaxel [m]
H1  = 0.40;        % Masscentrum relativit marknivå [m]
H2  = 0.10;        % Angeppspunkt för Fluft relativit masscentrum [m]
A   = 0.5;         % Projicerad area för luftmotstånd [m^2]
x2=0.5*10^-3;


% Data för bromsok/-cylinder
rb = 90 *10^-3;     % Mittradien för bromsbelägget på skivan [m]
dKok= 25.4 * 10^-3; % Diameter för kolv i bromsok [m] 
dKmcL=9/8*25.40*10^-3;    % Diameter för kolv i Mastercylinder[m]
dKmcS=5/8*25.4*10^-3;
Aok=(pi*dKok^2)/4;    % Kolvarea i ok [m^2]
AmcS=(pi*dKmcS^2)/4;    % Kolvarea i Mastercylinder [m^2]
AmcL=(pi*dKmcL^2)/4;

% Fysikaliska parametrar
mu  = 0.8;         % Friktionstal mellan däck och vägbana [-]
f   = 0.012;       % Rullmotståndskoefficient [-]
g   = 9.81;        % Tyngdacceleration [m/s^2]
c   = 0.6;         % Luftmotståndskoefficient [-]
rho = 1.22;        % Denstitet hos luft [kg/m^3]  
Mub=0.35;          % Friktionskoefficient mellan skiva och belägg [-]
%% Simulering av stoppbromsing
% Definiera tidsintervall och begynnelsevillkor
tstart = 0;                                                 % Starttid för simulering [s]
tend   = 10;                                                % Sluttid för simulering [s] 
vstart = 50/3.6;                                            % Utgångshastighet [m/s]
sstart = 0;                                                 % Utgångssträcka [m]
tspan  = [tstart tend];                                     % Start- och sluttid för simulering [s]
IC     = [vstart sstart];                                   % Sätter initialvillkor så att vi bromsar från utgångshastigheten och utgångsposition 
% Lös problem med lämplig ode-lösare
opt    = odeset('RelTol',1e-9);                             % Sätter toleranser på lösaren
[T,Y] = ode45('derivatabromsning',tspan,IC,opt);            % Anropar ode-lösare
% Dela upp tillståndsmatris på vektorer
V      = Y(:,1);                                            % Hastighet [m/s]
S      = Y(:,2);                                            % Sträcka [m]
% Beräkna kontaktkrafter och erforderlig bromsmoment
Frull  = f*m*g*(V>0);                                       % Rullmotstånd [N]
Fluft  = 0.5*c*rho*A*V.^2.*(V>0);                           % Luftmotstånd [N]
N2     = (m*g*L1-Frull*H1+Fluft*H2)./(L1+L2+mu*H1*(V>0));   % Normalkraft mellan bakdäck och vägbana [N]
N1     =  m*g-N2;                                           % Normalkraft mellan framdäck och vägbana [N]
F      = -mu*N2.*(V>0);                                     % Bromsande friktionskraft mellan bakdäck och vägbana [N]
M      =  mu*N2*r.*(V>0);                                   % Erforderligt bromsmoment för att låsa bakaxeln [Nm]

% Utvärdering av inbromsingsförlopp
i = find(V<1e-4,1);     % Hittar index till stillåstende
t   = T(i);             % Inbromsningstid [s]
s   = S(i);             % Bromssträcka [m]
Mb  = max(M);           % Maximalt erfoderligt bromsmoment [Nm]

Ff= (Mb/rb)/2;          % Friktionskraften på skivan från ett bromsbelägg [N]
Fk=Ff/Mub;              % Klämkraften från oket [N]

%% Presentation av resultat
% Grafer som visar hastighet och sträcka
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
% Grafer som visar kontaktkrafter mellan däck och vägbana samt erforderligt
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
disp(['Stoppsträcka: ' num2str(s) ' m'])
disp(['Bromsmoment:  ' num2str(Mb) ' Nm'])
disp(['Klämkraft:    ' num2str(Fk) ' N per sida'])


%% Positionering
d=610*10^-3; % Avstånd mellan lagern[m]
Ff=Fk*Mub; % Friktionskraft
Rx= @(alfa) Ff*cosd(alfa); Ry= @(alfa) Ff*sind(alfa);
alfavec=linspace(0,2*pi,1000);
dvec=linspace(0.1,d/2,1000);


FLvvec=zeros(1,1000);FLhvec=zeros(1,1000);

%for i=1:1000
i=1; % Variation av rotationsvinkeln har ingen inverkan på lagerkrafter
alfa=alfavec(i);
for ii=1:1000


Sv=dvec(ii); % Avstånd mellan bromsskiva och vänstra lagret
Sh=d-Sv; 
FLvx=(Sh/d)*Rx(alfa); FLvy=(Sh/d)*Ry(alfa); 
FLhx=(Sv/d)*Rx(alfa); FLhy=(Sv/d)*Ry(alfa);

FLvvec(ii)=sqrt(FLvx^2+FLvy^2);  FLhvec(ii)=sqrt(FLhx^2+FLhy^2);
end

[FLv,ind1]=min(FLvvec); [FLh,ind2]=min(FLhvec); 

FLvec=FLhvec+FLvvec; [FL,ind]=min(FLvec);

disp(['Lagerkraft:   ' num2str(FLvvec(ind)) ' N på vänsterlager'])
disp(['Lagerkraft:   ' num2str(FLhvec(ind)) ' N på högerlager'])







