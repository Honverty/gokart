%% LastfordelningAxlar
% Program för bestämning av lastfördelning mellan fram och bakaxel på
% gokarten under olika körfall
%
% Anders Söderberg, KTH-Maskinkonstruktion, 2018-08-31
%

%% Städar
clear all
close all
clc

%% Givna in data
% Massa och masscentrum för gokart inkl. förare
L1   = 0.658;
L2   = 0.600;
B1   = 0.410;               % Masscentrums läge relativit vänster hjul [m]
B2   = 0.410;               % Masscentrums läge relativit höger hjul [m]
H1   = 0.400;               % Masscentrums läge över marknivån [m]
m    = 120;                 % Massa inkl. förare [kg]
r    = 0.135;               % Hjulradie [mm]

% Data för drivlinan
u    = 2;                   % Utväxling [-]
eta  = 0.95;                % Verkningsgrad [-]
Mmax = 11.2;                % Maximalt motormoment enligt datablad [Nm]

% Fysikaliska parametrar
mu   = 0.8;                 % Friktionstal mellan däck och vägbana [ - ]
g    = 9.81;                % Tyngdacceleration [m/s^2]


%% Beräkning av lastfördelning mellan fram och bakaxel för olika körfall
% Acceleration med maximalt motormoment
F(1) = u*eta*Mmax/r;                    % Drivande kraft [N]
a(1) = F(1)/m;                          % Acceleration [m/s^2]
N_fram(1)= (m*g*L2 -F(1)*H1)/(L1+L2);   % Normalkrafter fram [N]
N_bak(1)= m*g -N_fram(1);               % Normalkrafter bak [N]

% Konstant hastighet
F(2) = 0;
a(2) = F(2)/m;
N_fram(2)= (m*g*L2 -F(2)*H1)/(L1+L2);   % Normalkrafter fram [N]
N_bak(2)= m*g -N_fram(2);               % Normalkrafter bak [N]

% Acceleration med maximalt motormoment
N_bak(3)= m*g*L1/(L1+L2+mu*H1) ;        % Normalkrafter bak [N]
N_fram(3)= m*g -N_bak(3);               % Normalkrafter fram [N]
F(3) = -mu*N_bak(3);
a(3) = F(3)/m;

% Tar fram lastfördelning mellan fram och bak
Last_fram = zeros(3,1); Last_bak = Last_fram;

for i = 1:3
    Last_fram(i) = N_fram(i)/(m*g);
    Last_bak(i) = N_bak(i)/(m*g);
end

% Tabell med lastfördelning
Fall_axel ={'Acceleration        ' 'Konstant hastighet  ' 'Inbromsning         '};
disp(['                    ' 'Fram   ' 'Bak'])
for i =1:3
    disp([Fall_axel{i} num2str(round(Last_fram(i)*100)) '%      ' num2str(round(Last_bak(i)*100)) '% '])
end
disp(' ')

% Tabell med normalkrafterna
disp(['                    ' 'Fram   ' 'Bak'])
for i =1:3
    disp([Fall_axel{i} num2str(round(N_fram(i))) ' N     ' num2str(round(N_bak(i))) ' N'])
end
disp(' ')



%%  LastfordelningSida.m
% Program för beräkning av lastfördelning och maximal hastighet vid kurvtagning
% 
% Anders Söderberg KTH Maskinkonstruktion 2018-08-31
disp('- - - - - - - - - - - - - - - - - - - - - - - - - - -')
%% Bestämning av masscentrums maximala höjd över marken för att undvika tippning

Bmin = min(B1,B2);              % Avgör minsta avstånd från masscentrum till hjul [m]
H1max = Bmin/mu;                % Beräknar maximalt tillåtet värde på H1

%% Bestämning av maximal kurvhastighet för att undvika sladd

R = [5:20];                     % Kurvradier [m]
v  = sqrt(mu*g*R);              % Maximal hastiget för att undvika sladd [m/s]

%% Bestämning av normalkrafter i olika körfall

% Vänsterkurva på gränsen till sladd
N_left(1) = m*g*(B2-mu*H1)/(B1+B2); % Normalkraft på vänster sida [N]
N_right(1) = m*g -N_left(1);             % Normalkraft på höger sida [N]

% Körning rakt fram
N_left(2) = m*g*B2/(B1+B2);         % Normalkraft på vänster sida [N]
N_right(2) = m*g -N_left(2);             % Normalkraft på höger sida [N]

% Högerkurva på gränsen till sladd
N_left(3) = m*g*(B2+mu*H1)/(B1+B2); % Normalkraft på vänster sida [N]
N_right(3) = m*g -N_left(3);             % Normalkraft på höger sida [N]

%% Presentation av resultat

% Skriver ut massecntrums maximalt tillåtna höjd över marken
disp(' ')
disp(['H1max = ' num2str(round(H1max*1e3)) ' mm']);
disp(' ')

% Graf som visar vmax som funktion av kurvradien
figure
plot(R,v*3.6)
xlabel('Kurvradie R [m]')
ylabel('Gränshastighet för sladd v_m_a_x [km/h]')
grid on

% Tar fram lastfördelningen mellan sidorna
Last_left = zeros(3,1); Last_right = Last_left;
for i = 1:3
    Last_left(i) = N_left(i)/(m*g);
    Last_right(i) = N_right(i)/(m*g);
end

% Tabell med lastfördelning
Fall_sida ={'Vänster kurva   ' 'Rakt            ' 'Höger kurva     '};
disp(['               ' 'Vänster   ' 'Höger'])
for i =1:3
    disp([Fall_sida{i} num2str(round(Last_left(i)*100)) '%      ' num2str(round(Last_right(i)*100)) '% '])
end
disp(' ')

disp(['               ' 'Vänster   ' 'Höger'])
for i =1:3
    disp([Fall_sida{i} num2str(round(N_left(i))) ' N     ' num2str(round(N_right(i))) ' N'])
end
disp(' ')

disp('- - - - - - - - - - - - - - - - - - - - - - - - - - -')
%%

% Tar fram normalkrafter mellan framdäck och vägbana, dessa krafter
% blir de som kommer verka i lagren??

% Skapar vektorer
Last_front_left = zeros(5,1); Last_front_right = Last_front_left;

% Lastfördelning fram för Vänster sida
Last_front_left(1) = Last_fram(2)*Last_left(1);
Last_front_left(2) = Last_fram(1)*Last_left(2);
Last_front_left(3) = Last_fram(2)*Last_left(2);
Last_front_left(4) = Last_fram(3)*Last_left(2);
Last_front_left(5) = Last_fram(2)*Last_left(3);

%Lastfördelning fram för Höger sida
Last_front_right(1) = Last_fram(2)*Last_right(1);
Last_front_right(2) = Last_fram(1)*Last_right(2);
Last_front_right(3) = Last_fram(2)*Last_right(2);
Last_front_right(4) = Last_fram(3)*Last_right(2);
Last_front_right(5) = Last_fram(2)*Last_right(3);

% Krafter fram
Kraft_left = zeros(5); Kraft_right = Kraft_left;
for i = 1:5
    Kraft_right(i) = Last_front_right(i)*m*g;
    Kraft_left(i) = Last_front_left(i)*m*g;
end



% Visar tabell över total lastfördelning fram
Fall_fram ={'Vänstersväng        ' 'Acceleration        ' 'Konstant hastighet  ' 'Inbromsning         ' 'Högersväng          '};
disp(['                    ' 'Vänster   ' 'Höger'])
for i =1:5
    disp([Fall_fram{i} num2str(round(Last_front_left(i)*100)) '%      ' num2str(round(Last_front_right(i)*100)) '% '])
end
disp(' ')

% Visar tabell över totala normalkrafter fram
disp(['                    ' 'Vänster   ' 'Höger'])
for i =1:5
    disp([Fall_fram{i} num2str(round(Kraft_left(i))) ' N     ' num2str(round(Kraft_right(i))) ' N'])
end
disp(' ')

disp('- - - - - - - - - - - - - - - - - - - - - - - - - - -')

%% Beräkning av lagerkrafter ?? och livslängder

% Data för däck
v_kon = 50/3.6;
r = 0.254;
omkrets = pi*2*r;
n_kon = 60*(v_kon/omkrets);

% Krav

Min_s = 3000;                       % Minimala sträckan gokarten ska gå [mil]
Min_n = Min_s*10000*10^-6/omkrets   % Minimalt antal miljoner 
                                    % varv hjul/lager snurrar för sträckan

U1 = 0.25;      % Andel den körs i konstant hastighet i högerkurva
U2 = 0.25;      % Andel den körs i konstant hastighet i vänsterkurva
U3 = 0.25;      % Andel den körs i maximal acceleration på raksträcka
U4 = 0.10;      % Andel den körs i konstant hastighet på raksträcka
U5 = 0.15;      % Andel av sträckan den bromsas in från 50km/h

U_vec = [U1 U2 U3 U4 U5];

% Lagerdata 629 spårkullager
C0 = 1960;
C = 4750;
f0 = 12;
p = 3;
d = 0.009;
D = 0.026;

%% Olika fall för körning rakt fram

% Livslängd för främre lager vid konstant hastighet [vänster lager]
Fa_kon = 0;
Fr_kon = Kraft_left(3);
P_kon = Ekv_kraft_fun(Fr_kon,Fa_kon,f0,C0);
L10_kon = L10_fun(C,p,P_kon);

% Livsläng för främre lager vid acceleration [vänster lager]
Fa_acc = 0;
Fr_acc = Kraft_left(2);
P_acc = Ekv_kraft_fun(Fr_acc,Fa_acc,f0,C0);
L10_acc = L10_fun(C,p,P_acc);

% Livslängd för främre lager vid inbromsning [vänster lager]
Fa_broms = 0;
Fr_broms = Kraft_left(4);
P_broms = Ekv_kraft_fun(Fr_broms,Fa_broms,f0,C0);
L10_broms = L10_fun(C,p,P_broms);


%% Olika fall för svängning

R = 6;
v_turn = 24.7/3.6;

% Livslängd för frömre lager vid vänstersväng [vänster lager]
Fa_turn_l = (m*v_turn^2)/R;
Fr_turn_l = Kraft_left(1);
P_turn_l = Ekv_kraft_fun(Fr_turn_l,Fa_turn_l,f0,C0);
L10_turn_l = L10_fun(C,p,P_turn_l);

% Livslängd för frömre lager vid högersväng [vänster lager]
Fa_turn_r = (m*v_turn^2)/R;
Fr_turn_r = Kraft_left(5);
P_turn_r = Ekv_kraft_fun(Fr_turn_r,Fa_turn_r,f0,C0);
L10_turn_r = L10_fun(C,p,P_turn_r);



%% Total livslängd

L10_vec = [L10_turn_l L10_turn_r L10_acc L10_kon L10_broms];

L10_tot = L10_tot_fun(L10_vec,U_vec)




%%
function L10 = L10_fun(C,p,P)
    L10 = (C/P)^p;
end

function L10h = L10h_fun(C,p,P,n)
    L10 = L10_fun(C,p,P);
    L10h = L10*(10^6)/(60*n);
end

function L10_tot = L10_tot_fun(L10_vec,U_vec)
    kvot = zeros(length(L10_vec),1);

    for i = 1:length(L10_vec)
        kvot(i) = U_vec(i)/L10_vec(i);
    end

    L10_tot = 1/sum(kvot);
end

function P = Ekv_kraft_fun(Fr,Fa,f0,C0)
    e_vec = [0.19 0.22 0.26 0.28 0.3 0.34 0.38 0.42 0.44];

    if Fa/Fr < max(e_vec)
        P = Fr;
    else
        [X,Y] = interpol_fun(Fa,f0,C0);
        P = X*Fr + Y*Fa;
    end

end

function [X,Y] = interpol_fun(Fa,f0,C0)

    x = (f0*Fa/C0);
    
    Y_vec = [2.3 1.99 1.71 1.55 1.45 1.31 1.15 1.04 1];
    x_vec = [0.172 0.345 0.689 1.03 1.38 2.07 3.45 5.17 6.89];
    X = 0.56;

    x_vec(10) = x;
    x_vec = sort(x_vec);
    i = find(x_vec==x);

    Y = Y_vec(i) + (Y_vec(i-1)-Y_vec(i))*(x_vec(i+1)-x_vec(i))/(x_vec(i+1)-x_vec(i-1));

end


































 