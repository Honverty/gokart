%% Produktutveckling
clear all;
clc
format long
%% I7 - Skruvförband

p = 13*10^6; %[Pa] Övertrycket i röret
D = 135*10^(-3); %[m] Diametern på röret
S_a = 12; %[1] Antal skruvar

%M10 ger
d_y = 10*10^(-3); %[m] skruvens ytterdiameter
d_m = 9.03*10^(-3); %[m] skruvens medeldiameter
d_i = 8.38*10^(-3); %[m] skruvens innerdiameter
N = 17*10^(-3); %[m] Skruvens nyckelvidd
D_h = 11*10^(-3); %[m] Skruvens håldiameter
stigning = 1.5*10^(-3); %[m] Stigning
%F_s = 37.1*10^3; %[N] Skruvens sträckkraft

%Hållfasthetsklass 10.9 ger
sigma_b = 10*100*10^6; %[Pa] Skruvens brottgräns
sigma_s = 0.9*sigma_b; %[Pa] Skruvens sträckgräns

L = 6*d_y; %[m] skruvens klämlängd

%a
A = pi*(D/2)^2; %[m^2] tvärsnittsarean på röret
F_p = p*A %[N] Kraften som drar isär flänshalvorna

%b
A_s = pi*(d_y/2)^2; %[m^2] medelarea för skruvens tvärsnitt
E_s = 206*10^9; %[N/m^2] Skruvens elastisitetsmodul
k_s = A_s * E_s / L %[N/m] Skruvens fjäderkonstant

%c
E_f = E_s; %[N/m^2] Flänsens elastisitetsmodul
A_f = pi/4 * ( (N+0.3*L)^2 - D_h^2) %[m^2] Flänsarean
k_f = A_f * E_f / L %[N/m] Flänsens fjäderkonstant

%d
sigma = 0.8*sigma_s; %[Pa] Spänningen i en skruv
A_i = pi * (d_i/2)^2; %[m^2] skruvens innerarea
F_i = sigma*A_i %[N] Förspänningskraften i en skruv

%e
F_smax = F_i + F_p/S_a*(1/(1+k_f/k_s))

%f
F_a = F_i;
tanalfa = stigning / (pi*d_m);
%alfa = atan(tanalfa);
beta = pi/3;
%mu_g = 0.12/cos(beta/2);
%mu_g = 0.12;
%theta = atan(cos(alfa)*tan(beta/2));
%mu_g = 0.12/cos(theta);
mu_u = 0.11;
D = (D_h + N)/2;
M1 = F_a*d_m/2 * ( (tanalfa + mu_g)/(1-mu_g*tanalfa) + mu_u*D/d_m)

%Härifrån har jag fått fel, bara så du vet.
%test! Testad mot uppg.9 i tenta 230607 och får rätt svar. 
%F_i = 24985;
%F_a = F_i;
%d_m = 7.19*10^(-3);
%stigning = 1.25 * 10^(-3);
%tanalfa = stigning / (pi*d_m);
%mu_g = 0.15/cos(beta/2);
%mu_u = 0.15;
%D_h = 9*10^(-3);
%N = 13*10^(-3);
%D = (D_h + N)/2;

%M1test = F_a*d_m/2 * ( (tanalfa + mu_g)/(1-mu_g*tanalfa) + mu_u*D/d_m)


%g)
F_a = F_i;
M2_g = F_a*d_m/2 * ( (tanalfa - mu_g)/(1+mu_g*tanalfa) - mu_u*D/d_m)

%h)
F_a = F_smax;
M2_h = F_a*d_m/2 * ( (tanalfa - mu_g)/(1+mu_g*tanalfa) - mu_u*D/d_m)
