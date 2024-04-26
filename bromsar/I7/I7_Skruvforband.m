%% Produktutveckling
clear all;
clc
format long
%% I7 - Skruvf�rband

p = 13*10^6; %[Pa] �vertrycket i r�ret
D = 135*10^(-3); %[m] Diametern p� r�ret
S_a = 12; %[1] Antal skruvar

%M10 ger
d_y = 10*10^(-3); %[m] skruvens ytterdiameter
d_m = 9.03*10^(-3); %[m] skruvens medeldiameter
d_i = 8.38*10^(-3); %[m] skruvens innerdiameter
N = 17*10^(-3); %[m] Skruvens nyckelvidd
D_h = 11*10^(-3); %[m] Skruvens h�ldiameter
stigning = 1.5*10^(-3); %[m] Stigning
%F_s = 37.1*10^3; %[N] Skruvens str�ckkraft

%H�llfasthetsklass 10.9 ger
sigma_b = 10*100*10^6; %[Pa] Skruvens brottgr�ns
sigma_s = 0.9*sigma_b; %[Pa] Skruvens str�ckgr�ns

L = 6*d_y; %[m] skruvens kl�ml�ngd

%a
A = pi*(D/2)^2; %[m^2] tv�rsnittsarean p� r�ret
F_p = p*A %[N] Kraften som drar is�r fl�nshalvorna

%b
A_s = pi*(d_y/2)^2; %[m^2] medelarea f�r skruvens tv�rsnitt
E_s = 206*10^9; %[N/m^2] Skruvens elastisitetsmodul
k_s = A_s * E_s / L %[N/m] Skruvens fj�derkonstant

%c
E_f = E_s; %[N/m^2] Fl�nsens elastisitetsmodul
A_f = pi/4 * ( (N+0.3*L)^2 - D_h^2) %[m^2] Fl�nsarean
k_f = A_f * E_f / L %[N/m] Fl�nsens fj�derkonstant

%d
sigma = 0.8*sigma_s; %[Pa] Sp�nningen i en skruv
A_i = pi * (d_i/2)^2; %[m^2] skruvens innerarea
F_i = sigma*A_i %[N] F�rsp�nningskraften i en skruv

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

%H�rifr�n har jag f�tt fel, bara s� du vet.
%test! Testad mot uppg.9 i tenta 230607 och f�r r�tt svar. 
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
