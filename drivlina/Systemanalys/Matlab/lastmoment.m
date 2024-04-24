function ML = lastmoment(v)
% Funktionsfil för beräkning av lastmoment på motoraxeln som funktion av
% gokartens hastighet i färdrikningen. Modellen tar hänsyn till
% rullmotstånd och luftmotstånd.
%
% Anders Söderberg, KTH Maskinkonstuktion, 2018-08-31

global f m g c rho A r

Frull  = f*m*g;                 % Rullmotstånd [N]
Fluft  = 0.5*c*rho*A*v.^2;      % Luftmotstånd [N]
ML     = (Frull + Fluft)*r;     % Lastmoment [Nm]