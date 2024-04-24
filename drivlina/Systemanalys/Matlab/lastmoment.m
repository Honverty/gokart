function ML = lastmoment(v)
% Funktionsfil f�r ber�kning av lastmoment p� motoraxeln som funktion av
% gokartens hastighet i f�rdrikningen. Modellen tar h�nsyn till
% rullmotst�nd och luftmotst�nd.
%
% Anders S�derberg, KTH Maskinkonstuktion, 2018-08-31

global f m g c rho A r

Frull  = f*m*g;                 % Rullmotst�nd [N]
Fluft  = 0.5*c*rho*A*v.^2;      % Luftmotst�nd [N]
ML     = (Frull + Fluft)*r;     % Lastmoment [Nm]