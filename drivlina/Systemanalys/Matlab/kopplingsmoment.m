function MK = kopplingsmoment(omega1)
% Funktionsfil för beräkning av maximalt överförbart kopplingsmoment som
% funktion av vinkelhastigheten hos medbringaren.
%
% OBS! Ekvationerna i funktionen ska bytas ut när ni kommer så
% långt att ni har tagit fram en modell som beskriver er egen koppling.
%
% Anders Söderberg, KTH Maskinkonstuktion, 2018-08-31

% b(1) =  -3.392;
% b(2) =  2.148e-004;

m=0.185;                        % Massa bromsblock [kg]
mu=0.3;                         % Friktionskoefficient block-trumma [-]
r=0.0331;                       % Tyngddpunktsavstånd bromsblock [m]
R=0.0545;                       % Trumradie [m]
Ffj=80;                         % Förspänningskraft per fjäder [N]

b(1)=-4*Ffj*R*mu;
b(2)=2*m*mu*R*r;

MK = b(1) + b(2)*omega1.^2;     % Maximalt överförbart kopplingsmoment [Nm]
MK = MK.*(MK>0);                % Ingen kontakt mellan block och trumma: MK=0 Nm
end