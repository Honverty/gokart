function yt = derivataacceleration(t,y)
% Funktionsfil f�r ode-anrop som beskriver acceleration och hastighet f�r
% drivlinan och gokarten vid en acceleration. Funktionen anropar i sin tur
% funktionsfilerna:
%   1. motormoment.m
%   2. kopplingsmoment.m
%   3. lastmoment.m
%
% Anders S�derberg, KTH Maskinkonstuktion, 2018-08-31

% Globala systemparametrar som �r definerade i huvudprogram
global eta u r JM JK1 JK2 JV1 JV2 JL

% Uppdeling av tillst�ndsvektor p� variabler
omega1  = y(1); % Varvtal p� motor [rad/s]
omega2  = y(2); % Varvtal p� kopplingstrumma [rad/s]
omega3  = y(3); % Varvtal p� bakaxel [rad/s]
v       = y(4); % Hastighet [rad/s]
s       = y(5); % F�rdstr�cka [rad/s]

% Ber�kning momentant motormoment [Nm]
MM     = motormoment(omega1);

% Ber�kning av momentant �verf�rt kopplingsmoment [Nm]
MK     = kopplingsmoment(omega1);   % Maximalt �verf�rbart friktionsmoment
if MK > MM
    MK = MM;                        % Hela motormomentet �verf�rs
end

% Ber�kning av momentant lastmoment [Nm]
ML = lastmoment(v);
if ML/(eta*u) > MK
    ML = u*eta*MK;                  % Lasten kan inte driva systemet
end

% Ber�kning av accelerationer och hastigheter
if omega1>omega2 || (omega1== 0 && omega2==0)       % Slirar
    omega1t = (MM-MK)/(JM+JK1);
    omega2t = (MK-ML/(eta*u))/(JK2+JV1+(JV2+JL)/(eta*u^2));
else                                                % Greppar
    omega1t = (MM-ML/(eta*u))/(JM+JK1+JK2+JV1+(JV2+JL)/(eta*u^2));
    omega2t = omega1t;
end
omega3t = omega2t/u;
vt      = omega3t*r;
st      = v;

% Sammanst�llning av vektor med derivator att returnera
yt(1,:)  = omega1t;
yt(2,:)  = omega2t;
yt(3,:)  = omega3t;
yt(4,:)  = vt;
yt(5,:)  = st;


%% EOF






