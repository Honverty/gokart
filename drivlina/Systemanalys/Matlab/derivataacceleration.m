function yt = derivataacceleration(t,y)
% Funktionsfil för ode-anrop som beskriver acceleration och hastighet för
% drivlinan och gokarten vid en acceleration. Funktionen anropar i sin tur
% funktionsfilerna:
%   1. motormoment.m
%   2. kopplingsmoment.m
%   3. lastmoment.m
%
% Anders Söderberg, KTH Maskinkonstuktion, 2018-08-31

% Globala systemparametrar som är definerade i huvudprogram
global eta u r JM JK1 JK2 JV1 JV2 JL

% Uppdeling av tillståndsvektor på variabler
omega1  = y(1); % Varvtal på motor [rad/s]
omega2  = y(2); % Varvtal på kopplingstrumma [rad/s]
omega3  = y(3); % Varvtal på bakaxel [rad/s]
v       = y(4); % Hastighet [rad/s]
s       = y(5); % Färdsträcka [rad/s]

% Beräkning momentant motormoment [Nm]
MM     = motormoment(omega1);

% Beräkning av momentant överfört kopplingsmoment [Nm]
MK     = kopplingsmoment(omega1);   % Maximalt överförbart friktionsmoment
if MK > MM
    MK = MM;                        % Hela motormomentet överförs
end

% Beräkning av momentant lastmoment [Nm]
ML = lastmoment(v);
if ML/(eta*u) > MK
    ML = u*eta*MK;                  % Lasten kan inte driva systemet
end

% Beräkning av accelerationer och hastigheter
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

% Sammanställning av vektor med derivator att returnera
yt(1,:)  = omega1t;
yt(2,:)  = omega2t;
yt(3,:)  = omega3t;
yt(4,:)  = vt;
yt(5,:)  = st;


%% EOF






