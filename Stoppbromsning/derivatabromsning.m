function yt = derivatabromsning(t,y)
% Funktionsfil f�r ode-anrop som beskriver acceleration och hastighet f�r 
% gokarten vid en inbromsing med l�st skivbroms
%
% Anders S�derberg, KTH Maskinkonstuktion, 2018-08-31 

% Globala systemparametrar som �r definerade i huvudprogram
global  L1 L2 H1 H2 m g f c rho A mu

% Uppdeling av tillst�ndsvektor p� variabler
v       = y(1); % Hastighet [m/s]
s       = y(2); % F�rdstr�cka [m]

% Ber�knar krafter som verkar p� systemet
if v>0
    Frull  = f*m*g;                                     % Rullmotst�nd
    Fluft  = 0.5*c*rho*A*v^2;                           % Luftmotst�nd
    N2     = (m*g*L1-Frull*H1+Fluft*H2)/(L1+L2+mu*H1);  % Normalkraft mellan bakd�ck och v�gbana
    F      = -mu*N2;                                    % Bromsande friktionskraft mellan bakd�ck och v�gbana
else
    Frull = 0;
    Fluft = 0;
    N2    = m*g*L1/(L1+L2);
    F     = 0;
end

% Best�mmer tidderivator av tillst�ndsvariabler
vt      = (1/m)*(F-Fluft-Frull);                    % Acceleration [m/s^2]
st      = v;                                        % Hastighet [m/s^2]

% Sammanst�llning av vektor med derivator att returnera
yt(1,:)  = vt;
yt(2,:)  = st;

%% EOF