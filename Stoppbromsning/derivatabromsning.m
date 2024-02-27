function yt = derivatabromsning(t,y)
% Funktionsfil för ode-anrop som beskriver acceleration och hastighet för 
% gokarten vid en inbromsing med låst skivbroms
%
% Anders Söderberg, KTH Maskinkonstuktion, 2018-08-31 

% Globala systemparametrar som är definerade i huvudprogram
global  L1 L2 H1 H2 m g f c rho A mu

% Uppdeling av tillståndsvektor på variabler
v       = y(1); % Hastighet [m/s]
s       = y(2); % Färdsträcka [m]

% Beräknar krafter som verkar på systemet
if v>0
    Frull  = f*m*g;                                     % Rullmotstånd
    Fluft  = 0.5*c*rho*A*v^2;                           % Luftmotstånd
    N2     = (m*g*L1-Frull*H1+Fluft*H2)/(L1+L2+mu*H1);  % Normalkraft mellan bakdäck och vägbana
    F      = -mu*N2;                                    % Bromsande friktionskraft mellan bakdäck och vägbana
else
    Frull = 0;
    Fluft = 0;
    N2    = m*g*L1/(L1+L2);
    F     = 0;
end

% Bestämmer tidderivator av tillståndsvariabler
vt      = (1/m)*(F-Fluft-Frull);                    % Acceleration [m/s^2]
st      = v;                                        % Hastighet [m/s^2]

% Sammanställning av vektor med derivator att returnera
yt(1,:)  = vt;
yt(2,:)  = st;

%% EOF