function MM = motormoment(omega1)

% Motorkoefficienter
a(1) = -1654;
a(2) =    24.35;
a(3) =    -0.0256;

MM = a(1)./omega1 + a(2) + a(3)*omega1; % Drivande motormoment [Nm]