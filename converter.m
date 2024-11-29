function [Theta1, Theta2] = convert_cartesian_to_torus(X, Y, Z, R, r)
% convert_cartesian_to_torus: Converts Cartesian coordinates (X, Y, Z) to toroidal coordinates (Theta1, Theta2)
% X, Y, Z: Cartesian coordinates
% R: Major radius of the torus (distance from center of hole to center of tube)
% r: Minor radius of the torus (radius of the tube)

% Calculate Theta1 from X and Y
Theta1 = atan2(Y, X);

% Calculate the radial distance from the torus center
radial_distance = sqrt(X.^2 + Y.^2) - R;

% Calculate Theta2 from radial distance and Z
Theta2 = atan2(Z, radial_distance);
end

% Example usage:
% Define parameters
R = 0.05; % Major radius
r = 0.1;  % Minor radius

% Define Cartesian coordinates
X = 0.0614807;
Y = 0.130859;
Z = 0.0324;

% Convert to toroidal coordinates
[Theta1, Theta2] = convert_cartesian_to_torus(X, Y, Z, R, r);

% Display the result
fprintf('Theta1: %f radians\n', Theta1);
fprintf('Theta2: %f radians\n', Theta2);
