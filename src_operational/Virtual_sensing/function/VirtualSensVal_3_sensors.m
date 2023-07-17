function [xp,qt] = VirtualSensVal_3_sensors(data,modeshapes,num_ms,im,ls)
% Inputs:
% data = data matrix containing all measurements for all locations
% modeshapes = containing all mode shapes for all locations
% num_ms = number of mode shapes included in approximation <= length(im)
% im = index of measured locations
% ls = location of sensors
% Output:
% xp = displacements at predicted locations

% Total number of sensors
ns = size(data,1);

% Total number of time steps
% nt = size(data,2);

% Indicies of predicted DOFs
% ip = 1:ns;
% for i = 1:length(im)
%     ip = ip(ip~=im(i));
% end

% Measured data
xm = data(ls,:);

% Measured mode shapes used in approximation
phi_m = modeshapes(im,num_ms);

% Predicted mode shapes
phi_p = modeshapes(1:ns,num_ms);

% Calculate the psuedo-inverse of the measured mode shapes
phi_minv = (phi_m'*phi_m)^-1*phi_m';
% phi_minv = phi_m'*inv(phi_m*phi_m')
% phi_minv = pinv(phi_m)

% Calculate modal coordinates
qt = phi_minv*xm;

% Calculate the displacements at the predicted locations for each t
% xp = zeros(ns,nt);
% for i = 1:nt
%     xp(:,i) = phi_p*phi_minv*xm(:,i);
% end
xp = phi_p*qt;

end