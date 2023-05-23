function [xp] = VirtualSensVal(data,modeshapes,num_ms,im)

% Total number of sensors
ns = size(data,1);

% Total number of time steps
nt = size(data,2);

% Indicies of predicted DOFs
ip = 1:ns;
for i = 1:length(im)
    ip = ip(ip~=im(i));
end

% Measured data
xm = data(im,:);

% Measured mode shapes used in approximation
phi_m = modeshapes(im,1:num_ms);

% Predicted mode shapes
phi_p = modeshapes(ip,1:num_ms);

% Calculate the psuedo-inverse of the measured mode shapes
phi_minv = pinv(phi_m);

% Calculate the displacements at the predicted locations for each t
xp = zeros(length(ip),nt);
for i = 1:nt
    xp(:,i) = phi_p*phi_minv*xm(:,i);
end
end