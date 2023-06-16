function mAc=Mac(Phi1,Phi2)
% This function is to compute and plot Modal Assurance Criterion (MAC) matrix between identified mode shapes
% rectangle around the peaks.
% Programmer: Mohammad Farshchin, Ph.D candidate at The UofM
% Email: Mohammad.Farshchin@gmail.com
% Last modified: 9/23/2015

% This function calculates mac between phi1 and phi2
mAc= (abs(Phi1'*Phi2))^2/((Phi1'*Phi1)*(Phi2'*Phi2));
end