function J=costfun(omegaNUM,omegaOMA)
% omegasNUM: numerical natural frequencies
% omegasOMA: OMA natural frequencies
J = sum((omegaNUM-omegaOMA).^2);
end