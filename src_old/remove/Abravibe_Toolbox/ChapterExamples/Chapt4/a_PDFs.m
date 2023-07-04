% a_PDFs    Probability Density Function examples
%           This example creates some Gaussian and non-Gaussian data, and then 
%           calculates the PDF of the data, overlaid with the PDF of true Gaussian
%           data with the same mean and standard deviation.


% This file is part of the examples for the ABRAVIBE Toolbox for NVA which 
% is an accompanying toolbox for the book
% Brandt, Anders: "Noise and Vibration Analysis: Signal Analysis and
% Experimental Procedures," Wiley 2011. ISBN: 13-978-0-470-74644-8.
% Copyright 2011, Anders Brandt.


close all
clear

% Some settings to easily change plot appearance
%--------------------------------------------------
FontSize=9;
FontName='Times New Roman';
LineWidth=1;
LineType={'-k','--k','-.k',':k'};
%--------------------------------------------------

%================================================================
% Gaussian data
N=40;
x=randn(100000,1);
[Pdf, XAx, G, mu, sigma] = apdf(x, N, 1);
title('Gaussian Data')
[Pdf, XAx, G, mu, sigma] = apdf(x, N, 2);
title('Gaussian Data')

% non-Gaussian
y=x+0.1*x.*abs(x);
[Pdf, XAx, G, mu, sigma] = apdf(y, N, 1);
title('non-Gaussian Data')
[Pdf, XAx, G, mu, sigma] = apdf(y, N, 2);
title('non-Gaussian Data')
