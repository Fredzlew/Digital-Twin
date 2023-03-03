function abrainst
% ABRAINST  Run once during installation, after unpacking zip files, 
%           to update path list etc.
%
% This function also 'compensates' for differences between MATLAB and GNU
% Octave.
%
% NOTE! This command MUST be run from the directory where this file is
% located, i.e. you should run something similar to (inside MATLAB)
% >> cd c:\program files\matlab\R2015b\toolbox\abravibe_toolbox
% >> abrainst

% Copyright (c) 2009-2018 by Anders Brandt
% Email: abra@iti.sdu.dk
% Version: 1.0 2011-06-23   
%          1.1 2018-03-09 Changed checksw procedure as very slow
%
% This file is part of ABRAVIBE Toolbox for NVA


% Check if run from MATLAB or GNU/OCTAVE
v=ver;
sw=0;
for n = 1:length(v)
    Name=v(n).Name;
    if strcmp(upper(Name),'MATLAB')
        sw='MATLAB';
        save(strcat('abravibe',filesep,'Platform'),'sw');
    elseif strcmp(upper(Name),'OCTAVE')
        sw='OCTAVE';
        save(strcat('abravibe',filesep,'Platform'),'sw');
    end
end
if sw == 0
    error('Something went wrong establishing the platform (MATLAB or OCTAVE)!')
    error('Installation aborted')
end

% Check installation directory
InstDir=pwd;

if strcmp(sw,'MATLAB')
    addpath(strcat(InstDir,filesep,'abravibe'));
    addpath(strcat(InstDir,filesep,'animate'));
    addpath(strcat(InstDir,filesep,'ImpactGui'));
    savepath
elseif strcmp(sw,'OCTAVE')
	fprintf('Currently you need to install manually for Octave. See ABRAVIBE user''s manual');
end
