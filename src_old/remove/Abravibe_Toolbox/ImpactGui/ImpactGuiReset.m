function ImpactGuiReset
% ImpactGuiReset    Deletes the file ImpactGuiSave.mat
%
% Execute this command if ImpactGui does not start due to a conflict with
% a directory that has been moved or deleted

warning off
F=which('ImpactGuiSave.mat');
delete(F)
warning on