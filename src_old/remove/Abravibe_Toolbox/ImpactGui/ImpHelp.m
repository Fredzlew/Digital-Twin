D=which('ImpactGui.m');
idx=findstr(D,'\');    % Clear file path from filename
ImpactGuiPath=D(1:idx(end)-1);;

open(fullfile(ImpactGuiPath,'Help\help.htm'))