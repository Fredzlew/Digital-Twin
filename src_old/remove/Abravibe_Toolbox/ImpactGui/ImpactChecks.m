function [status,OutTrigIdx,idx] = ImpactChecks(x,TrigIdx,N)
% IMPACTCHECKS  Checks that combination of triggers and blocksize do not
% exahust data
%
%       [status,OutTrigIdx,Idx] = ImpactChecks(x,TrigIdx,N);
%
%       status=0 means everything is okay
%       status=-1 means one or more blocks overlap
%       status=-2 means TrigIdx is empty after editing
%
%       TrigIdx         New indeces that do not cause overlap, and are not 
%                       exhausting data
%       Idx             Index into TrigIdx which cause blocks to overlap
%                       Only nonempty if status == -1
%
% This is an internal function for ImpactGui.m. If it finds trigger indeces 
% that exhaust data, or make blocks overlap, those triggers are removed

% Copyright (c) 2009-2014 by Anders Brandt
% Email: abra@iti.sdu.dk
% Version: 1.0 2014-07-15
% This file is part of ABRAVIBE Toolbox for NVA

status=0;

% First check if any index plus blocksize is outside x
% If so, remove that/those indeces
idx=find(TrigIdx+N-1 > length(x));
if ~isempty(idx)
    if length(TrigIdx) > length(idx)
        OutTrigIdx=TrigIdx(1:idx(1)-1);
    else
        OutTrigIdx='';
        status = -2;
    end
end

% Next, check if any blocks overlap
D=diff(TrigIdx);
idx=find(D <= N);
if ~isempty(idx)
    status = -1;
    OutTrigIdx=setxor(TrigIdx,TrigIdx(idx));
else
    OutTrigIdx=TrigIdx;
end


