function [H,f,C,Tff, ImpactNos] = Imp2FrfGui(x,y,fs,N,TrigIdx,DIdx,FWinLength,ExpWinPar,Type,flo,fhi)
% IMP2FRFGUI  Calculate FRF(s) from impact time data recording (ImpactGui version)
%
%        [H,f,C,Tff, ImpactNos] = imp2frf2(x,y,fs,N,TrigIdx,DIdx,FWinLength,FWinLength,ExpWinPar,Type,flo,fhi)
%
%           H           Frequency response(s), N/2+1-by-D
%           f           Frequency axis, N/2+1-by-1
%           C           Coherence function(s), N/2+1-by-D
%           Tff         Transient spectrum of force signal
%           ImpactNos   For automatic processing types this parameter
%                       contains a cell array with vectors of the impact 
%                       numbers that was used to process the FRFs
%
%           x           Force vector
%           y           Response signal(s) in D column(s)
%           fs          Sampling frequency
%           N           Block size for FFT
%           TrigIdx     Trigger indeces from imptrig
%           DIdx        Double impact indeces from imptrig
%           FWinLength  Force window in percent of N
%           ExpWinPar   End value of exponential window in percent
%           Type        Process type can be:
%                       FrequencyManual         % 'Normal' impact testing
%                       TimeSynchronousManual   % Averaging is done in time domain
%                       BestCoherenceOptimization % Automatic mode 1
%           flo         Low frequency limit for optimization methods
%           fhi         High frequency limit for optimization methods
%
% This is an internal, special version of imp2frf for use in ImpactGui.
% It does not prompt for user selection but uses predefined trigger events.
%
% See also IMPSETUP IMPMASSCAL IMPPROC IMPTRIG AFORCEW AEXPW
%

% Copyright (c) 2009-2011 by Anders Brandt
% Email: abra@iti.sdu.dk
% Version: 1.0 2014-07-06
% This file is part of ABRAVIBE Toolbox for NVA

% References
% Brandt, A. "Noise and Vibration Analysis - Signal Analysis and
% Experimental Procedures," John Wiley and Sons, 2011
% Specially for the Best coherence optimization: This method finds the two
% impacts that give the best coherence.
% Brandt, A. & Brincker, R. Impact Excitation Processing for Improved 
% Frequency Response Quality Proc. 28th International Modal Analysis 
% Conference, Jacksonville, FL, 2010

% Check consistency
% Find lengths and dimensions
D=length(y(1,:));
if length(x(1,:)) > 1
    x=x(:,1);
end
if length(y(1,:)) > 1
    D=length(y(1,:));       % Number of responses
end

ImpactNos=[];

if strcmp(upper(Type),'FREQUENCYMANUAL')
    % Use the selected impacts to calculate FRF etc for each response (column)
    % in y
    yin=y; clear y
    for k = 1:D
        y=yin(:,k);
        L=N/2+1;
        Tyf=zeros(L,1);
        Tyy=zeros(L,1);
        Tff=zeros(L,1);
        fw=aforcew(N,FWinLength);
        ew=aexpw(N,ExpWinPar);
        for n=1:length(TrigIdx)
            F=x(TrigIdx(n):TrigIdx(n)+N-1);
            F=F-F(1);
            F=F.*fw.*ew;
            Ff=fft(F);
            Ff=Ff(1:L);
            Y=y(TrigIdx(n):TrigIdx(n)+N-1);
            Y=Y.*ew;
            Yf=fft(Y);
            Yf=Yf(1:L);
            Tyf=Tyf+(Yf.*conj(Ff));
            Tyy=Tyy+abs(Yf).^2;
            Tff=Tff+abs(Ff).^2;
        end
        H(:,k)=Tyf./Tff;
        f=(0:fs/N:fs/2)';
        if length(TrigIdx) == 1     % Undefined coherence
            C=0;
        else
            C(:,k)=abs(Tyf).^2./Tff./Tyy;
        end
    end
    Tff=Tff/n/fs;               % Scaled as transient average spectrum
elseif strcmp(upper(Type),'TIMESYNCHRONOUSMANUAL')
    % Use the selected impacts to calculate time synchronous averages of
    % each signal (force, and each response). Then produce FRF by simply
    % dividing spectra
    xavg=zeros(N,1);
    yavg=zeros(N,D);
    for n=1:length(TrigIdx)
        xtemp=x(TrigIdx(n):TrigIdx(n)+N-1);
%         xtemp=xtemp-xtemp(1);       % Remove value at time 0
        xavg=xavg+xtemp;
        ytemp=y(TrigIdx(n):TrigIdx(n)+N-1,:);
%         ytemp=ytemp-ones(N,1)*ytemp(1,:);
        yavg=yavg+ytemp;
    end
    xavg=xavg/length(TrigIdx);
    yavg=yavg/length(TrigIdx);
    % Now produce FFT results and compute outputs:
    fw=aforcew(N,FWinLength);
    ew=aexpw(N,ExpWinPar);
    L=N/2+1;
    F=xavg.*fw.*ew;        % Force and exponential window on force, still
    Ff=fft(F);
    Ff=Ff(1:L);
    y=yavg.*(ew*ones(1,D));
    Yf=fft(y);
    Yf=Yf(1:L,:);
    H=Yf./(Ff*ones(1,D));
    f=(0:fs/N:fs/2)';
    C=[];                         % Undefined coherence
    Tff=abs(Ff)/fs;               % Scaled as transient average spectrum
elseif strcmp(upper(Type),'BESTTWO')
    minCSum=1e10;
    minPair=[0 0];
    CSumIdx=1;
    IdxPair=[];
    % Find the pair of impacts that gives the "best coherence" (i.e.
    % minimum df(1-g^2)
    for ch=1:D;
        for l=1:length(TrigIdx)
            Tidx1=TrigIdx(l);
            for k=l+1:length(TrigIdx)
                Tidx2=TrigIdx(k);
                [Ht,ft,Ct] = imp2frf2(x,y(:,ch),fs,N,TrigIdx([l k]),0,FWinLength,ExpWinPar,0);
                idx1=min(find(ft>=flo));
                idx2=min(find(ft>=fhi));
                fidx=idx1:idx2;
                CSum(CSumIdx)=sum(1-Ct(fidx));
                if CSum(CSumIdx) < minCSum;
                    minCSum=CSum(CSumIdx);
                    minPair=[Tidx1 Tidx2];
                end
                CSumIdx=CSumIdx+1;
            end
        end
        idx=find(TrigIdx~=minPair(1) & TrigIdx~=minPair(2));
        for k=1:length(idx)
            if ismember(minPair,TrigIdx(k)) == zeros(size(minPair))
                Triggers=[minPair TrigIdx(k)];
                [Ht,ft,Ct,Tfft] = imp2frf2(x,y(:,1),fs,N,Triggers,0,FWinLength,ExpWinPar,0);
                CSum=sum(1-Ct(fidx));
                if CSum < minCSum;
                    minCSum=CSum;
                    minPair=Triggers;
                end
            end
        end
        Triggers=minPair;           % Better name
        ImpactNos{ch} = find(ismember(TrigIdx,Triggers));
        [H(:,ch),f,C(:,ch),Tff] = imp2frf2(x,y(:,ch),fs,N,Triggers,0,FWinLength,ExpWinPar,0);
    end
end

f=f(:);