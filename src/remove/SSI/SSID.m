function [Result]=SSID(output,fs,ncols,nrows,cut)
%Input:
%output: output data of size (No. output channels, No. of data)
%fs: Sampling frequency 
%ncols: The number of columns in hankel matrix (more than 2/3 of No. of data)
%nows: The number of rows in hankel matrix (more than 20 * number of modes)
%cut: cutoff value=2*no of modes

%Outputs :
%Result  : A structure consist of the below components

%Parameters.NaFreq : Natural frequencies vector
%Parameters.DampRatio : Damping ratios vector
%Parameters.ModeShape : Mode shape matrix

%Indicators.EMAC : Extended Modal Amplitude Coherence
%Indicators.MPC  : Modal Phase Collinearity
%Indicators.CMI  : Consistent Mode Indicator
%Indicators.partfac : Participation factor

%Matrices A,C: Discrete A and C matrices
%--------------------------------------------------------------------------

[outputs,npts]=size(output);       % Computes the size of matrix output

if outputs > npts             % Check if output matrix size is proper or should be transposed
    output=output';
    [outputs,npts]=size(output);
end

%--------------------------------------------------------------------------
% Find block sizes

brows=fix(nrows/outputs);     % brows = how many output blocks.
nrows=outputs*brows;          % Redefine the row numbers.
bcols=fix(ncols/1);           % bcols = how many time steps.
ncols=1*bcols;                % Redefine the column numbers.
m=1;                          % inputs number.
q=outputs;                    % outputs number.


%--------------------------------------------------------------------------
% Form the Hankel matriices Yp,Yf.

Yp=zeros(nrows/2,ncols); %Past Output
Yf=zeros(nrows/2,ncols); %Future Output

for ii=1:brows/2
    for jj=1:bcols
        Yp([(ii-1)*q+1:ii*q] ,[(jj-1)*m+1:jj*m] )=output(:,((jj-1)+ii-1)*m+1:((jj)+ii-1)*m); 
    end
end

for ii=brows/2+1:brows
    i=ii-brows/2;
    for jj=1:bcols
        Yf([(i-1)*q+1:i*q] ,[(jj-1)*m+1:jj*m] )=output(:,((jj-1)+ii-1)*m+1:((jj)+ii-1)*m); 
    end
end

%--------------------------------------------------------------------------
% Projection
O=Yf*Yp'*pinv(Yp*Yp')*Yp;



%--------------------------------------------------------------------------
% Decompose the data matrix

[R1,Sigma1,S1]=svd(O,0);
sv=diag(Sigma1);

%--------------------------------------------------------------------------
% Truncate the matrices using the cutoff

D=diag(sqrt(sv(1:cut)));      % build square root of the singular values.
Dinv=inv(D);				  % (sigma)^(-1/2)
Rn=R1(:,1:cut);               % use only the principal eigenvectors
Sn=S1(:,1:cut);               % use only the principal eigenvectors

Obs=Rn*D;                     % Observability matrix

%--------------------------------------------------------------------------
% Calculate the realization of A and find eigenvalues and eigenvectors

A=pinv(Obs(1:nrows/2-q,:))*Obs(q+1:nrows/2,:);  % build A 
C=Obs(1:q,:);            % build C 
clear Yp Yf;

%--------------------------------------------------------------------------
% Extract the modal frequencis , damping ratios and natural frequencis

[Vectors,Values]=eig(A);       % Eigenvalues and Eigenvectors
Lambda=diag(Values);           % roots in the Z-plane
s=log(Lambda).*fs;             % Laplace roots 
zeta=-real(s)./abs(s)*100;     % damping ratios
fd=(imag(s)./2./pi);            
                               % damped natural freqs:
shapes=C*Vectors;              % Mode shapes.

%--------------------------------------------------------------------------
%  Calculate  Modal Participation factors

InvV=inv(Vectors);
partfac=std((InvV*D*Sn')')';           %InvV involved in transforming states to the modal amplitudes by uncoupling A matrix
%--------------------------------------------------------------------------
% Calculate EMAC values

partoutput22=(C*A^(brows/2-1)*Vectors)';               %Predicted Last Block in Obs*Vectors Matrix
partoutput2=(Rn(nrows/2-q+1:nrows/2,:)*D*Vectors)';    %Actual Last Block in Obs*Vectors Matrix (Vectors involved in transforming states to the modal amplitudes by uncoupling A matrixx)
for i=1:1:cut
    
sum1=0;
sum2=0; 
    
for k=1:1:q
            
Rik=abs(partoutput22(i,k))/abs(partoutput2(i,k));

if Rik>1
Rik=1/Rik;
end

Wik=1-(4/pi)*abs(angle(partoutput22(i,k)/partoutput2(i,k)));
            
if Wik<0
Wik=0;
end
            
EMACout=Rik*Wik;
            
EMACijk=EMACout;
sum1=sum1+EMACijk*(abs(partoutput2(i,k)))^2*100;
sum2=sum2+(abs(partoutput2(i,k)))^2;
end
        
EMACC(i)=sum1/sum2;
end

%--------------------------------------------------------------------------
% Calculate MPC values
 
for ii=1:1:cut
    
cd=sum(shapes(:,ii))/q;
cj=shapes(:,ii);

Sxx=(real(cj-cd))'*(real(cj-cd));
Sxy=(real(cj-cd))'*(imag(cj-cd));
Syy=(imag(cj-cd))'*(imag(cj-cd));

mu=(Syy-Sxx)/(2*Sxy);
Beta=mu+sign(Sxy)*sqrt(mu^2+1);
Tau=atan(Beta);

l1=Sxx+(Sxy*(2*(mu^2+1)*(sin(Tau))^2-1))/mu;
l2=Syy-(Sxy*(2*(mu^2+1)*(sin(Tau))^2-1))/mu;

MPC(ii)=100*(2*(l1/(l1+l2)-0.5))^2;

end
%--------------------------------------------------------------------------
% Sort into ascending order

[fd,I]=sort(fd);
zeta=zeta(I);
shapes=shapes(:,I);
s=s(I); 
partfac=partfac(I,:);
EMACC=EMACC(I);
MPC=MPC(I);

%--------------------------------------------------------------------------
% Remove the negative frequencies and frequencies>fs/2

lower=1;
upper=cut;
for ii=1:cut
    if fd(ii) <= 0
        lower=ii+1;
    end
    if fd(cut-ii+1) >= 0.499*fs
        upper=cut-ii;
    end
end

fd1=fd(lower:upper); 
zeta1=zeta(lower:upper);
fd1=fd1./sqrt(1-(zeta1/100).^2); % Calculate the undamped natural frequency
Lambda1=s(lower:upper);
shapes=shapes(:,lower:upper);
partfac=partfac(lower:upper,:);
EMACC1=EMACC(lower:upper);
MPC1=MPC(lower:upper);

% Sort again
%--------------------------------------------------------------------------

[fd1,ii]=sort(fd1); 
zeta1=zeta1(ii); 
shapes=shapes(:,ii); 
EMACC1=EMACC1(ii);
MPC1=MPC1(ii);
shapes=shapes(:,ii); 
partfac=partfac(ii,:);
Lambda1=Lambda1(ii);

% Normalize modeshapes and process them
%--------------------------------------------------------------------------
Phi=shapes;
HH = size(Phi);

[C1,II] = max(abs(Phi));
    for jj = 1:HH(2)
        b = -angle(Phi(II(jj),jj));
        ModeShapeS(:,jj) = real(Phi(:,jj)*exp(1i*b));
        ModeShapeS(:,jj) = ModeShapeS(:,jj)/norm(ModeShapeS(:,jj));
    end
shapes=ModeShapeS;

%--------------------------------------------------------------------------
EMAC=EMACC1';
MPC=MPC1';
CMI=EMAC.*MPC/100;

%--------------------------------------------------------------------------
Result.Parameters.NaFreq=fd1;
Result.Parameters.DampRatio=zeta1;
Result.Parameters.ModeShape=shapes;

Result.Indicators.EMAC=EMAC;
Result.Indicators.MPC=MPC;
Result.Indicators.CMI=CMI;
Result.Indicators.partfac=partfac;

Result.Matrices.A=A;
Result.Matrices.C=C;

end