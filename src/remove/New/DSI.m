function [Result]=DSI(output,input,fs,cut)

%Input:
%output: output data of size (No. output channels, No. of data)
%input: input data of size (No. input channels, No. of data)
%fs: Sampling frequency 
%cut: cutoff value=2*no of modes

%Outputs :
%Result  : A structure consist of the below components
%Parameters :  NaFreq : Natural frequencies vector
%DampRatio : Damping ratios vector
%ModeShape : Mode shape matrix
%Matrices A,B,C,D: Discrete A,B,C and D matrices

%--------------------------------------------------------------------------
Dt=1/fs;
DAT = iddata(output',input',Dt);     % Input Output data
TH1 = n4sid(DAT,cut,'Disturbance','none');%DSI

A = TH1.A;                           % A matrix
B = TH1.B;                           % B matrix
C = TH1.C;                           % C matrix
D = TH1.D;                           % D matrix

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
% Sort into ascending order

[fd,I]=sort(fd);
zeta=zeta(I);
shapes=shapes(:,I);
s=s(I); 

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

% Sort again
%--------------------------------------------------------------------------

[fd1,ii]=sort(fd1); 
zeta1=zeta1(ii); 
shapes=shapes(:,ii); 
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



Result.Parameters.NaFreq=fd1;
Result.Parameters.DampRatio=zeta1;
Result.Parameters.ModeShape=shapes;

Result.Matrices.A=A;
Result.Matrices.B=B;
Result.Matrices.C=C;
Result.Matrices.D=D;

end