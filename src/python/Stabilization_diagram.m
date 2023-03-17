% parameters
clc; clear; close all;
addpath(genpath('data'),genpath('functions'),genpath('OMA'),genpath('python'))
% Loading modal parameters from OMA 
SSIstab2 = readNPY('stab.npy');
% definition:
% freq = SSIstab(:,1);
% modelorder = SSIstab(:,2);
% label = SSIstab(:,3);
% damp = SSIstab(:,4);

% Only use frequency up to 20 
k=1;

for i = 1:length(SSIstab2)
 if SSIstab2(i,1) <= 20
     SSIstab(k,:) = SSIstab2(i,:);
     k=k+1;
 end
end

% sortering out from label
m = 1;
n = 1;
s = 1;
q = 1;
b = 1;
for i = 1:length(SSIstab)
    if SSIstab(i,3) == 0 % label = 0
        freq0(m) = SSIstab(i,1);
        modelorder0(m) = SSIstab(i,2);
        damp0(m) = SSIstab(i,4);
        m = m+1;
    elseif SSIstab(i,3) == 1 % label = 1
        freq1(n) = SSIstab(i,1);
        modelorder1(n) = SSIstab(i,2);
        damp1(n) = SSIstab(i,4);
        n = n+1;
    elseif SSIstab(i,3) == 2 % label = 2
        freq2(s) = SSIstab(i,1);
        modelorder2(s) = SSIstab(i,2);
        damp2(s) = SSIstab(i,4);
        s = s+1;
    elseif SSIstab(i,3) == 3 % label = 3
        freq3(q) = SSIstab(i,1);
        modelorder3(q) = SSIstab(i,2);
        damp3(q) = SSIstab(i,4);
        q = q+1;    
    elseif SSIstab(i,3) == 4 % label = 4
        freq4(b) = SSIstab(i,1);
        modelorder4(b) = SSIstab(i,2);
        damp4(b) = SSIstab(i,4);
        b = b+1;
    end
end

% plotting
plot3(damp0,freq0,modelorder0,'ro',damp1,freq1,modelorder1,'mo',damp2,freq2,modelorder2,'yo',damp3,freq3,modelorder3,'bo',damp4,freq4,modelorder4,'go')
legend('Unstable pole','Stable for frequency','Stable for frequency and damping','Stable for frequency and mode shape','Stable pole','FontSize', 12)
xlabel('Damping ratio [-]','FontSize', 14)
ylabel('Frequency [Hz]','FontSize', 14)
zlabel('Model order','FontSize', 14)
xlim([0 0.05])