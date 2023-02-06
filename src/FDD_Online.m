clear all;clc;

% using the analyze FDD tpo find the frequency and eigenvectors
[Frq,phi_FDD] = FDD('data_1_2_1.xlsx',100)

% taking the data from the numerical model
filename = load('modelprop.mat'); % Loads mass, stiffness matrices, height and mode shape vector

% plotting the mode shapes from numerical and FDD to compare
nm = 5; % number of modes
x = [0, H];
phi = [zeros(1,length(U)); U];
fig = figure;
fig.Position=[100 100 1600 700];
for i=1:nm
    subplot(1,nm,i)
    hold on
    plot(phi(:,i),x,'-b')
    plot(phi_FDD(:,i),x,'-g')
    plot(phi(2:end,i),x(2:end),'b.','markersize',30)
    plot(phi_FDD(2:end,i),x(2:end),'b.','markersize',30)
%     title(['f = ' num2str(fn(i)) ' Hz'],sprintf('Mode shape %d',i),'FontSize',14)
    xline(0.0,'--')
    xlim([-1.1,1.1])
    ylim([0,x(end)])
end

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Height [m]','FontSize',14);
xlabel(han,'Deflection [-]','FontSize',14);