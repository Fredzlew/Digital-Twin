clear all;clc;

% using the analyze FDD tpo find the frequency and eigenvectors
[Frq,Us] = FDD('data_1_2_1.xlsx',100)

% normalizing FDD mode shapes
MVec_x = max(Us); % start normalization
mVec_x = min(Us);
for j = 1:5
    if abs(MVec_x(j)) > abs(mVec_x(j))
        mxVec_x(j) = MVec_x(j);
    else
        mxVec_x(j) = mVec_x(j);
    end
    for l = 1:5
        phi_FDD(l,j) = Us(l,j)/mxVec_x(j);
    end
end % end normalization

% taking the data from the numerical model
filename = load('modelprop.mat'); % Loads mass, stiffness matrices, height and mode shape vector

% plotting the mode shapes from numerical and FDD to compare
nm = 5; % number of modes
x = [0, filename.H];
phi = [zeros(1,length(filename.U)); filename.U];
fig = figure;
fig.Position=[100 100 1600 700];
for i=1:nm
    subplot(1,nm,i)
    hold on
    if i == 5
        plot(phi(:,i),x,'-m')
        plot(-[0;phi_FDD(:,i)],x,'-g')
        plot(phi(2:end,i),x(2:end),'b.','markersize',30)
        plot(-phi_FDD(1:end,i),x(2:end),'b.','markersize',30)
        title(['f = ' num2str(Frq(i)) ' Hz'],sprintf('Mode shape %d',i),'FontSize',14)
        xline(0.0,'--')
        xlim([-1.1,1.1])
        ylim([0,x(end)])   
    else
        plot(phi(:,i),x,'-m')
        plot([0;phi_FDD(:,i)],x,'-g')
        plot(phi(2:end,i),x(2:end),'b.','markersize',30)
        plot(phi_FDD(1:end,i),x(2:end),'b.','markersize',30)
        title(['f = ' num2str(Frq(i)) ' Hz'],sprintf('Mode shape %d',i),'FontSize',14)
        xline(0.0,'--')
        xlim([-1.1,1.1])
        ylim([0,x(end)])
        if i==1
            legend('Numerical','Approximation','Location','northwest')
        else
        end
    end
end
sgtitle(sprintf('Numerical vs FDD Online')) 
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Height [m]','FontSize',14);
xlabel(han,'Deflection [-]','FontSize',14);