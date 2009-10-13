function plotMatlabOutputs(fileName,Type)
% plot the MATLAB outputs
%
% Written by T.Y. Yang on 09/12/2009

% load the Matlab outputs
load(fileName);

% select the plot cases
switch Type
    %%%%%%%%%%%
    case 'DMNR'
        % % plot the element hysteresis
        % figure;
        % for j=1:numElem
        %     subplot(numElem,1,j);
        %     plot(u(j,:),pr(j,:),'b-');
        %     xlabel(['u' num2str(j)]);
        %     ylabel(['pr' num2str(j)]);
        %     grid
        % end
        
        % plot the element hysteresis
        figure;
        for j=1:numElem
            subplot(numElem,1,j);
            plot(uall(j,1:count-1),prall(j,:),'b-');
            xlabel(['uall' num2str(j)]);
            ylabel(['prall' num2str(j)]);
            grid
        end
        
        % % plot the element force history
        % figure;
        % for j=1:numElem
        %     subplot(numElem,1,j);
        %     plot(t,pr(j,:),'b-');
        %     xlabel('Time [sec]')
        %     ylabel(['pr' num2str(j)]);
        %     grid
        % end
        %
        % % plot the element displacement history
        % figure;
        % for j=1:numElem
        %     subplot(numElem,1,j);
        %     plot(t,u(j,:),'b-');
        %     xlabel('Time [sec]')
        %     ylabel(['u' num2str(j)]);
        %     grid
        % end
        %
        % % plot the Node response history
        % figure;
        % for j=1:ndf
        %     subplot(3,ndf,j);
        %     plot(t,U(j,:),'b-');
        %     xlabel('Time [sec]')
        %     ylabel(['U' num2str(j)]);
        %     grid
        %     subplot(3,ndf,j+ndf);
        %     plot(t,Udot(j,:),'b-');
        %     xlabel('Time [sec]')
        %     ylabel(['Udot' num2str(j)]);
        %     grid
        %     subplot(3,ndf,j+2*ndf);
        %     plot(t,Udot(j,:),'b-');
        %     xlabel('Time [sec]')
        %     ylabel(['Udot' num2str(j)]);
        %     grid
        % end
        
        figure;
        for j=1:ndf
            subplot(ndf,1,j);
            plot(R(j,:));
            xlabel(['[-]' num2str(j)]);
            ylabel(['R' num2str(j)]);
            grid
        end
        
        % % error
        % figure;
        % plot(t,errorNorms)
        % ylabel('Error')
        % xlabel('Time [sec]')
        % grid
        
        % Iterations
        figure;
        plot(t,iters)
        ylabel('Error')
        xlabel('Time [sec]')
        grid
        
        % experimental element trial displacement
        figure;
        eD1 = load('ElementDisp1.txt');
        plot(eD1,'b-')
        ylabel('trialDisp')
        xlabel('Step [-]')
        grid
        
        % plot all experimental force
        figure;
        plot(prall(1,2:end),'b-')
        ylabel('prall')
        xlabel('Step [-]')
        grid
        
        
        % % ground motion
        % figure;
        % plot(t,ag(1:length(t)),'b-')
        % ylabel('Ag')
        % xlabel('Time [sec]')
        % grid
        
    %%%%%%%%%%%%%%%%%%%%%%
    case 'DMNRLimit'
        % plot the element hysteresis
        figure;
        for j=1:numElem
            subplot(numElem,1,j);
            plot(u(j,:),pr(j,:),'b-');
            xlabel(['u' num2str(j)]);
            ylabel(['pr' num2str(j)]);
            grid
        end
        
        % plot the element force history
        figure;
        for j=1:numElem
            subplot(numElem,1,j);
            plot(t,pr(j,:),'b-');
            xlabel('Time [sec]')
            ylabel(['pr' num2str(j)]);
            grid
        end
        
        % plot the element displacement history
        figure;
        for j=1:numElem
            subplot(numElem,1,j);
            plot(t,u(j,:),'b-');
            xlabel('Time [sec]')
            ylabel(['u' num2str(j)]);
            grid
        end
        
        % plot the Node response history
        figure;
        for j=1:ndf
            subplot(3,ndf,j);
            plot(t,U(j,:),'b-');
            xlabel('Time [sec]')
            ylabel(['U' num2str(j)]);
            grid
            subplot(3,ndf,j+ndf);
            plot(t,Udot(j,:),'b-');
            xlabel('Time [sec]')
            ylabel(['Udot' num2str(j)]);
            grid
            subplot(3,ndf,j+2*ndf);
            plot(t,Udot(j,:),'b-');
            xlabel('Time [sec]')
            ylabel(['Udot' num2str(j)]);
            grid
        end
        
        % error
        figure;
        plot(t,errorNorms,'b-')
        ylabel('Error')
        xlabel('Time [sec]')
        grid
        
        %         % experimental element trial displacement
        %         figure;
        %         eD1 = load('ElementDisp1.txt');
        %         plot(eD1,'b-')
        %         ylabel('trialDisp')
        %         xlabel('Step [-]')
        %         grid
        
        % % ground motion
        % figure;
        % plot(t,ag(1:length(t)),'b-')
        % ylabel('Ag')
        % xlabel('Time [sec]')
        % grid
    %%%%%%%%%%%%%%%%%%%%%%
    case 'FMNR'
        % % plot the element hysteresis
        figure;
        for j=1:numElem
            subplot(numElem,1,j);
            plot(u(j,:),pr(j,:),'r-');
            xlabel(['u' num2str(j)]);
            ylabel(['pr' num2str(j)]);
            grid
        end
        
        % plot the element hysteresis
        figure;
        for j=1:numElem
            subplot(numElem,1,j);
            plot(uall(j,:),prall(j,1:count-1),'r-');
            xlabel(['uall' num2str(j)]);
            ylabel(['prall' num2str(j)]);
            grid
        end
        
        % plot the element force history
        figure;
        for j=1:numElem
            subplot(numElem,1,j);
            plot(t,pr(j,:),'r-');
            xlabel('Time [sec]')
            ylabel(['pr' num2str(j)]);
            grid
        end
        
        % plot the element displacement history
        figure;
        for j=1:numElem
            subplot(numElem,1,j);
            plot(t,u(j,:),'r-');
            xlabel('Time [sec]')
            ylabel(['u' num2str(j)]);
            grid
        end
        
        % % plot the Node response history
        % figure;
        % for j=1:ndf
        %     subplot(3,ndf,j);
        %     plot(t,U(j,:),'r-');
        %     xlabel('Time [sec]')
        %     ylabel(['U' num2str(j)]);
        %     grid
        %     subplot(3,ndf,j+ndf);
        %     plot(t,Udot(j,:),'r-');
        %     xlabel('Time [sec]')
        %     ylabel(['Udot' num2str(j)]);
        %     grid
        %     subplot(3,ndf,j+2*ndf);
        %     plot(t,Udot(j,:),'r-');
        %     xlabel('Time [sec]')
        %     ylabel(['Udot' num2str(j)]);
        %     grid
        % end
        
        figure;
        for j=1:ndf
            subplot(ndf,1,j);
            plot(R(j,:),'r-');
            xlabel(['[-]' num2str(j)]);
            ylabel(['R' num2str(j)]);
            grid
        end
        
        % error
        figure;
        plot(t,errorNorms,'r-')
        ylabel('Error')
        xlabel('Time [sec]')
        grid
        
        % Iteration
        figure;
        plot(t,iters,'r-')
        ylabel('Error')
        xlabel('Time [sec]')
        grid
        
        % plot all experimental displace
        figure;
        plot(uall(1,2:end),'r-')
        ylabel('uall')
        xlabel('Step [-]')
        grid
        
        % experimental element trial force
        figure;
        eF1 = load('ElementForce1.txt');
        plot(eF1,'r-')
        ylabel('trialForce')
        xlabel('Step [-]')
        grid
        
        % % ground motion
        % figure;
        % plot(t,ag(1:length(t)),'r-')
        % ylabel('Ag')
        % xlabel('Time [sec]')
        % grid
    %%%%%%%%%%%%%%%%
    case 'FMNRLimit'
        % plot the element hysteresis
        figure;
        for j=1:numElem
            subplot(numElem,1,j);
            plot(u(j,:),pr(j,:),'r-');
            xlabel(['u' num2str(j)]);
            ylabel(['pr' num2str(j)]);
            grid
        end
        
        % plot the element force history
        figure;
        for j=1:numElem
            subplot(numElem,1,j);
            plot(t,pr(j,:),'r-');
            xlabel('Time [sec]')
            ylabel(['pr' num2str(j)]);
            grid
        end
        
        % plot the element displacement history
        figure;
        for j=1:numElem
            subplot(numElem,1,j);
            plot(t,u(j,:),'r-');
            xlabel('Time [sec]')
            ylabel(['u' num2str(j)]);
            grid
        end
        
        % plot the Node response history
        figure;
        for j=1:ndf
            subplot(3,ndf,j);
            plot(t,U(j,:),'r-');
            xlabel('Time [sec]')
            ylabel(['U' num2str(j)]);
            grid
            subplot(3,ndf,j+ndf);
            plot(t,Udot(j,:),'r-');
            xlabel('Time [sec]')
            ylabel(['Udot' num2str(j)]);
            grid
            subplot(3,ndf,j+2*ndf);
            plot(t,Udot(j,:),'r-');
            xlabel('Time [sec]')
            ylabel(['Udot' num2str(j)]);
            grid
        end
        
        % error
        figure;
        plot(t,errorNorms,'r-')
        ylabel('Error')
        xlabel('Time [sec]')
        grid
        
        % experimental element trial force
        figure;
        eF1 = load('ElementForce1.txt');
        plot(eF1,'r-')
        ylabel('trialForce')
        xlabel('Step [-]')
        grid
        
        % % ground motion
        % figure;
        % plot(t,ag(1:length(t)),'r-')
        % ylabel('Ag')
        % xlabel('Time [sec]')
        % grid
end
