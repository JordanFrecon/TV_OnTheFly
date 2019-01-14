% (November 09, 2016)
% 
% Author:
% Jordan Frecon (jordan.frecon@ens-lyon.fr) 
% --> Contact author Nelly Pustelnik (nelly.pustelnik@ens-lyon.fr)
% 
% Contributors:
% Nelly Pustelnik (nelly.pustelnik@ens-lyon.fr)
% Patrice Abry (patrice.abry@ens-lyon.fr)
% Laurent Condat (laurent.condat@gipsa-lab.grenoble-inp.fr)
% 
% This software is governed by the CeCILL license under French law and
% abiding by the rules of distribution of free software.  You can  use,
% modify and/ or redistribute the software under the terms of the CeCILL
% license as circulated by CEA, CNRS and INRIA at the following URL
% "http://www.cecill.info".
% 
% As a counterpart to the access to the source code and  rights to copy,
% modify and redistribute granted by the license, users are provided only
% with a limited warranty  and the software's author,  the holder of the
% economic rights,  and the successive licensors  have only  limited
% liability.
% 
% In this respect, the user's attention is drawn to the risks associated
% with loading,  using,  modifying and/or developing or reproducing the
% software by the user in light of its specific status of free software,
% that may mean  that it is complicated to manipulate,  and  that  also
% therefore means  that it is reserved for developers  and  experienced
% professionals having in-depth computer knowledge. Users are therefore
% encouraged to load and test the software's suitability as regards their
% requirements in conditions enabling the security of their systems and/or
% data to be ensured and,  more generally, to use and operate it in the
% same conditions as regards security.
% 
% The fact that you are presently reading this means that you have had
% knowledge of the CeCILL license and that you accept its terms.
%
%--------------------------------------------------------------------------
% TV On-The-Fly: On-the-fly approximation of the multivariate TV problem  
%                                                                         
% For theoretical aspects please refer to :                               
% J. Frecon, N. Pustelnik, P. Abry, L. Condat                             
% On-The-Fly Approximation of Multivariate Total Variation Minimization   
% IEEE Transactions on Signal Processing,                                 
% Vol. 64, Issue 9, pp. 2355-2364, May. 2016                              
%--------------------------------------------------------------------------
%
% demo_TV_OnTheFly.m
%
% Provides one denoising example of a multivariate signal (with M components)
% 
% Two display settings are available by changing 'param.disp' :
% 1) if false (Offline setting), shows the solution for a given observation.
% 2) if true (Online setting), shows the solution as the observation stream incomes. 
% 
% The main function is TV_OnTheFly.m.
%
%%

mydir  = which('demo_TV_OnTheFly.m');
idcs   = strfind(mydir,'/');
newdir = mydir(1:idcs(end));       
addpath(genpath(newdir));
grey = .6.*[1 1 1];


% Signal
M = 3;
X = [ones(M,250), 3*ones(M,200), 2*ones(M,300), 0*ones(M,250)];
N = length(X);
s = randn(M,1)*ones(1,N);
Y = X + s.*randn(M,N);




%% On-The-Fly TV Algorithm


% TV Parameters
lambda          = 40;           % Regularization parameter
Q               = 20;           % The higher Q, the higher the quality and execution time


% Optional Parameters
param.partial   = [400 800];    % if mentionned, additionally outputs solution on partial observations of 'Y'
param.disp      = true;         % if true, display the solution as the data stream incomes
param.X         = X;            %   if mentionned, additionally plots 'X'
param.online    = true;         %   if true, additionally plots the online solution
param.fps       = 2^8;          %   #frames per second


% Computation
[Xhat,Rhat,Zhat]    = TV_OnTheFly(Y,lambda,Q,param);



%% DISPLAY

% Display final solution
figure(11);clf
for mm=1:M
    subplot(M,1,mm)
    plot(Y(mm,:),'color',grey); hold on;
    plot(X(mm,:),'k','linewidth',3);
    plot(Xhat.otf(mm,:),'-.r','linewidth',3);
    plot(Xhat.online(mm,:),'-b','linewidth',2);
    axis([0 N min(Y(mm,:)) max(Y(mm,:))]);
    grid on;
    
    if mm == 1
        legend('Observation','True Signal','On-the-fly Solution','Online Soluton');
        title('Final solution(s)');
    end
    set(gca,'fontsize',20);
end
drawnow;


% Display partial solutions
if isfield(param,'partial') && ~isempty(param.partial)
    figure(12);clf
    for mm=1:M
        subplot(M,1,mm)
        plot(Y(mm,:),'color',grey); hold on;
        plot(X(mm,:),'k','linewidth',3);
        for ii=1:length(Xhat.partial)
            plot(Xhat.partial{ii}(mm,:),'-.','linewidth',3);
        end
        axis([0 N min(Y(mm,:)) max(Y(mm,:))]);
        grid on;
        
        if mm == 1
            legend('Observation','True Signal');
            title('Partial solution(s)');
        end
        set(gca,'fontsize',20);
    end
    drawnow;
end

% List of discontinuities
disp('True positions of discontinuities')
disp(find(abs(diff(X(1,:)))>0))
disp('Estimated positions')
disp(Rhat.otf)
disp('When they are detected on-the-fly')
disp(Rhat.detect)
disp('Discontinuities detected while computing the online solution');
disp(Rhat.online)



