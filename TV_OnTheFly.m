% (October 25, 2016)
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
% [Xhat,Rhat,Zhat] = TV_OnTheFly(Y, Lambda, Q, param)
%
% Computes an on-the-fly approximate solution 'X_otf' to the multivariate
% total variation (TV) problem applied to the multivariate observation 'Y'.
%
% Input:
%   - 'Y' (MxN vector) Multivariate observation
%   - 'Lambda' (positive scalar) TV regularization parameter
%   - 'Q' (integer) Quality parameter for the approximated algorithm
%   - 'param' Optional parameters
%       - 'param.partial' (line-vector) additionally outputs solution on partial observations of 'Y' whose length are defined by 'param.partial'
%       - 'param.disp' (boolean) if true, display the solution as the data stream incomes (default: false)
%       - 'param.X' (MxN vector) if mentionned, additionally plots the true piecewise constant signal to recover
%       - 'param.online' (boolean) if true, additionally plots the online solution (default: false)
%       - 'param.fps' number of frames per second (default: 2^8)
%
% Output:
%   - 'Xhat' Solutions
%       - 'Xhat.otf' On-the-fly solution
%       - 'Xhat.online' Online solution
%       - 'Xhat.partial' Partial solutions (if param.partial is mentionned)
%   - 'Rhat' List of instants
%       - 'Rhat.otf' List of discontinuities of 'Xhat.otf'
%       - 'Rhat.detect' List of instant where the algorithm stops and goes backward to introduce a discontinuity
%       - 'Rhat.online' List of discontinuities encountered while computing 'Xhat.otf' online
%   - 'Z' Estimation of zeta (see article)
%
% Dependency:
%   - TV_OnTheFly_aux
%       - changePointDetection
%       - changePointDetection_bound
%       - zetaEstimation            
%
function [Xhat,Rhat,Zhat] = TV_OnTheFly(Y, Lambda, Q, param)



if nargin <= 3
    param = struct;
end

if ~isfield(param,'X');
    param.X = NaN*ones(size(Y));
end

if ~isfield(param,'fps');
    param.fps = 2^8;
end

if ~isfield(param,'online');
    param.online = false;
end

if ~isfield(param,'disp')
    param.disp = false;
end

if ~isfield(param,'partial')
    K_partial = [];
    X_partial = [];
else
    K_partial = param.partial;
    X_partial = cell(1,sum((K_partial<length(Y)).*(K_partial > 0)));
end
    




detectRupt = @(z) find(abs(diff(z(1,:)))>0);


M           = size(Y,1);
N           = length(Y);
X_otf       = [];
Zhat           = [];
k0          = 1;
R_online    = [];
R_detect    = [];
if M == 1
    Init    = 1;
else
    Init    = 0;
end


% - Approximated Solution
[temp, aux, k0,kOut, Ztmp]  = TV_OnTheFly_aux(Y(:,k0), Lambda, Q, k0, k0, Init);
X_otf                       = horzcat(X_otf,temp);
Zhat                           = horzcat(Zhat,Ztmp);
kp                          = kOut;
aux.sigma2                  = ones(Q,1)*std(Y,0,2)';







%---------------------< Initialization Display >---------------------%
if param.disp
    X           = param.X;
    fps         = param.fps;
    
    
    % Plots options
    grey        = .6.*[1 1 1];
    red         = [1 0 0];
    black       = [0 0 0];
    orange      = [0.929 0.694 0.125];
    blue        = [0 0 1];
    lightblue   = [0 0.447 0.741];
    markers = {grey,black,red,blue,orange,lightblue};
    options = {{},{'LineWidth', 3},{'LineWidth', 3,'LineStyle','-.'},{'LineWidth', 2,'LineStyle','-'},{},{}};
    
    
    % Initialize figure
    Fig=figure;
    clf
    for isp = 1 : M
        ax(isp) = subplot (M, 1, isp);
        hold on;
        for ipl = 1 : length(markers)
            hp(isp, ipl) = plot (NaN, NaN, 'color',markers{ipl}, options{ipl}{:});
        end
        axis(ax(isp),[0 N min(Y(isp,:)) max(Y(isp,:))]);
        hold off;
        grid on;
        
        set(gca,'fontsize',20);
    end
    refreshdata(Fig,'caller');
    drawnow;
    xx = 1:N;
end
%---------------------</ Initialization Display >---------------------%



timer = tic;
for(k=2:N)
        
    while(kp <= k);
 
        %---------------------< Online Solution >---------------------%
        k_bnd       = kp;
        k0_bnd      = k0;
        aux_bnd     = aux;
        X_otf_bnd   = X_otf;
        while(k_bnd <= kp)
            
            [temp_bnd, aux_bnd, k0_bnd, k_bnd] = TV_OnTheFly_aux(Y(:,k_bnd), Lambda, Q, k_bnd, k0_bnd, Init, aux_bnd, kp);
            X_otf_bnd = horzcat(X_otf_bnd, temp_bnd);
            
            R_otf_bnd_tmp = detectRupt(X_otf_bnd);
            R_online     = union(R_online,R_otf_bnd_tmp);
        end
        
        
        if(kp == k)
            X_online(:,k) = X_otf_bnd(:,end);
        end
        
        
        if ~isempty(find(K_partial==kp-1, 1)) 
            X_partial{K_partial==kp-1} = X_otf_bnd;
        end

        %---------------------</ Online Solution >---------------------%
        
        
        
        
        %-------------------< On-The-Fly Solution >--------------------%
        if(kp < N)
            [temp, aux, k0,kOut, Ztmp]    = TV_OnTheFly_aux(Y(:,kp), Lambda, Q, kp, k0, Init, aux);
            R_detect    = horzcat(R_detect,aux.kdetect);
        else
            [temp, aux, k0,kOut, Ztmp]    = TV_OnTheFly_aux(Y(:,kp), Lambda, Q, kp, k0, Init, aux, N);
        end
        X_otf       = horzcat(X_otf, temp);
        Zhat           = horzcat(Zhat, Ztmp);
        %-------------------</ On-The-Fly Solution >--------------------%
        
        
        
        
        
        %-------------------------< Display >--------------------------%
        if param.disp
            for isp = 1 : M
                
                hp(isp, 1).XDataSource = 'xx(1:k)';
                hp(isp, 1).YDataSource = 'Y(isp,1:k)';
                
                hp(isp, 2).XDataSource = 'xx(1:k)';
                hp(isp, 2).YDataSource = 'X(isp,1:k)';
                
                hp(isp, 3).XDataSource = 'xx(1:length(X_otf_bnd))';
                hp(isp, 3).YDataSource = 'X_otf_bnd(isp, :)';
                
                if param.online
                hp(isp, 4).XDataSource = 'xx(1:k)';
                hp(isp, 4).YDataSource = 'X_online(isp, 1:k)';
                end
                
                hp(isp, 5).XDataSource = '[kp kp]';
                hp(isp, 5).YDataSource = 'ylim(ax(isp))';
                
                hp(isp, 6).XDataSource = '[k0 k0]';
                hp(isp, 6).YDataSource = 'ylim(ax(isp))';
                
            end
            
            elap = toc (timer);
            if elap >= 1 / fps
                for isp=1:M
                    for ilm=1:5
                        refreshdata(hp(isp,ilm),'caller');
                    end
                end
                if k>=3
                drawnow limitrate
                end
                timer = tic;
            end
        end
        %-------------------------</ Display >--------------------------%
        
        
        
       
        kp          = kOut;
    end
    
end


X_otf(:,end+1)  = X_otf(:,end);
R_otf           = detectRupt(X_otf);


Xhat.otf        = X_otf;
Xhat.online     = X_online;
Xhat.partial    = X_partial;
Rhat.otf        = R_otf;
Rhat.online     = R_online;
Rhat.detect     = R_detect;