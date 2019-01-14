% (April 19, 2016)
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
% [X, aux, k0, kOut] = TV_OnTheFly_aux(Y, Lambda, Q, k, k0, Init, aux, k_bnd)
%
% Auxiliary function used in TV_OnTheFly.
%
% Input:
%   - 'Y' Multivariate observation
%   - 'Lambda' TV regularization parameter
%   - 'Q' Quality parameter for the approximated algorithm
%   - 'k' Current location
%   - 'k0' Starting location of the current constant signal to assigned
%   - 'Init' [additional]
%       (Set the type of initialization (see article))
%       / Init = 0 : Forget the past (default)
%       / Init > 0 : Take into account the past
%   - 'aux' Auxiliary variables used to detect change-point locations
%   - 'k_bnd' [additional] If k==k_bnd, it enforces boundary conditions
%
% Output:
%   - 'X_otf' Solution of the proposed algorithm
%   - 'aux' Auxiliary variables used to detect change-point locations
%   - 'k0' Starting location of the current constant signal to assigned
%   - 'kOut' Next location
%
% Dependency:
%   - changePointDetection
%   - changePointDetection_bound
%   - zetaEstimation

function [X, aux, k0, kOut, Z] = TV_OnTheFly_aux(Y, Lambda, Q, k, k0, Init, aux, k_bnd)



if(nargin<=7)
    bound = 0;
else
    bound = 1;
end




% - Si c'est la premi?re observation
if(nargin<=6)
    aux.M       = size(Y,1);
    if aux.M == 1
        aux.Q   = 1;
    else
        aux.Q   = Q;
    end
    aux.z       = zetaGeneration(Lambda, aux.M, aux.Q );
    aux.c0      = zeros(1,aux.M);
    aux.z0      = zeros(1,aux.M);
    aux.lambda  = Lambda;
    aux.xmin    = NaN*ones(aux.Q,aux.M);
    aux.xmax    = NaN*ones(aux.Q,aux.M);
    aux.umin    = NaN*ones(aux.Q,aux.M);
    aux.umax    = NaN*ones(aux.Q,aux.M);
    
end



if(k==k0) % - Initialisation des variables

    aux.k0      = k0;
    aux.krupt   = k0.*ones(aux.Q,2^aux.M);
    aux.k       = NaN.*ones(1,aux.Q);
    aux.x       = NaN.*ones(aux.Q,aux.M);
    
    if(Init==0)
        aux.xmin   = ones(aux.Q,1)*Y' - aux.z ;
        aux.xmax   = ones(aux.Q,1)*Y' + aux.z ;
    else
        aux.xmin   = ones(aux.Q,1)*Y' - aux.z + ones(aux.Q,1)*(-aux.c0.*aux.z0);
        aux.xmax   = ones(aux.Q,1)*Y' + aux.z + ones(aux.Q,1)*(-aux.c0.*aux.z0);
    end
    
    aux.umin   = + aux.z;
    aux.umax   = - aux.z;
    

        
elseif(k>k0) % - Mise ? jour des variables

    
    Qlist = find(isnan(aux.x(:,1)))';
    aux.umin(Qlist,:) = aux.umin(Qlist,:) + ones(length(Qlist),1)*Y' - aux.xmin(Qlist,:);
    aux.umax(Qlist,:) = aux.umax(Qlist,:) + ones(length(Qlist),1)*Y' - aux.xmax(Qlist,:);
    
    
    if (bound && k==k_bnd)
        aux = changePointDetection_bound(aux,k,k0);  % /!\ A optimiser -> parallelisation /!\
    else
        aux = changePointDetection(aux,k,k0);        % /!\ A optimiser -> parallelisation /!\
    end
        
end



%sum(isnan(aux.x(:,1))) >= 1

if( ~isempty(find(isnan(aux.x(:,1)), 1))) % - Si une discontinuit? n'a pas ?t? observ?e pour chaque q=1...aux.nQ
    
    Z           = [];
    X           = [];
    kOut        = k+1;
    aux.kdetect = [];
    
else % - Sinon, on assigne le segment d?tect?
    
    qstar       = zetaEstimation( aux );
    aux.kdetect = k;
    aux.c0      = aux.c(qstar,:);
    aux.z0      = aux.z(qstar,:);
    
    N           = aux.k(qstar) - k0;
    k0          = aux.k(qstar);     % (Car il faut ensuite initialiser en k0)
    kOut        = k0;
    X           = [aux.x(qstar,:)]'*ones(1,N);
    Z           = [aux.z(qstar,:)]'*ones(1,N);

end


    
    
    
