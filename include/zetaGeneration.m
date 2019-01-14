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
% [ zeta ] = zetaGeneration(Lambda, M, Q )
%
% Generates 'Q' positive values of zeta drawn randomly on a
% hypersphere of dimension 'M' with radius sqrt('Lambda').
%
% Input:
%   - 'Lambda' TV regularization parameter
%   - 'M' Number of components of the multivariate observation
%   - 'Q' Quality parameter for the approximated algorithm
% 
% Output:
%   - 'zeta' auxiliary variable
%
%
function [ zeta ] = zetaGeneration(Lambda, M, Q )

if M == 1
    zeta = Lambda;
else
    
    psi             = rand(M,Q)*pi/2;
    zeta            = cos(psi);
    zeta(end,:)     = 1;
    zeta            = Lambda.*zeta;
    
    for k=1:M-1
        Tmp         = circshift(sin(psi),[k 0]);
        Tmp(1:k,:)  = 1;
        zeta        = zeta.*Tmp;
    end
    
    zeta = zeta';
end

end

