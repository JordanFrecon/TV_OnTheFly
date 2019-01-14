TV_OnTheFly

*****************************************************************************************************************
* author: Jordan Frecon  											*
* institution: Univ Lyon, Ens de Lyon, Univ Claude Bernard, CNRS, Laboratoire de Physique, F-69342 Lyon, France *
* date: April 19 2016     	              									*
* License CeCILL-B                                    								*
*****************************************************************************************************************


*********************************************************
* RECOMMENDATIONS:                                   	*
* This toolbox is designed to work with Matlab 2015.a   *
*********************************************************

------------------------------------------------------------------------------------------------------------------------
DESCRIPTION:
This toolbox provides an efficient on-the-fly (yet approximate) solution of the $\ell_{1,2}$-TV minimization problem. 
The corresponding method is based on the local validation of the Karush-Kuhn-Tucker conditions on the dual problem. 

This toolbox consists of 1 subfolder containing MATLAB functions designed for the proposed algorithm.

------------------------------------------------------------------------------------------------------------------------
SPECIFICATIONS for using TV_OnTheFly:

One demo file 'demo_TV_OnTheFly.m' is proposed.
It provides one denoising example of a multivariate signal (with M components)
Two display setting are proposed depending on the variable ‘param.disp’, namely:
1) Offline setting, showing the solution for a given observation.
2) Online setting, showing the solution as the observation stream incomes.

The main function is TV_OnTheFly.m.

------------------------------------------------------------------------------------------------------------------------
RELATED PUBLICATION:

# J. Frecon, N. Pustelnik, P. Abry, L. Condat
On-The-Fly Approximation of Multivariate Total Variation Minimization
IEEE Transactions on Signal Processing, Vol. 64, Issue 9, pp. 2355-2364, May. 2016

------------------------------------------------------------------------------------------------------------------------