%% fixedPoint.m
% This function computes a fixed point of the contraction Psi by successive approximations

%{
	The function |fixedPoint| computes the fixed point $U$ of the Bellman-like operator $\Psi$, using the method of successive approximations. It is easy to show that $\Psi$ is a contraction, so that the method of successive approximations converges linearly and globally to a point within a positive maximum absolute distance |tolFixedPoint| (see Appendix \ref{contractions.computation}).
%}
function [capU0,capU1] = fixedPoint(u0,u1,capPi,rho,tolFixedPoint,bellman,capU0,capU1)
%{
	It requires the following input arguments:
	\begin{dictionary}
	\item{|u0|} a $K\times 2$ matrix of which the $(i,j)$th entry is $u_0(x^i,j-1)$;
	\item{|u1|} a $K\times 2$ matrix of which the $(i,j)$th entry is $u_1(x^i,j-1)$;
	\item{|capPi|} the $K\times K$ Markov transition matrix $\Pi$ for $\{X_t\}$, with typical element $\Pi_{ij}=\Pr(X_{t+1}=x^j|X_t=x^i)$; 
	\item{|rho|} a scalar with the value of the discount factor $\rho$;
	\item{|tolFixedPoint|} a scalar tolerance level that is used to determine convergence of the successive approximations;
	\item{|bellman|} the handle to the function |[capU0,capU1] = bellman(capU0,capU1,u0,u1,capPi,rho)| that iterates once on $\Psi$;
	\item{|capU0|} a $K\times 2$ matrix of which the $(i,j)$th entry is a starting value for $U_0(x^i,j-1)$ (optional; set to |[]| to select default starting value);
	\item{|capU1|} a $K\times 2$ matrix of which the $(i,j)$th entry is a
	starting value for $U_1(x^i,j-1)$ (optional; set to |[]| to select default starting value).
	\end{dictionary}
	It returns
	\begin{dictionary}
	\item{|capU0|} a $K\times 2$ matrix of which the $(i,j)$th entry is $U_0(x^i,j-1)$; and
	\item{|capU1|} a $K\times 2$ matrix of which the $(i,j)$th entry is $U_1(x^i,j-1)$.
	\end{dictionary}	
The function |fixedPoint| first stores the number $K$ of elements of |supportX| in a scalar |nSuppX|.
	%}
nSuppX = size(capPi,1);
%{
The starting values for $U_0$ and $U_1$ are set to $0$ if the input arguments |capU0| and |capU1| are empty.
%}
if isempty(capU0) 
	capU0 = zeros(nSuppX,2);
end
if isempty(capU1)
	capU1 = zeros(nSuppX,2);
end
%{
	The $K\times 2$ matrices |inU0| and |inU1| store the values of $U$ that are fed into the operator $\Psi$, for comparison with the value of $\Psi(U)$ that is subsequently stored in |capU0| and |capU1|. They are initialized to deviate from |capU0| and |capU1| by more than |tolFixedPoint|, so that the |while| statement allows at least one iteration of $\Psi$, and stops as soon as $\max\{\max|\Psi_0(U)-U_0|,\max|\Psi_1(U)-U_1|\}$ no longer exceeds the tolerance level in |tolFixedPoint|.
%}	
inU0 = capU0+2*tolFixedPoint;
inU1 = capU1+2*tolFixedPoint;
while (max(max(abs(inU0-capU0)))>tolFixedPoint) || (max(max(abs(inU1-capU1)))>tolFixedPoint);
	inU0 = capU0;
	inU1 = capU1;
	[capU0,capU1] = bellman(inU0,inU1,u0,u1,capPi,rho);
end
%{
You can replace |fixedPoint| by another function if you want to use alternative methods, such as Newton methods, for computing the fixed point $U$, or if you want to work with finite horizon problems.  
%}