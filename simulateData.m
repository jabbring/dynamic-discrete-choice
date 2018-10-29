%% simulateData.m
% Empirical Industrial Organization II

%{
Suppose that we have computed $\Delta U(x,a)\equiv U_1(x,a)-U_0(x,a)$ for all $(x,a)\in{\cal X}\times\{0,1\}$, using e.g. |flowpayoffs| and |fixedPoint|, and that we have specified the Markov transition matrix $\Pi$. Then, the function |simulateData| can be used to simulate $N$ independent histories $\{(X_1,A_1),\ldots,(X_T,A_T)\}$, taking $A_0=0$ as given and drawing $X_1$ from the stationary distribution of $\{X_t\}$.
%}
function [choices,iX] = simulateData(deltaU,capPi,nPeriods,nFirms)
%{
The function |simulateData| requires the following input arguments:
\begin{dictionary}
\item{|deltaU|} a $K\times 2$ matrix of which the $(i,j)$th entry is $\Delta U(x^i,j-1)$;
\item{|capPi|} the $K\times K$ Markov transition matrix $\Pi$ for $\{X_t\}$, with typical element $\Pi_{ij}=\Pr(X_{t+1}=x^j|X_t=x^i)$;
\item{|nPeriods|} the scalar number $T$ of time periods to simulate data for; and
\item{|nFirms|} the scalar number $N$ of firms to simulate data for.
\end{dictionary}
It returns
\begin{dictionary}
\item{|choices|} a $T\times N$ matrix with simulated choices, with each each column containing an independent simulation of $(A_1,\ldots,A_T)$; and
\item{|iX|} a $T\times N$ matrix with simulated states, with each each column containing the \emph{indices} (in ${\cal X}$) of an independent simulation of $(X_1,\ldots,X_T)$. For example, if $X_1=x^3$ for some firm, then 3, not $x^3$, is stored in |iX|. 
\end{dictionary}	
The function |simulateData| first stores the number $K$ of elements of |supportX| in a scalar |nSuppX|.
%}
nSuppX = size(capPi,1);
%{
	Next, it assumes that $\{X_t\}$ is ergodic and initializes the simulation by drawing $X_1$, $N$ independent times, from the stationary distribution of $\{X_t\}$. To this end, it first solves
	\[\left[\begin{array}{cccc}
					 1-\Pi_{11}  &-\Pi_{21}  &\cdots            &-\Pi_{K1}\\
					-\Pi_{12}     &1-\Pi_{22}&\cdots            &-\Pi_{K2}\\
					 \vdots      &          &\ddots            &\vdots\\
					 -\Pi_{1(K-1)}&\cdots    &~~1-\Pi_{(K-1)(K-1)}~~&-\Pi_{K(K-1)}\\
					 1           &\cdots    &1                 &1
					 \end{array}\right]P^\infty=
			\left(\begin{array}{c}
				  0\\0\\ \vdots\\ 0\\1
			\end{array}\right)
	  \]
	for the $K\times 1$ vector $P^\infty$ with the stationary probabilities $P^\infty_k=\lim_{t\rightarrow\infty}\Pr(X_t=x^k)$ and stores the result in |pInf|.
%}
oneMinPi = eye(nSuppX)-capPi';
pInf     = [oneMinPi(1:nSuppX-1,:);ones(1,nSuppX)]\[zeros(nSuppX-1,1);1];
%{
Then, it uses the auxiliary function |randomDiscrete| (see Appendix \ref{misc}) and the values stored in |pInf| to simulate a $1\times N$ vector of values of $X_1$ from the stationary distribution $P^\infty$ and stores their indices in |iX|.
%}
iX = randomDiscrete(pInf*ones(1,nFirms));
%{
Using these $N$ simulated values of $X_1$, and $N$ simulated values of $-\Delta\varepsilon_1\equiv\varepsilon_1(0)-\varepsilon_1(1)$ that are stored in |deltaEpsilon|, it simulates $N$ values of the first choice by using that $A_1=1$ if $\Delta U(X_1,0)>-\Delta\varepsilon_1$ and $A_1=0$ otherwise. These are stored in the $1\times N$ vector |choices|.
%}		
deltaEpsilon = random('ev',zeros(1,nFirms),ones(1,nFirms))-random('ev',zeros(1,nFirms),ones(1,nFirms));
choices  = deltaU(iX,1)' > deltaEpsilon;
%{
	Finally, $N$ values of $X_t$ are simulated, using the transition matrix $\Pi$ and |randomDiscrete|, and their indices added as a row to the bottom of the $(t-1)\times N$ matrix |iX|; and $N$ values of $A_t$ are simulated, using that  $A_t=1$ if $\Delta U(X_t,A_{t-1})>-\Delta\varepsilon_t$ and $A_t=0$ otherwise, and stored as a row at the bottom of the $(t-1)\times N$ matrix |choices|; recursively for $t=2,\ldots,T$.
%}
for t = 2:nPeriods
	iX = [iX;randomDiscrete(capPi(iX(end,:),:)')];
	deltaEpsilon  = random('ev',zeros(1,nFirms),ones(1,nFirms))-random('ev',zeros(1,nFirms),ones(1,nFirms));
	choices = [choices;(deltaU(iX(end,:)+nSuppX*choices(end,:)) > deltaEpsilon)];
end