%% randomDiscrete.m
%  This function simulates from a discrete distribution

%{		
	The function |randomDiscrete| returns a random draw from the distribution of $(Y_1,\ldots,Y_n)$; with $Y_1,\ldots,Y_n$ independently discretely distributed with (not necessarily identical) distributions on $\{1,2,\ldots,k\}$ for some $2\leq k<\infty$.
%}
function y = randomDiscrete(p)
%{	
It requires one input argument:
\begin{dictionary}
\item{|p|} a $k\times n$ matrix with $(k,n)$th element $\Pr(Y_n=k)$.
\end{dictionary}
It returns
\begin{dictionary}
\item{|y|} a $1\times n$ vector with a random draw $(y_1,\ldots,y_N)$ from $(Y_1,\ldots,Y_N)$.
\end{dictionary}	
The function |randomDiscrete| first stores the number $k$ of support points in a scalar |nSupp| and the number $n$ of random variables in a scalar |nVar|.
%}
nSupp = size(p,1);
nVar = size(p,2);
%{
Then, it creates a $(k-1)\times n$ matrix |uniformDraws| with $k-1$ identical rows, containing $n$ independent draws from a standard uniform distribution, and computes the $k\times n$ matrix |cumulativeP| with $(k,n)$th element $\Pr(Y_n\leq k)$.
%}		
uniformDraws = ones(nSupp-1,1)*random('unif',zeros(1,nVar),ones(1,nVar));
cumulativeP = cumsum(p);
%{
	Finally, for each $n$, it sets $y_n$ equal to 1 plus the number of elements of $\{\Pr(Y_n\leq 1), \ldots, \Pr(Y_n\leq k-1)$\} that are weakly smaller than the uniform random draw in the $n$th column of |uniformDraws|.
%}
y = sum([ones(1,nVar);cumulativeP(1:nSupp-1,:)<=uniformDraws]);
