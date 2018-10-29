%% estimatePi.m
% This function computes estimates of the Markov transition matrix of the observed state variables in the basic firm entry and exit model used as an example in CentER's Empirical Industrial Organization II. 

%{
The function |estimatePi| computes a frequency estimate of the Markov transition matrix $\Pi$ from state transition data.
%}
function piHat = estimatePi(iX,nSuppX)
%{
	It requires the following input arguments:
	\begin{dictionary}
	\item{|iX|} a $T\times N$ matrix with indices of observed states $x_{tn}$ in ${\cal X}$ (for example, if $x_{11}=x^3$, then the first element of |iX| is 3, not $x^3$); 
	\item{|nSuppX|} the scalar number $K$ of support points of the profit state $X_t$ (the number of elements of $\cal{X}$)
	\end{dictionary}
	and returns
	\begin{dictionary}
	\item{|piHat|} an estimate $\hat\Pi$ of the $K\times K$ Markov transition matrix $\Pi$ for $\{X_t\}$, with typical element $\hat\Pi_{ij}$ equal to the sample frequency of transitions to state $j$ among all transitions from state $i$.
	\end{dictionary}	
The function |estimatePi| first stores the number of time periods $T$ in a scalar |nPeriods|.
%}
nPeriods=size(iX,1);
%{
Then, for each pair $(i,j)\in\{1,\ldots,K\}\times\{1,\ldots,K\}$, it estimates the probability $\Pi_{ij}=\Pr(X_{t+1}=x^j|X_t=x^i)$ by the appropriate sample frequency, the number of transitions from $i$ to $j$ divided by the total number of transitions from $i$ in the data |iX|.
%}
for i=1:nSuppX
    for j=1:nSuppX
        piHat(i,j)=sum(sum((iX(2:nPeriods,:)==j)&(iX(1:nPeriods-1,:)==i)))/sum(sum((iX(1:nPeriods-1,:)==i)));
    end
end
%{ 
Note that |estimatePi| requires a positive number of transition observations from each state. More generally, the frequency estimator that it implements only performs well with samples that are large relative to the state space. With relatively small samples, the frequency estimator should be replaced by one that smoothes across support points.
%}