%% flowpayoffs.m
%  This function computes the flow payoffs for the basic firm entry and exit model used as an example in CentER's Empirical Industrial Organization II

%{		
The function |flowpayoffs| computes the mean (over $\varepsilon$) flow payoffs, $u_0(x,a)$ and $u_1(x,a)$, for each profit and past choice pair $(x,a)\in{\cal X}\times\{0,1\}$. 
%}
function [u0,u1] = flowpayoffs(supportX,beta,delta)
%{	
It requires the following input arguments:
	\begin{dictionary}
	\item{|supportX|} a $K\times 1$ vector with the support points of the profit state $X_t$ (the elements of $\cal{X}$, consistently ordered with the Markov transition matrix $\Pi$);
	\item{|beta|} a $2\times 1$ vector that contains the intercept ($\beta_0$) and profit state slope ($\beta_1$) of the net payoffs to choice $1$; and
	\item{|delta|} a $2\times 1$ vector that contains the firm's exit ($\delta_0$) and entry ($\delta_1$) costs.
	\end{dictionary}
	It returns
	\begin{dictionary}
	\item{|u0|} a $K\times 2$ matrix of which the $(i,j)$th entry is $u_0(x^i,j-1)$ and 
	\item{|u1|} a $K\times 2$ matrix of which the $(i,j)$th entry is $u_1(x^i,j-1)$ .
	\end{dictionary}	
That is, the rows correspond to the support points of $X_t$, and the columns to the choice in the previous period, $A_{t-1}$. 

The function |flowpayoffs| first stores the number $K$ of elements of |supportX| in a scalar |nSuppX|.
%}
nSuppX = size(supportX,1);
%{
	Then, it constructs a $K\times 2$ matrix |u0| with the value of
	\[
	  \left[\begin{array}{ccc}
		u_0(x^1,0)&~~~&u_0(x^1,1)\\
		\cdot&~~~&\cdot\\ \cdot&~~~&\cdot\\ \cdot&~~~&\cdot\\
		u_0(x^{K},0)&~~~&u_0(x^{K},1)
		\end{array}\right]
	  =
	  \left[\begin{array}{ccc}
		0&~~~&-\delta_0\\
		\cdot&~~~&\cdot\\ \cdot&~~~&\cdot\\ \cdot&~~~&\cdot\\
		0&~~~&-\delta_0
		\end{array}\right]
	\]
	and a $K\times 2$ matrix |u1| with the value of
	\[
	  \left[\begin{array}{ccc}
		u_1(x^1,0)&~~~&u_1(x^1,1)\\
		\cdot&~~~&\cdot\\ \cdot&~~~&\cdot\\ \cdot&~~~&\cdot\\
		u_1(x^K,0)&~~~&u_1(x^{K},1)
		\end{array}\right]
	  =
	  \left[\begin{array}{ccc}
		\beta_0+\beta_1 x^1-\delta_1&~~~&\beta_0+\beta_1 x^1\\
		\cdot&~~~&\cdot\\ \cdot&~~~&\cdot\\ \cdot&~~~&\cdot\\
		\beta_0+\beta_1 x^{K}-\delta_1&~~~&\beta_0+\beta_1 x^{K}
		\end{array}\right].
	\]
%}
u0 = [zeros(nSuppX,1) -delta(1)*ones(nSuppX,1)];
u1 = [ones(nSuppX,1) supportX]*beta*[1 1]-delta(2)*ones(nSuppX,1)*[1 0];
%{
	You can change the specification of the flow profits by editing these two lines of code.
%}