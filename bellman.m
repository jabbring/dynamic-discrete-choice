%% bellman.m
%  This function iterates once on the Bellman-like operator for the basic firm entry and exit model used as an example in CentER's Empirical Industrial Organization II

%{		
The expected discounted profits net of $\varepsilon_{t}(a)$ immediately following choice $a$ at time $t$, $U_a(X_t,A_{t-1})$, satisfy a recursive system $U_0=\Psi_0(U_0,U_1)$ and $U_1=\Psi_1(U_0,U_1)$ or, with $U\equiv(U_0,U_1)$ and $\Psi\equiv(\Psi_0,\Psi_1)$, simply $U=\Psi(U)$ (see e.g. \cite{nh94:rust}). Here, $\Psi$ is a Bellman-like operator that embodies the model's dynamic specification and that depends on the flow payoffs $u_0$ and $u_1$, the Markov transition matrix $\Pi$, and the discount factor $\rho$. Its elements are the operators $\Psi_0$ and $\Psi_1$, with $\Psi_a$ implicitly defined by the right hand side of
	\begin{equation}\label{eq:bellman}
		U_a(X_t,A_{t-1})=u_a(X_t,A_{t-1})+\rho\mathbb{E}\left[R_a(X_{t+1})|X_t\right],
	\end{equation}
where $R_a(x)\equiv\mathbb{E}\left[\max\{U_0(x,a)+\varepsilon_{t+1}(0),U_1(x,a)+\varepsilon_{t+1}(1)\}\right]$ is \cite{McFadden_1981_probabilistic}'s social surplus for the binary choice problem with utilities $U_0(x,a)+\varepsilon_{t+1}(0)$ and $U_1(x,a)+\varepsilon_{t+1}(1)$. The first term in the right-hand side of (\ref{eq:bellman}) equals period $t$'s flow profits following choice $a$, net of $\varepsilon_t(a)$. The second term equals the expected discounted profits from continuing into period $t+1$ with choice $a$ (note that, given the choice $a$, these continuation payoffs do not depend on the past choice $A_{t-1}$). The (mean-zero) extreme value assumptions on $\varepsilon_{t+1}(0)$ and $\varepsilon_{t+1}(1)$ imply that
	\begin{equation}\label{eq:surplus}	
		R_a(x)=\log\left\{\exp\left[U_0(x,a)\right]+\exp\left[U_1(x,a)\right]\right\}.
	\end{equation}
The function |bellman| applies the operator $\Psi$ once to input values of $U$ and returns $\Psi(U)$.
%}
function [capU0,capU1] = bellman(capU0,capU1,u0,u1,capPi,rho)
%{
It requires the following input arguments:
\begin{dictionary}
\item{|capU0|} a $K\times 2$ matrix of which the $(i,j)$th entry is $U_0(x^i,j-1)$;
\item{|capU1|} a $K\times 2$ matrix of which the $(i,j)$th entry is $U_1(x^i,j-1)$;
\item{|u0|} a $K\times 2$ matrix of which the $(i,j)$th entry is $u_0(x^i,j-1)$;
\item{|u1|} a $K\times 2$ matrix of which the $(i,j)$th entry is $u_1(x^i,j-1)$;
\item{|capPi|} the $K\times K$ Markov transition matrix $\Pi$ for $\{X_t\}$, with typical element $\Pi_{ij}=\Pr(X_{t+1}=x^j|X_t=x^i)$; and
\item{|rho|} a scalar with the value of the discount factor $\rho$. 
\end{dictionary}
It returns
\begin{dictionary}
\item{|capU0|} a $K\times 2$ matrix of which the $(i,j)$th entry is $U_0(x^i,j-1)$;
\item{|capU1|} a $K\times 2$ matrix of which the $(i,j)$th entry is
$U_1(x^i,j-1)$.
\end{dictionary}	
	To this end, |bellman| first computes the surpluses $R_0(x)$ and $R_1(x)$ in (\ref{eq:surplus}) for all $x\in\cal{X}$ and stacks these in $K\times 1$ vectors |R0| and |R1|.
%}
r0 = log(exp(capU0(:,1))+exp(capU1(:,1)));
r1 = log(exp(capU0(:,2))+exp(capU1(:,2)));
%{
Then, it applies (\ref{eq:bellman}) to compute new values of |capU0| and |capU1|.
%}
capU0 = u0 + rho*capPi*r0*[1 1];
capU1 = u1 + rho*capPi*r1*[1 1];
%{
Here, the conditional expectation over $X_{t+1}$ in (\ref{eq:bellman}) is taken by premultiplying the vectors |r0| and |r1| by the Markov transition matrix |capPi|. The vectors |r0| and |r1| are postmultiplied by |[1 1]| because the surpluses, and therefore the continuation payoffs, are independent of the past choice that indexes the columns of |capU0| and |capU1|.
	
The logit assumption only affects the operator $\Psi$, and therefore the function |bellman|, through the specification of the surpluses $R_0$ and $R_1$ in (\ref{eq:surplus}). If you want to change the logit assumption, you should change the computation of |r0| and |r1| (and make sure to adapt the computation of choice probabilities and inverse choice probabilities elsewhere as well).
%}
