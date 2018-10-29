%% negLogLik.m
% The function computes minus the log partial likelihood and, optionally, the corresponding minus score and information matrix estimate for conditional choice data on the basic firm entry and exit model used as an example in CentER's Empirical Industrial Organization II.

%{
The function |negLogLik| computes minus the log partial likelihood for the conditional choice part of the data. Optionally, it also returns minus the corresponding score vector and an estimate of the information matrix for the parameter (sub)vector $\theta\equiv(\beta_0,\beta_1,\delta_1)'$ (the scores are specific to the estimation example in Section \ref{script}'s script and should be adapted for inference on other parameters).
%}
function [nll,negScore,informationMatrix] = ...
         negLogLik(choices,iX,supportX,capPi,beta,delta,rho,flowpayoffs,bellman,fixedPoint,tolFixedPoint)
%{
The function |negLogLik| requires the following input arguments:
	\begin{dictionary}
	\item{|choices|} a $T\times N$ matrix with choice observations $a_{tn}$;
	\item{|iX|} a $T\times N$ matrix with indices of observed states $x_{tn}$ in ${\cal X}$ (for example, if $x_{11}=x^3$, then the first element of |iX| is 3, not $x^3$); 
	\item{|supportX|} a $K\times 1$ vector with the support points of the profit state $X_t$ (the elements of $\cal{X}$, consistently ordered with the Markov transition matrix $\Pi$);
	\item{|capPi|} the (possibly estimated) $K\times K$ Markov transition matrix $\Pi$ for $\{X_t\}$, with typical element $\Pi_{ij}=\Pr(X_{t+1}=x^j|X_t=x^i)$;
	\item{|beta|} a $2\times 1$ vector that contains the intercept ($\beta_0$) and profit state slope ($\beta_1$) of the flow payoffs to choice $1$;
	\item{|delta|} a $2\times 1$ vector that contains the firm's exit ($\delta_0$) and entry ($\delta_1$) costs;	
	\item{|rho|} a scalar with the value of the discount factor $\rho$;
	\item{|flowpayoffs|} a handle of a function |[u0,u1]=flowpayoffs(supportX,beta,delta)| that computes the mean flow payoffs $u_0$ and $u_1$;
	\item{|bellman|} a handle of a function |[capU0,capU1] = bellman(capU0,capU1,u0,u1,capPi,rho)| that iterates once on $\Psi$;
	\item{|fixedPoint|} a handle of a function |[capU0,capU1] = fixedPoint(u0,u1,capPi,rho,tolFixedPoint,bellman,capU0,capU1)| that computes the fixed point $U$ of $\Psi$; and
	\item{|tolFixedPoint|} a scalar tolerance level that is used to determine convergence of the successive approximations of the fixed point $U$ of $\Psi$.
	\end{dictionary}
	It returns
	\begin{dictionary}
	\item{|nll|} a scalar with minus the log partial likelihood for the conditional choices
	\end{dictionary}
    and optionally
	\begin{dictionary}
	\item{|negScore|} a $3\times 1$ vector with minus the partial likelihood score for $\theta$ and
    \item{|informationMatrix|} a $3\times 3$ matrix with the sum of the $N$ outer products of the individual contributions to the score for $\theta$.
	\end{dictionary}
The function |negLogLik| first stores the number $K$ of elements of |supportX| in a scalar |nSuppX|.
%}
nSuppX = size(supportX,1);
%{
Next, it computes the flow payoffs $u_0$ (|u0|) and $u_1$ (|u1|), the choice-specific net expected discounted values $U_0$ (|capU0|) and $U_1$ (|capU1|), their contrast $\Delta U$ (|deltaU|), and the implied probabilities $1/\left[1+\exp(\Delta U)\right]$ of not serving the market (|pExit|) for the inputted parameter values. Note that this implements the NFXP procedure's inner loop.
%}
[u0,u1] = flowpayoffs(supportX,beta,delta);
[capU0,capU1] = fixedPoint(u0,u1,capPi,rho,tolFixedPoint,bellman,[],[]);
deltaU = capU1-capU0;
pExit = 1./(1+exp(deltaU));
%{
\paragraph{Log Partial Likelihood}
The contribution to the likelihood of firm $n$'s choice in period $t$ is the conditional choice probability 
	\[p(a_{tn}|x_{tn},a_{(t-1)n})=a_{tn}+\frac{1-2 a_{tn} }{1+\exp\left[\Delta U(x_{tn},a_{(t-1)n})\right]},\] 
with $a_{0n}=0$. The function |negLogLik| first computes these probabilities for each firm $n$ and period $t$ and stores them in a $T\times N$ matrix |p|. Then, it returns minus the sum of their logs, the log partial likelihood for the conditional choices, in |nll|. 
%}
laggedChoices = [zeros(1,size(choices,2));choices(1:end-1,:)];
p = choices + (1-2*choices).*pExit(iX+nSuppX*laggedChoices);
nll = -sum(sum(log(p)));
%{
\paragraph{Score}
If two or more output arguments are demanded from |negLogLik|, it computes and returns minus the partial likelihood
score for $\theta$ (the derivative of minus the log partial likelihood with respect to $\theta$), in |negScore|.
%}
if nargout>=2
    %{
    Firm $n$'s contribution to the score equals
        \[\frac{\partial\log\left[\prod_{t=1}^Tp(a_{tn}|x_{tn},a_{(t-1)n})\right]}{\partial \theta} 
            = -\sum_{t=1}^T\left(1-2 a_{tn}\right)\left[1-p(a_{tn}|x_{tn},a_{(t-1)n})\right] \frac{\partial\Delta U(x_{tn},a_{(t-1)n})}{\partial\theta}.
        \]
    Its calculation requires that we compute $\partial\Delta U/\partial\theta$. Recall that $U$, and therewith $\Delta U$, is only implicitly given by $U=\Psi(U)$. Note that $U=(U_0,U_1)$ is defined on a set with $2K$ points, so that $U$ can be represented by a $4K\times 1$ vector and $\Psi$ by a mapping from $\mathbb{R}^{4K}$ into $\mathbb{R}^{4K}$. Specifically, $U$ can be represented by the $4K\times 1$ vector $\bar U$ that lists the values of $U_0$ and $U_1$ on their domain,
        \[
        \bar U=\left[U_0(x^1,0),\ldots,U_0(x^K,0),U_0(x^1,1),\ldots,U_0(x^K,1),U_1(x^1,0),\ldots,U_1(x^K,0),U_1(x^1,1),\ldots,U_1(x^K,1)\right]',
        \] 
    and that satisfies		
        \[\bar U= \tilde\Psi_\theta(\bar U),\]		
    with $\tilde\Psi_\theta:\mathbb{R}^{4K}\rightarrow\mathbb{R}^{4K}$ an appropriately rearranged version of $\Psi$ (note that we made its dependence on $\theta$ explicit). With this alternative representation of $U$ and $\Psi$ in place, we can solve 
        \[\left[I_{4K} - \frac{\partial \tilde \Psi_\theta(\bar U)}{\partial\bar U'}\right]\frac{\partial\bar U}{\partial \theta'} = \frac{\partial \tilde\Psi_\theta(\bar U)}{\partial\theta'}\]
    for $\partial\bar U/\partial\theta'$ (see also
    \cite{nh94:rust}, p.3110), where $I_i$ is a $i\times i$ unit matrix; 
        \[\frac{\partial \tilde\Psi_\theta(\bar U)}{\partial\bar U'}=
            \left(\begin{array}{cccc}
			  ~D_{00}~	&~O_K~	&~D_{10}~	&~O_K~\\												 
			  ~D_{00}~	&~O_K~	&~D_{10}~	&~O_K~\\
			  ~O_K~		&~D_{01}~	&~O_K~	&~D_{11}~\\												 
			  ~O_K~		&~D_{01}~	&~O_K~	&~D_{11}~												 
			\end{array}\right),
        \]
    with	
        \[D_{a'a}=\rho\Pi\left[\begin{array}{cccc}
						  p(a'|x^1,a)	&0			&\cdots            &0\\
						  0				&p(a'|x^2,a)&				   &\vdots\\
						  \vdots		&			&\ddots            &0\\
						  0				&\cdots		&0                 &p(a'|x^K,a)
						  \end{array}\right],
        \]
    and $O_i$ is a $i\times i$ matrix of zeros; and 
        \[\frac{\partial \tilde \Psi_\theta(\bar U)}{\partial\theta'}=
        \left[
			\begin{array}{cccc}
			~0~ & ~0~ & ~0~\\
			&\cdot&\\
			&\cdot&\\
			&\cdot&\\
			~0~ & ~0~ & ~0~\\
    		~0~ & ~0~ & ~0~\\
			&\cdot&\\
			&\cdot&\\
			&\cdot&\\
			~0~ & ~0~ & ~0~\\
			~1~ & ~x^1~ & ~-1~\\
			&\cdot&\\
			&\cdot&\\
			&\cdot&\\
			~1~ & ~x^K~ & ~-1~\\
			~1~ & ~x^1~ & ~0~\\
			&\cdot&\\
			&\cdot&\\
			&\cdot&\\
			~1~ & ~x^K~ & ~0~
			\end{array}
			\right].\]
    The function |negLogLik| first computes $D_{00}$ (|d00|), $D_{01}$ (|d01|), $D_{10}$ (|d10|), $D_{11}$ (|d11|) and constructs $\partial \tilde\Psi_\theta(\bar U)/\partial\bar U'$ (|dPsi_dUbar|).
    %}
    d00 = rho*capPi*diag(pExit(:,1));
    d01	= rho*capPi*diag(pExit(:,2));
    d10	= rho*capPi-d00;
    d11	= rho*capPi-d01;
    dPsi_dUbar = [[d00;d00;zeros(2*nSuppX,nSuppX)] [zeros(2*nSuppX,nSuppX);d01;d01] ...
                  [d10;d10;zeros(2*nSuppX,nSuppX)] [zeros(2*nSuppX,nSuppX);d11;d11]];
    %{
	It then computes $\partial\tilde\Psi_\theta(\bar U)/\partial\theta'$ (|dPsi_dTheta|; 
    note that the following line is the only code in |negLogLik| that needs to be changed if other parameters than those in $\theta$ are estimated). 
    %}
	dPsi_dTheta = [[zeros(2*nSuppX,1);ones(2*nSuppX,1)] [zeros(2*nSuppX,1);supportX;supportX] [zeros(2*nSuppX,1);-ones(nSuppX,1);zeros(nSuppX,1)]];
    %{
    Next, it computes $\partial\bar U/\partial\theta'$ (|dUbar_dTheta|) and $\partial\Delta U/\partial\theta'$ (|dDeltaU_dTheta|).
    %}
	dUbar_dTheta   = (eye(4*nSuppX)-dPsi_dUbar)\dPsi_dTheta;
	dDeltaU_dTheta    = dUbar_dTheta(2*nSuppX+1:4*nSuppX,:)-dUbar_dTheta(1:2*nSuppX,:);	
    %{
    Finally, it computes the $1\times 3$ vector $-\partial\log\left[\prod_{t=1}^T p(a_{tn}|x_{tn},a_{(t-1)n})\right]/\partial \theta'$ for each $n$, stacks these individual (minus) score contributions in the $N\times 3$ matrix |negFirmScores|, and sums them to compute minus the score vector, |negScore|.   
    %}			
    nTheta	= size(dUbar_dTheta,2);    
    negFirmScores = repmat((1-2*choices).*(1-p),[1 1 nTheta]);
    for i=1:nTheta
        negFirmScores(:,:,i) = negFirmScores(:,:,i).*dDeltaU_dTheta(iX+nSuppX*laggedChoices+2*(i-1)*nSuppX);
    end
    negFirmScores=squeeze(sum(negFirmScores,1));
    negScore = sum(negFirmScores)';
end
%{
\paragraph{Information Matrix}
If three output arguments are demanded, |negLogLik| computes the sum of the $N$ outer products of the individual score contributions, 
\[
	  \sum_{n=1}^N\frac{\partial\log\left[\prod_{t=1}^Tp(a_{tn}|x_{tn},a_{(t-1)n})\right]}{\partial\theta}\cdot\frac{\partial\log\left[\prod_{t=1}^Tp(a_{tn}|x_{tn},a_{(t-1)n})\right]}{\partial\theta'}, 
\]
and returns it in |informationMatrix|.
%}
if nargout==3
    informationMatrix = zeros(nTheta,nTheta);
    for n=1:size(negFirmScores,1)
        informationMatrix = informationMatrix + negFirmScores(n,:)'*negFirmScores(n,:);
    end
end
%{
When evaluated at an estimate of $\theta$, this provides an estimate of the expected partial likelihood information matrix for $\theta$. This estimate can in turn be used to estimate the variance-covariance matrix, and thus the standard errors, of the maximum partial likelihood estimator $\hat\theta$ of $\theta$. 
%}