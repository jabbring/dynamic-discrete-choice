%{
\documentclass{article}
	\title{Dynamic Discrete Choice Models: Methods, Matlab Code, and Exercises\thanks{We wish to thank Jeffrey R. Campbell for help with \textsc{Komments++} and to Nan Yang and our students for comments on earlier versions of this code. First draft: February 8, 2011. \copyright 2020 Jaap H. Abbring and Tobias J. Klein.}}
	\author{Jaap H. Abbring\thanks{Department of Econometrics \& OR, Tilburg University, P.O. Box
		90153, 5000 LE Tilburg, The Netherlands. E-mail: \url{mailto:J.H.Abbring@uvt.nl}{J.H.Abbring@uvt.nl}. \url{http://jaap.abbring.org/}{http://jaap.abbring.org/}.}}
			\author{Tobias J. Klein\thanks{Department of Econometrics \& OR, Tilburg University, P.O. Box
		90153, 5000 LE Tilburg, The Netherlands. E-mail: \url{mailto:T.J.Klein@uvt.nl}{T.J.Klein@uvt.nl}. \url{http://www.tobiasklein.ws/}{http://www.tobiasklein.ws/}.}}
			\date{This version: November 12, 2020}
	\begin{abstract}
			This document supports the first \textsc{Matlab} computing sessions in our PhD elective course Empirical Industrial Organization 2 in CentER Tilburg's Research Master in Economics program (230323). It contains some notes on the theory of dynamic discrete choice models and on methods for their computation and estimation. It is centered around some basic \textsc{Matlab} code for solving, simulating, and empirically analyzing a simple dynamic discrete choice model. Student exercises ask students to extend this code to apply different and more advanced computational and econometric methods to a wider range of models.
\end{abstract}


\paragraph{Viewing and Using this File}

This file documents the source code that is available from the GitHub repository \url{https://github.com/jabbring/dynamic-discrete-choice}{jabbring/dynamic-discrete-choice}. You can also download a \url{http://jabbring.github.io/dynamic-discrete-choice/dynamicDiscreteChoice.zip}{\textsc{Zip} archive} that includes both this file and all code. We welcome the use of this documentation and code under an MIT license (see the \textsc{LICENSE} file for details). 
			
This file (|dynamicDiscreteChoice.m.html|) was generated from the \textsc{Matlab} script |dynamicDiscreteChoice.m| using the \textsc{Komments++} package, which was created and generously provided to us by Jeffrey R. Campbell. It documents how you can run the script |dynamicDiscreteChoice.m| with \textsc{Matlab} to specify, simulate, and estimate an empirical discrete-time version of the model of firm entry and exit under uncertainty by \cite{jpe89:dixit}. These computations will use, and therefore exemplify the use of, various tailor-made \textsc{Matlab} functions that are documented in this file. Thus, this file also is a guide to these functions and the way they can be adapted and used in your own exercises.
    
In Safari and Firefox, you can switch between the default view of this document, which displays the working code with all its documentation, and an alternative view that shows \emph{only} the code by pressing ``c''.

\paragraph{Software Requirements}
	
Running the code documented in this file requires a recent version of \url{http://www.mathworks.com/products/matlab/}{\textsc{Matlab}} with its \url{http://www.mathworks.nl/help/optim/index.html}{\textsc{Optimization Toolbox}}. The code can easily be adapted to use the \url{https://www.artelys.com/solvers/knitro/}{\textsc{Knitro}} solver instead of \textsc{Matlab}'s \textsc{Optimization Toolbox} (see Section \ref{script}).\footnote{The code has been tested with \textsc{Matlab} 2014a and later, with its \textsc{Optimization Toolbox}, on iMacs with OS X 10.9 and later. A free trial version of \textsc{Knitro} is available from \url{https://www.artelys.com/solvers/knitro/}{Artelys}.}

\paragraph{Road Map}

The remainder of this document proceeds as follows. The next section covers two functions that define the decision problem, |flowpayoffs|  and |bellman|. The versions of these functions handed out to the students, and documented here, define a very basic entry and exit problem with sunk costs and ongoing uncertainty. By changing these functions, different and more extensive decision problems can be studied with the procedures documented in later sections. Section \ref{solveDP} discusses |fixedPoint|, which solves for the optimal choice-specific expected discounted values (and therewith for the optimal decision rule) by finding the unique fixed point of the contraction mapping defined by |bellman|. Section \ref{simulate} documents how |simulateData| can be used to simulate data for the decision problem set up by |flowpayoffs| and |bellman|. Section \ref{estimate} documents functions for estimating the transition matrix $\Pi$ and for computing the log partial likelihood for conditional choice data. Section \ref{script} documents the script that uses these functions to set up the model, simulate data, and estimate the model with \cite{ecta87:rust}'s nested fixed point (NFXP) maximum likelihood method. Section \ref{exercises} contains a range of student exercises. Appendix \ref{contractions} provides mathematical background by discussing some theory of contractions. Appendix \ref{misc} documents procedures that are called by the main routines, but that are of limited substantial interest, such as |randomDiscrete|.

\section{A Firm Entry and Exit Problem\label{DP}}
	
Consider a simple, discrete-time version of the model of firm entry and exit under uncertainty in \cite{jpe89:dixit}. At each time $t\in\mathbb{N}$, a firm can choose to either serve a market (choice 1) or not (choice 0). The payoffs from either choice depend on the firm's choice $A_{t-1}$ in the previous period, because the firm may have to incur a sunk cost to enter or exit the market (for definiteness, suppose that the firm is initially inactive, $A_0=0$). Payoffs also depend on externally determined (choice-independent) observed scalar state variables $X_t$, as well as unobserved state variables $\varepsilon_{t}(0)$ and $\varepsilon_{t}(1)$ that are revealed to the firm before it makes its choice in period $t$. Specifically, in period $t$, its flow profits from choice 0  are
			\[u_0(X_t,A_{t-1})+\varepsilon_{t}(0)= -A_{t-1}\delta_0+\varepsilon_{t}(0),\]		
	and its flow profits from choice 1 are
			\[u_1(X_t,A_{t-1})+\varepsilon_{t}(1)=\beta_0+\beta_1 X_t-(1-A_{t-1})\delta_1+\varepsilon_{t}(1).\]
	Note that $\delta_0$ and $\delta_1$ are exit and entry costs that the firm only pays if it changes state. Gross of these costs, an active firm makes profits $\beta_0+\beta_1 X_t+\varepsilon_{t}(1)$ and an inactive firm has profits $\varepsilon_{t}(0)$.  
	
	The state variable $X_t$ has finite support ${\cal X}\equiv\{x^1,\ldots,x^K\}$. From its random initial value $X_0$, it follows a first-order Markov chain with $K\times K$ transition probability matrix $\Pi$, with typical element $\Pi_{ij}=\Pr(X_{t+1}=x^j|X_t=x^i)$, independently of the firm's choices. The profit shocks $\varepsilon_{t}(a)$ are independent of $\{X_t\}$ and across time $t$ and choices $a$, with type I extreme value distributions centered around 0. Like $X_t$, they may affect but are not affected by the firm's choices.	Thus, $X_t$ and $\varepsilon_t\equiv[\varepsilon_t(0),\varepsilon_t(1)]$ are ``exogenous'' state variables of which the evolution is externally specified. Consequently, the firm controls the evolution of the state $(X_t,A_{t-1},\varepsilon_t)$ only through its choice of $A_{t-1}$.
	
The firm has rational expectations about future states and choices. It chooses the action $A_t$ that maximizes its expected flow of profits, discounted at a factor $\rho < 1$. 
	
Two functions, |flowpayoffs| and |bellman|, together code up this simple model. If you wish to experiment with other functional forms for the flow profits, you should edit |flowpayoffs.m|. The model's dynamic specification can be changed in |bellman.m|.
	 	
\subsection{Flow Payoffs\label{flowpayoffs}}
	\input[4..end]{flowpayoffs.m}
	
\subsection{Bellman Operator\label{bellman}}
	\input[4..end]{bellman.m}

\section{Solving the Decision Problem\label{solveDP}}
	\input[4..end]{fixedPoint.m}

\section{Data Simulation\label{simulate}}
	\input[4..end]{simulateData.m}

\section{Estimation\label{estimate}}
	Suppose that we have a random sample $\left\{\left[(x_{1n},a_{1n}),\ldots,(x_{Tn},a_{Tn})\right]; n=1,\ldots,N\right\}$ from the population of state and choice histories $\{(X_1,A_1),\ldots,(X_T,A_T)\}$. Following \cite{nh94:rust}, we provide functions for implementing a two-stage procedure in which, first, the state transition matrix $\Pi$ is estimated directly from observed state transitions and, second, the remaining parameters are estimated by maximizing the partial likelihood for the observed choices. 	
	
	\subsection{First Stage: State Transitions\label{Pi}}	
	\input[4..end]{estimatePi.m}
	\subsection{Second Stage: Choices\label{likelihood}}
	\input[4..end]{negLogLik.m}

\section{The Script that Puts It All Together\label{script}}

	The script in |dynamicDiscreteChoice.m| simulates data and computed maximum likelihood estimates using the nested fixed point (NFXP) method of \cite{ecta87:rust} and \cite{nh94:rust}. It takes ${\cal
	X},\delta_0,\rho$ to be known, either takes $\Pi$ to be known or estimates it in a first stage, and focuses on maximum partial likelihood estimation of the remaining parameters; $\beta_0,\beta_1,\delta_1$; from conditional choice probabilities. 
	
\subsection{Simulating Data}
	
	First, we set the number of time periods (|nPeriods|) and firms (|nFirms|) that we would like to have in our sample.
%}
nPeriods = 100
nFirms = 1000
%{
	We also set the tolerance |tolFixedPoint| on the fixed point $U$ of $\Psi$ that we will use to determine the simulation's entry and exit rules. This same tolerance will also be used when solving the model in the inner loop of the NFXP procedure.
%}
tolFixedPoint = 1e-10
%{
Next, we specify the values of the model's parameters used in the simulation: 
	\begin{dictionary}
	\item{|nSuppX|} the scalar number $K$ of elements of ${\cal X}$;
	\item{|supportX|} the $K\times 1$ vector ${\cal X}$ with the support points of $X_t$;	
	\item{|capPi|} the $K\times K$ Markov transition matrix $\Pi$ for $\{X_t\}$, with typical element $\Pi_{ij}=\Pr(X_{t+1}=x^j|X_t=x^i)$;
	\item{|beta|} the $2\times 1$ vector $\beta$ with the parameters of the flow profit of active firms;
	\item{|delta|} the $2\times 1$ vector of exit and entry costs $\delta$; and
	\item{|rho|} the scalar discount factor $\rho$.
	\end{dictionary}
	%}													
nSuppX = 5;
supportX = (1:nSuppX)'
capPi = 1./(1+abs(ones(nSuppX,1)*(1:nSuppX)-(1:nSuppX)'*ones(1,nSuppX)));
capPi = capPi./(sum(capPi')'*ones(1,nSuppX))
beta = [-0.1*nSuppX;0.2]
delta = [0;1]
rho = 0.95	
%{
For these parameter values, we compute the flow payoffs $u_0$ (|u0|) and $u_1$ (|u1|), the choice-specific expected discounted values $U_0$ (|capU0|) and $U_1$ (|capU1|), and their contrast $\Delta U$ (|deltaU|).
%}
[u0,u1] = flowpayoffs(supportX,beta,delta); 
[capU0,capU1] = fixedPoint(u0,u1,capPi,rho,tolFixedPoint,@bellman,[],[]);
deltaU = capU1-capU0;
%{
	With $\Delta U$ computed, and $\Pi$ specified, we proceed to simulate a $T\times N$ matrix of choices |choices| and a $T\times N$ matrix of states |iX| (recall from Section \ref{simulate} that |iX| contains indices that point to elements of ${\cal X}$ rather than those values themselves).
%}
[choices,iX] = simulateData(deltaU,capPi,nPeriods,nFirms);
%{
\subsection{Nested Fixed Point Maximum Likelihood Estimation}

First, suppose that $\Pi$ is known. We use |fmincon| from \textsc{Matlab}'s \textsc{Optimization Toolbox} to maximize the partial likelihood for the choices (the code can easily be adapted to use other optimizers and packages, because these have a very similar \url{http://www.mathworks.nl/help/optim/ug/fmincon.html}{syntax}; see below). Because |fmincon| is a minimizer, we use minus the log likelihood as its objective. The function |negLogLik| computes this objective, but has input arguments other than the vector of model parameters to be estimated. Because \url{http://www.mathworks.nl/help/optim/ug/passing-extra-parameters.html}{the syntax of |fmincon| does not allow this}, we define a function handle |objectiveFunction| to an anonymous function that equals |negLogLik| but does not have this extra inputs.
%}
objectiveFunction = @(parameters)negLogLik(choices,iX,supportX,capPi,parameters(1:2),[delta(1);parameters(3)],...
                                           rho,@flowpayoffs,@bellman,@fixedPoint,tolFixedPoint)
%{
Before we can put |fmincon| to work on this objective function, we first have to set some of its other input arguments. We specify a $3\times 1$ vector |startvalues| with starting values for the parameters to be estimated, $(\beta_0,\beta_1,\delta_1)'$.
%}
startvalues = [-1;-0.1;0.5];
%{
    We also set a lower bound of 0 on the third parameter, $\delta_1$, and (nonbinding) lower bounds of $-\infty$ on the other two parameters (|lowerBounds|). There is no need to specify upper bounds.\footnote{Note that |fmincon|, but also its alternatives discussed below, allow the user to specify bounds on parameters; if another function is used that does not allow for bounds on the parameters, you can use an alternative parameterization to ensure that parameters only take values in some admissible set (for example, you can specify $\delta_1=\exp(\delta_1^*)$ for $\delta_1^*\in\mathbb{R}$ to ensure that $\delta_1>0$). Minimizers like |fmincon| also allow you to impose more elaborate constraints on the parameters; you will need this option when implementing the MPEC alternative to NFXP of \cite{ecta12:juddsu} (see Section \ref{exercises}).}
%}
lowerBounds = -Inf*ones(size(startvalues));
lowerBounds(3) = 0;
%{
    Finally, we pass some options, including tolerances that specify the criterion for the outer loop convergence, to |fmincon| through the structure |OptimizerOptions| (recall that we have already set the inner loop tolerance in |tolFixedPoint|). We use the function |optimset| from the \textsc{Optimization Toolbox} to assign values to specific fields (options) in |OptimizerOptions| and then call |fmincon| to run the NFXP maximum likelihood procedure (to use \textsc{Knitro} instead, simply replace |fmincon| by |knitromatlab|, |knitrolink|, or |ktrlink|, depending on the packages installed\footnote{|fmincon| requires \textsc{Matlab}'s \textsc{Optimization Toolbox}, |knitromatlab| is included in \textsc{Knitro} 9.0, |knitrolink| uses both, and |ktrlink| can be used if the \textsc{Optimization Toolbox} is installed with an earlier version of \textsc{Knitro}.}).
%}
OptimizerOptions = optimset('Display','iter','Algorithm','interior-point','AlwaysHonorConstraints','bounds',...
                            'GradObj','on','TolFun',1E-6,'TolX',1E-10,'DerivativeCheck','off','TypicalX',[beta;delta(2)]);
[maxLikEstimates,~,exitflag] = fmincon(objectiveFunction,startvalues,[],[],[],[],lowerBounds,[],[],OptimizerOptions)
%{
This gives maximum partial likelihood estimates of $(\beta_0,\beta_1,\delta_1)$. To calculate standard errors, we call |negLogLik| once more to estimate the corresponding Fisher information matrix and store this in |informationMatrix|. Its inverse is an estimate of the maximum likelihood estimator's asymptotic variance-covariance matrix.
%}
[~,~,informationMatrix] = objectiveFunction(maxLikEstimates);
standardErrors = diag(sqrt(inv(informationMatrix)));
%{
The resulting parameter estimates and standard errors are displayed (third and fourth columns), together with the parameters' true (first column) and starting values (second column).
%}
disp('Summary of Results');
disp('--------------------------------------------');
disp('      true     start     estim      ste.');
disp([[beta;delta(2)] startvalues maxLikEstimates standardErrors]);
%{
\subsection{Extension to an Unknown Markov Transition Matrix for the State}
Finally, consider the more realistic case that $\Pi$ is not known. In this case, \cite{nh94:rust} suggests a two-stage procedure. In the first stage, we estimate $\Pi$ using |estimatePi| and store the results in a $K\times K$ matrix |piHat|.
%}
piHat = estimatePi(iX,nSuppX)
%{ 
In the second stage, maximum partial likelihood estimates of $(\beta_0,\beta_1,\delta_1)$ can be computed using the NFXP procedure, with |piHat| replacing |capPi|. This, with the question how the first-stage sampling error affects the precision of the second-stage estimator of $(\beta_0,\beta_1,\delta_1)$, is left for the exercises.

\section{Student Exercises\label{exercises}}
	
\subsection{Theoretical Exercises}

    \paragraph{Theoretical Exercise 1} Prove that $\Psi$ is a contraction.
    
    \paragraph{Theoretical Exercise 2} Is $\rho$ identified in the example model? What about $\delta_0$? (To the extent that these are identified, you may want to extend your numerical experiments below to include their estimation.) See \cite{are10:abbring} for results and references.

    \paragraph{Theoretical Exercise 3} Show the equivalence of policy iteration
    and the Newton-Kantorovich method for solving $U-\Psi(U)=0$ for $U$ (see
    Appendices \ref{contractions.computation} and \ref{contractions.DP}). Read about the 
    properties of policy iteration and discuss what these imply about the
    convergence properties of the Newton-Kantorovich method. Does
    policy iteration converge in a finite number of steps? Why (not)?
		
	
\subsection{Computational Exercises}

The computational exercises typically ask you to modify the code so that it can handle alternative econometric procedures and models. This often requires that you adapt the analytic computation of the derivatives of the log likelihood. It may be convenient to first (or only) work with numerical derivatives by setting the option "GradObj" in |fmincon| or its alternative to "off"; this allows you to find the estimates without using the analytical gradient in the optimization step.
Once you have coded up the analytical derivatives (in order to use them in the optimization step and to obtain an estimate of the standard errors), you are advised to verify them against numerical derivatives, for example by setting the option "DerivativeCheck" in |fmincon| to "on" (however, make sure to set this option to "off" once you have verified the derivatives and start using the code for repeated, e.g. Monte Carlo, calculations).

    \paragraph{Computational Exercise 1} Extend the script in |dynamicDiscreteChoice.m| with a single simulation to a Monte Carlo experiment in which estimates are computed for a sequence of simulated data sets (for now, as in |dynamicDiscreteChoice.m|, take $\Pi$ to be known). Keep track of the estimates and report statistics such as their means and standard deviations. Compare the latter to asymptotic standard errors. Make sure to set the seed of the random number generator so that your numerical results can be replicated. Note that you can also use this Monte Carlo setup in later questions to study the finite-sample performance of the various procedures studied.

    \paragraph{Computational Exercise 2} For the case of estimating the demand for differentiated products using NFXP maximum likelihood, as in \cite{ecta95:berryetal}, \cite{ecta12:dubeetal} claim that it is important to carefully balance the convergence tolerances used in the inner (contraction) loop and the outer (likelihood maximization) loop. In particular, they argue that the inner loop needs to be computed at a much higher precision than the tolerance used in the outer loop. Experiment with the values of |tolFixedPoint| on the one hand and the tolerances in the |tolX| and |tolFun| fields of |OptimizerOptions| on the other hand to investigate this issue in the context of this package's firm entry and exit problem.

    \paragraph{Computational Exercise 3} Implement the MPEC approach to maximum likelihood estimation of our structural model, as advocated by \cite{ecta12:juddsu}.
		\begin{itemize} 
		\item First, modify |negLogLik| so that it it takes $(U_0,U_1)$ or $\Delta U$ as an input argument (instead of solving the model for them) and computes the log (partial) likelihood directly from the choice probabilities implied by $\Delta U$.
        \item Then, extend the script in |dynamicDiscreteChoice.m| so that it alternatively maximizes the log likelihood with respect to both the parameters of interest and the values of $\Delta U$, subject to the constraint on $\Delta U$ implied by $U=\Psi(U)$ (which can be specified using the function |bellman|).
		\end{itemize}
		Would you expect the NFXP and MPEC approaches to give the same estimates
		of the parameters of interest (up to numerical precision)? How is the
		relative numerical performance of both procedures? Which one is faster?
		Compare your results to those in \cite{ecta15:iskhakovetal} and explain.
    
    \paragraph{Computational Exercise 4} Compute the two-stage maximum partial likelihood estimates of $(\beta_0,\beta_1,\delta_1)$ for the case that $\Pi$ is not known. Are the estimators of the standard errors that assume $\Pi$ known consistent if $\Pi$ is estimated instead? Read \cite{nh94:rust} and explain why the two-stage estimator of $(\beta_0,\beta_1,\delta_1)$ is not efficient. The full information maximum likelihood (FIML) estimator estimates all parameters, including the ones in $\Pi$, simultaneously by maximizing the full likelihood for the observed state transitions and choices. Write an alternative for |negLogLik| that codes up this likelihood function and adapt the estimation script so that it obtains the FIML estimates and their estimated standard errors. Do not code up the gradient and do not use it in the likelihood maximization. Perform a Monte Carlo study of the two-stage and FIML estimators of $(\beta_0,\beta_1,\delta_1)$. For both estimators, report the means and standard deviations of the estimates of $(\beta_0,\beta_1,\delta_1)$ and the average estimates of the standard errors across Monte Carlo simulations. Compare and discuss. Finally, discuss the three-stage approach suggested by \cite{nh94:rust}.

    \paragraph{Computational Exercise 5} Implement a simulation-based two-step estimator along the lines of \cite{res94:hotzetal}, \cite{nh94:rust}, and \cite{ecta07:bajarietal}.
		\begin{itemize} 
		\item Add a new function for nonparametrically estimating the choice probabilities $p(a|X_t,A_{t-1})\equiv\Pr(A_t=a|X_t,A_{t-1})$ over the support of $(X_t,A_{t-1})$. With enough data and a finite state space, you can simply use the conditional sample frequency. 
		\item Add (to this same function or in a new function) a procedure for estimating $\Delta U$ by inverting the estimated choice probabilities (as proposed by \cite{res93:hotzmiller} and \cite{res94:hotzetal}). 
		\item Extend |negLogLik| (or add functions) so that it optionally takes nonparametric estimates of $\Delta U$ as inputs, computes the corresponding estimates of the entry and exit rules, and uses these estimates and the model to forwardly simulate $U_0(x,a)$ and $U_1(x,a)$, and then $\Delta U(x,a)$, for each point $(x,a)$ in the sample. As possible objectives to be minimized, both implement a weighted distance between the nonparametric estimates of $\Delta U$ and the simulated values of $\Delta U$, and minus a log pseudo-likelihood based on the choice probabilities implied by the simulated values of $\Delta U$.
			\item Finally, extend the script in |dynamicDiscreteChoice.m| so that it successively calls these functions to implement a two step procedure for estimating the model, as in \cite{res93:hotzmiller} and \cite{res94:hotzetal}. Analyze the numerical and statistical performance of this procedure with both of the possible objective functions implemented and compare. Discuss the theoretical relation between both approaches; see \cite{res08:pesendorferschmidtdengler}.
		\end{itemize}
	See the course slides for a brief and general description of this procedure and \cite{nh94:rust} for a detailed algorithm, with discussion. Note that the algorithm described in \cite{nh94:rust} is similar to that of \cite{ecta07:bajarietal} for games that we will discuss later in the course. Both build on the ideas of \cite{res94:hotzetal}, who use the special logit structure to simplify the simulation (see the discussion in \cite{ecta07:bajarietal}). 

    \paragraph{Computational Exercise 6} Extend the code for simulation and NFXP estimation so that it can handle a finite number of unobserved types, as in the work of \cite{Eckstein_Wolpin_1990_Econometrica}, \cite{Keane_Wolpin_1997_JPE}, and \cite{Eckstein_Wolpin_1999_Econometrica}. Specifically, suppose that there are two types in the population, one with entry cost $\delta_1^1$ and the other with entry cost $\delta_1^2$. The firms are distributed across both unobserved entry cost types independently from all other variables in the model. We would like to estimate $\beta_0,\beta_1,\delta_1^1,\delta_1^2$, and the share of agents with entry cost $\delta_1^1$. To this end:
		\begin{itemize} 
		\item Extend |simulateData| so that it randomly draws an entry cost from a distribution with two points of support, $\delta_1^1$ and $\delta_1^2$, and simulates choice data corresponding to the entry cost drawn, for each firm, independently across firms. 
		\item Extend |negLogLik| so that it computes each firm's likelihood contribution as the mixture (expectation) of the likelihood contributions for each entry cost type (which are given by the current specification of the likelihood for homogeneous firms) over the distribution of these two types in the population of firms.
		\item Finally, extend the script in |dynamicDiscreteChoice.m| so that it successively calls these functions. Check whether you can recover the entry cost distribution, and the other parameters, well.  
		\end{itemize}
    
    \paragraph{Computational Exercise 7} Implement a version of the two-step estimator of \cite{res94:hotzetal} that similarly allows for a finite number of unobserved types. To this end, combine this estimator with the EM algorithm as in \cite{ecta11:arcidiaconomiller} (see also \cite{are11:arcidiaconoellickson}).

    \paragraph{Computational Exercise 8} Perform a Monte Carlo study of all the estimators that you have implemented. Evaluate and compare their numerical performance as measured in terms of convergence success and computation time; and their finite sample statistical performance in terms of, among other things, mean squared errors. Experiment with different choices of the numerical design parameters, such as the inner and outer loop tolerances used in the NFXP procedure.

    \paragraph{Computational Exercise 9} Extend |fixedPoint| so that it has the option to compute $U$ using the Newton-Kantorovich method instead of successive approximations, or a hybrid method that combines both (as in \cite{ecta87:rust}). Note that, with a finite state space as in our example, the Newton-Kantovorich method reduces to the Newton-Raphson method and requires the Jacobian $\partial\tilde\Psi_\theta(\bar U)/\partial\bar U'$ computed in an intermediate step of |negLogLik|. Explore the numerical performance of the various methods. What is your preferred method?


    
\begin{bibliography}
    \bibitem[Abbring (2010)]{are10:abbring} Abbring, Jaap H. (2010). Identification of dynamic discrete choice models. \textit{Annual Review of Economics} 2, 367-394.
    \bibitem[Aguirregabiria and Mira (2002)]{ecta02:aguirregabiriamira} Aguirregabiria, Victor and Pedro Mira (2002). Swapping the nested fixed point algorithm: A class of estimators for discrete Markov decision models. \textit{Econometrica} 70(4), 1519-1543.  
    \bibitem[Arcidiacono and Ellickson (2011)]{are11:arcidiaconoellickson} Arcidiacono, Peter and Paul B. Ellickson (2011).  Practical methods for estimation of dynamic discrete choice models. \textit{Annual Review of Economics} 3, 363-394.
    \bibitem[Arcidiacono and Miller (2011)]{ecta11:arcidiaconomiller} Arcidiacono, Peter and David A. Miller (2011).  Conditional choice probability estimation of dynamic discrete choice models with unobserved heterogeneity. \textit{Econometrica} 79, 1823-1867.
    \bibitem[Bajari, Benkard, and Levin (2007)]{ecta07:bajarietal}Bajari, Patrick, Lanier Benkard, and Jonathan Levin (2007). Estimating dynamic models of imperfect competition. \textit{Econometrica} 75(5), 1331-1370.
    \bibitem[Berry, Levinsohn, and Pakes (1995)]{ecta95:berryetal} Berry, Steven, James Levinsohn, and Ariel Pakes (1995). Automobile prices in market equilibrium. \textit{Econometrica} 63, 841-890.
    \bibitem[Dixit (1989)]{jpe89:dixit}Dixit, Avinash K. (1989). Entry and exit decisions under uncertainty. \textit{Journal of Political Economy} 97(3), 620-638.
    \bibitem[Dub&eacute;, Fox, and Su (2012)]{ecta12:dubeetal}Dub&eacute;, Jean-Pierre, Jeremy T. Fox, and Che-Lin Su (2012). Improving the numerical performance of static and dynamic aggregate discrete choice random coefficients demand estimation. \textit{Econometrica} 80, 2231-2267.
    \bibitem[Eckstein and Wolpin (1990)]{Eckstein_Wolpin_1990_Econometrica} Eckstein, Zvi and Kenneth I. Wolpin (1990). Estimating a market equilibrium search model from panel data on individuals. \textit{Econometrica} 58(4), 783-808.
    \bibitem[Eckstein and Wolpin (1999)]{Eckstein_Wolpin_1999_Econometrica} Eckstein, Zvi and Kenneth I. Wolpin (1999). Why youths drop out of high school: The impact of preferences, opportunities, and abilities. \textit{Econometrica} 67(6), 1295-1339.
    \bibitem[Hotz and Miller (1993)]{res93:hotzmiller}Hotz, V. Joseph and David A. Miller (1993). Conditional choice probabilities and the estimation of dynamic models. \textit{Review of Economic Studies} 60, 497-529.
    \bibitem[Hotz et al. (1994)]{res94:hotzetal}Hotz, V. Joseph, David A. Miller, Seth Sanders and Jeffrey Smith (1994). A simulation estimator for dynamic models of discrete choice. \textit{Review of Economic Studies} 61(2), 265-289.
    \bibitem[Iskhakov et al. (2015)]{ecta15:iskhakovetal} Iskhakov, Fedor, Jinhyuk Lee, John Rust, Bertel Schjerning and Kyoungwon Seo (2015). Constrained optimization approaches to estimation of structural models: Comment. \textit{Econometrica}, forthcoming.
    \bibitem[Judd (1998)]{mit98:judd}Judd, Kenneth L. 1998. \textit{Numerical Methods in Economics}.  MIT Press. Cambridge, MA.
    \bibitem[Judd and Su (2012)]{ecta12:juddsu} Judd, Kenneth L. and Che-Lin Su (2012). Constrained optimization approaches to estimation of structural models. \textit{Econometrica} 80(5), 2213-2230.
    \bibitem[Keane and Wolpin (1997)]{Keane_Wolpin_1997_JPE} Keane, Michael P. and Kenneth I. Wolpin (1997). The career decisions of young men. \textit{Journal of Political Economy} 105(3), 473-522.
    \bibitem[Ljungqvist and Sargent (2000)]{mit00:ljungqvistsargent} Ljungqvist, Lars and Thomas J. Sargent (2000). \textit{Recursive Macroeconomic Theory}, \url{http://www.econ.yale.edu/smith/econ510a/sargent3.pdf}{Second edition}. MIT Press. Cambridge, MA.
    \bibitem[McFadden (1981)]{McFadden_1981_probabilistic} McFadden, Daniel (1981). Econometric models of probabilistic choice. In C. Manski and D. McFadden (Eds.). \textit{Discrete Data with Econometric Applications}. MIT Press. Cambridge, MA.
    \bibitem[Pesendorfer and Schmidt-Dengler (2008)]{res08:pesendorferschmidtdengler} Pesendorfer, Martin and Philipp Schmidt-Dengler (2008). Asymptotic least squares estimators for dynamic games. \textit{Review of Economic Studies} 75(3), 901-928.
    \bibitem[Puterman and Brumelle (1979)]{mor79:putermanbrumelle} Puterman, Martin L. and Shelby L. Brumelle (1979). On the convergence of policy iteration in stationary dynamic programming. \textit{Mathematics of Operations Research} 4.1, 60-69.
    \bibitem[Rust (1987)]{ecta87:rust} Rust, John (1987). Optimal replacement of GMC bus engines: An empirical model of Harold Zurcher. \textit{Econometrica} 55, 999-1033.
    \bibitem[Rust (1994)]{nh94:rust} Rust, John (1994). Structural estimation of Markov decision processes. In R. Engle and D. McFadden (Eds.). \textit{Handbook of Econometrics} 4, 3081-3143. North-Holland. Amsterdam.
    \bibitem[Stokey, Lucas, and Prescott (1989)]{harvard1989StokeyLucasPrescott} Stokey, Nancy L. and Robert E. Lucas (with Edward C. Prescott) (1989). \textit{Recursive Methods in Economic Dynamics}. Harvard University Press. Cambridge, MA.
\end{bibliography}
	
	
\appendix
\section{Contractions and Dynamic Programming\label{contractions}}
	
See e.g. \cite{harvard1989StokeyLucasPrescott} and \cite{mit98:judd} for a rigorous and comprehensive treatment of the topics in this appendix and their economic applications.

\subsection{Definitions\label{contractions.definitions}}

A \emph{metric space} is a set $\cal{U}$ with a metric $m:\cal{U}\times\cal{U}\rightarrow\mathbb{R}$ such that, for any $U,U',U''\in\cal{U}$,
\begin{enumerate}
	\item $m(U,U')=0$ if and only if $U=U'$,
	\item $m(U,U')=m(U',U)$, and 
	\item $m(U,U'')\leq m(U,U')+m(U',U'')$ (triangle inequality).
\end{enumerate}

A \emph{Cauchy sequence} in $({\cal U},m)$ is a sequence $\{U_n\}$ in ${\cal U}$ such that, for each $\epsilon>0$, there exists some $n(\epsilon)\in\mathbb{N}$ such that $m(U_n,U_{n'})<\epsilon$ for all $n,n'\geq n(\epsilon)$. 
		
A metric space $({\cal U},m)$ is \emph{complete} if every Cauchy sequence in $({\cal U},m)$  converges to a point in ${\cal U}$.
		
A map $\Psi:{\cal U}\rightarrow{\cal U}$ is a \emph{contraction with modulus $\rho\in(0,1)$} if $m\left[\Psi(U),\Psi(U')\right]\leq \rho m(U,U')$ for all $U,U'\in{\cal U}$.
			 
\subsection{Contraction Mapping (Banach Fixed Point) Theorem\label{contractions.Banach}}
																										
If $({\cal U},m)$ is a complete metric space and $\Psi:{\cal U}\rightarrow{\cal U}$ is a contraction, then there exists a unique $U\in{\cal U}$ such that $U=\Psi(U)$. Furthermore, for any $U_0\in{\cal U}$, the sequence $\{U_n\}$ with $U_{n+1}=\Psi(U_{n})$, $n\in\mathbb{N}$, converges to $U$.
																															
\paragraph{Sketch of Proof}
Suppose $U,U'\in{\cal U}$ are such that $U=\Psi(U)$ and $U'=\Psi(U')$. Then, $0\leq m(U,U')=m\left[\Psi(U),\Psi(U')\right]\leq\rho m(U,U')$ for some $\rho\in(0,1)$. Consequently, $m(U,U')=0$ and $U=U'$.

Because $\Psi$ is a contraction, any  $\{U_n\}$ it generates is a Cauchy sequence. Because $({\cal U},m)$ is complete, it converges to some $U\in{\cal U}$. Because contractions are continuous, $U=\Psi(U)$.

\subsection{Computing the Fixed Point of a Contraction\label{contractions.computation}}
The method of \emph{successive approximations} directly applies the Contraction Mapping Theorem and approximates the fixed point $U$ of a contraction $\Psi$ with the sequence $\{U_n\}$  generated from $U_{n+1}=\Psi(U_n)$ on a finitely discretized state space. This method has global but linear convergence.

The \emph{Newton-Kantorovich} method searches for a zero of $I-\Psi$, where $I:U\in{\cal U}\mapsto U$ is the identity mapping. It approximates $U$ with the sequence generated from \[U_{n+1}=U_n-\left[I-\Psi'_{U_n}\right]^{-1}\left[U_n-\Psi(U_n)\right],\] with $I-\Psi'_{U_n}$ the Fr&eacute;chet derivative of $I-\Psi$ at $U_n$ (with a discrete state space, this is simply the
    linear map defined by the appropriate finite dimensional Jacobian, as in |negLogLik|).
   
\subsection{Blackwell's Sufficient Conditions for a Contraction\label{contractions.Blackwell}}

Let ${\cal U}$ be the space of functions $U:{\cal X}\rightarrow\mathbb{R}$ such that $\sup|U|<\infty$ ($U$ bounded), with metric $m(U,U')=\sup|U-U'|$. Suppose that, for some $\rho\in(0,1)$, $\Psi:{\cal U}\rightarrow{\cal U}$ satisfies
\begin{enumerate}
\item (monotonicity) $\Psi(U)\leq \Psi(U')$ and
\item (discounting) $\Psi(U+\gamma)\leq \Psi(U)+\rho\gamma$
\end{enumerate}
for all $U,U'\in{\cal U}$ such that $U\leq U'$ and all $\gamma\in(0,\infty)$. Then, $\Psi$ is a contraction with modulus $\rho$.

\paragraph{Sketch of Proof}
Let $U,U'\in{\cal U}$ be such that  $U\leq U'$. Note that $U\leq U'+m(U,U')$, so that $\Psi(U)\leq \Psi(U')+\rho m(U,U')$ by monotonicity and discounting. Similarly,   $\Psi(U')\leq \Psi(U)+\rho m(U,U')$. Thus, $m\left[\Psi(U),\Psi(U')\right]\leq\rho m(U,U')$.

 \subsection{Application to Dynamic Programming\label{contractions.DP}} 

The choice-specific value function $U$ is a fixed point of $\Psi$, which satisfies Blackwell's sufficient conditions with $\rho$ equal to the discount factor. 
\begin{itemize}
\item Under regularity conditions such that $U$ is bounded,  this ensures that it is the unique fixed point of a contraction on the complete (Banach) space $({\cal U},m)$ of bounded functions. 
\item If ${\cal U}'\subset{\cal U}$ is closed and  $\Psi(U')\in{\cal U}'$ for all $U'\in{\cal U}'$, then $\Psi$ is a contraction on the complete subspace $({\cal U}',m)$, so that $U\in{\cal U}'$ (examples: monotonicity, continuity, etcetera).
\item $U$ can be computed by successive approximations, the Newton-Kantorovich method, or a hybrid algorithm as in \cite{ecta87:rust}.
\end{itemize}
In this context, the method of successive approximations is alternatively
referred to as \emph{value iteration}. Moreover, the Newton-Kantorovich method
is closely related to an alternative method called \emph{policy iteration}. Each of this method's iterations takes a value function as input, determines the corresponding
optimal policy, and then gives the values of applying that policy forever as
output (see e.g. \cite{mit00:ljungqvistsargent}, Section 3.1.1, which refers to
this method as \emph{Howard's policy improvement algorithm}, or \cite{nh94:rust}, Section 2.5). 

\cite{mor79:putermanbrumelle} show that policy iteration is generally
equivalent to applying the Newton-Kantorovich method to finding the fixed point
of the Bellman equation (in the sense that both produce the same sequence of value functions when starting from the same value). \cite{mit00:ljungqvistsargent}, Section 4.4, present and develop this result for the special case of a finite state space. 

Finally, for our example's special case of a finite ${\cal X}$,
\cite{ecta02:aguirregabiriamira}'s Proposition 1 establishes that policy
iteration is equivalent to applying the Newton-Kantorovich method to finding a
fixed point of the Bellman-like operator $\Psi$ (which \cite{ecta02:aguirregabiriamira} refer to as a \emph{smoothed} Bellman operator).

\section{Miscellaneous Utilities\label{misc}}
																																	
This section documents some miscellaneous utilities that are called by the main functions. At this point, it only includes |randomDiscrete|.
																																	
\input[4..end]{randomDiscrete.m}
																										

%}
