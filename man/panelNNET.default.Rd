\name{panelNNET.default}
\alias{panelNNET.default}
\alias{panelNNET}
\alias{panelNNET.est}
\alias{panelNNET-package}
\alias{plot.panelNNET}
\alias{summary.panelNNET}
\alias{print.panelNNET}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Semi-parametric fixed-effects models for panel data, using neural networks
}
\description{
Function fits fixed effects models in which some variables entering linearly and others entering non-parametrically.  Nonparametric component of the model represented by a feed-forward neural network.  
}
\usage{
panelNNET.default(  y, X, hidden_units, fe_var
  , maxit = 100, lam = 0, time_var = NULL, param = NULL
  , parapen = rep(0, ncol(param)), parlist = NULL, verbose = FALSE
  , para_plot = FALSE, report_interval = 100
  , save_each_iter = FALSE, path = NULL, tag = "", gravity = 1.01
  , convtol = 1e-8, bias_hlayers = TRUE, RMSprop = TRUE, start_LR = .01
  , activation = 'tanh', inference = TRUE, doscale = TRUE
  , treatment = NULL, interact_treatment = TRUE, batchsize = nrow(X)
  , maxstopcounter = 10, OLStrick = FALSE, useOptim = FALSE
  , optimMethod = 'BFGS', initialization = 'HZRS', ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
The response data
}
  \item{X}{
Covariates to enter non-parametrically
}
  \item{hidden_units}{
Numeric vector of hidden units within hidden layers.  First entry is the lowest layer, subsequent entries are higher layers.
}
  \item{fe_var}{
A factor indicating the cross-sectional unit.  At present only one cross-sectional unit is supported.
}
  \item{maxit}{
Maximum number of iterations
}
  \item{lam}{
Lambda, the penalty (or "weight decay") factor 
}
  \item{time_var}{
Numeric vector of the time variable
}
  \item{param}{
Terms to enter parametrically
}
  \item{parapen}{
Numeric vector multiplying the penalty terms on the the penalties for the parametric terms.  Defaults to a vector of zeros; parametric terms are unpenalized.  This is useful when inference on these terms is desired.  
}
  \item{parlist}{
A list of starting values for the parameters.  Chosen randomly if omitted (see the "initialization" argument).  Useful when re-starting an interrupted gradient descent. 
}
  \item{verbose}{
If true will print progress to console, and make plots of the algorithm's progress.  May work poorly on windows and/or in RStudio.
}
  \item{gravity}{
The learning rate will be multiplied by this factor after each step.
}
  \item{convtol}{
Convergence tolerance.  When <<maxstopcounter>> successive iterations fail to improve MSE by this amount, gradient descent exits
}
  \item{bias_hlayers}{
Add bias (i.e. intercept) terms to hidden layers (not top layer)
}
  \item{RMSprop}{
Gradient descent by RMSprop.  If FALSE, by "vanilla" gradient descent
}
  \item{start_LR}{
The initial learning rate a.k.a. step size
}
  \item{activation}{
"tanh", "logistic", "relu", or "lrelu" (for "leaky ReLU")
}
  \item{doscale}{
Should the input data be scaled?  Generally it should.
}
  \item{\dots}{
Not used.
}
  \item{treatment}{
For heterogeneous treatment effects models.  A binary vector of treatment assignments.  NOTE: HTE estimation isn't yet working smoothly.
}
  \item{interact_treatment}{
Logical.  Should the treatment be interacted with the upper-level derived variables?  NOTE: HTE estimation isn't yet working smoothly.
}
  \item{batchsize}{
Size of batches for minibatch gradient descent.  Defaults to nrow(X), which is batch gradient descent.
}
  \item{maxstopcounter}{
How many times should the learning rate be halved after an epoch increases MSE, before panelNNET exits?
}
  \item{OLStrick}{
Find the closed-form solution on the top layer of the network that minimizes the penalized loss function -- augmenting gradient descent.  
}
  \item{useOptim}{
Logical.  Use R's optim function for gradient descent?  This doesn't work very well with the OLS trick.
}
  \item{optimMethod}{
If useOptime == TRUE, an optimizer method.
}
  \item{initialization}{
If "HRZS", the weight initialization scheme proposed by He, Zhang, Ren, and Jian (2015).  If "XG", the weight initialization scheme proposed by Glorot and Bengio (2010).  If "enforce_normalization", uniform random numbers are scaled such that their variance equals one.  Otherwise, draws from a uniform distribution with bounds of -.7 and .7, following recommndations in Hastie, Tibshirani, and Friedman.
}
}
\details{
Function fits a model of the form

y_it = u_i + P_it B + Z_it C + e_it

Z_it = activation(X_it D)

in the single-layer case, and generalized in the multi-layer cased.  Parameters [B, C, and D] are fit by one of two methods of gradient descent, subject regularization by the parameter lam.  Given that the top-layer is linear, estimation is facilitated by the "within" transformation -- subtraction of the group-wise mean, in order to eliminate the fixed effects u_i.

}
\value{
  \item{yhat }{The fitted values}
  \item{weights }{The estimated parameters}
  \item{fe }{Estimates of the fixed effects, for each observation}
  \item{converged }{TRUE or FALSE}
  \item{mse }{mean squared error (in-sample)}
  \item{lam }{The supplied penalty}
  \item{hidden_layers }{The network architecture}
  \item{time_var }{The time variable supplied}
  \item{X }{The data supplied to enter non-parametrically}
  \item{y }{The supplied outcome}
  \item{param }{The data supplied to enter linearly}
  \item{fe_var }{The supplied cross-sectional unit}
  \item{hidden units }{The final estimate of the hidden units}
  \item{final_improvement }{The last improvement to MSE at exit}
  \item{msevec }{The evolution of MSE over the iterations}
  \item{RMSprop }{Whether RMSprop was used}
  \item{convtol }{The convergence tolerance used}
  \item{call }{The call}
  \item{grads }{The gradients at exit}
  \item{activation }{The activation function used}
  \item{parapen }{The factor that multiplies lambda for the penalty terms}
  \item{doscale }{Whether the data was scaled.  Note that the scaling factors can be gotten from the attributes of the returned X and param matrices}
  \item{treatment}{The user-supplied treatment assignments}
  \item{interact_treatment}{Logical, user-supplied}
  \item{batchsize}{User-supplied}
  \item{usedOptim}{User-supplied}
  \item{optimMethod}{User-supplied}
  \item{initialization}{User-supplied}
  \item{edf}{Effective degrees of freedom.  Only returned if inference set to true.}
  \item{vcs}{Various sorts of parameter covariance matrix.  Only returned if the function "do_inference" is run.  For example, output <- do_inference(output)}}
  \item{J}{The Jacobian.  Only returned if inference set to true.}
}

\references{
Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. The elements of statistical learning. Vol. 1. Springer, Berlin: Springer series in statistics, 2001.
}
\author{
Andrew Crane-Droesch
}
\note{
This package is in active development.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  set.seed(1)
  #Fake dataset
  N <- 1000
  p <- 20
  X <- as.data.frame(mvrnorm(N, rep(0, p), diag(rep(1, p))))
  id <- factor(0:(N-1)\%\%20+1)
  id.eff <- rnorm(nlevels(id), sd = 5)
  time <- 0:(N - 1)\%/\%20+1
  u <- rnorm(N, sd = 5)
  y <- sin(3*X$V1) - cos(4*X$V2) + 3*tanh((-2*X$V1+X$V2+X$V4)*X$V3) + X$V6/(X$V7+8) + id.eff[id] +
     .5*time - .005*time^2 + u
  hist(y)


  #Parametric and nonparametric terms
  X <- X
  P <- cbind(time, time^2)

  #Traiing and test set
  tr <- time<35
  te <- tr == FALSE

  #Fitting a two-layer neural net with 5 and 3 hidden units
  pnn <- panelNNET(y[tr], X[tr,], hidden_units = c(5,3)
    , fe_var = id[tr], lam = 1
    , time_var = time[tr], param = P[tr,], verbose = FALSE
    , bias_hlayers = TRUE, gravity = 1.01
    , RMSprop = TRUE, convtol = 1e-5, maxit = 10000
    , activation = 'tanh', doscale = TRUE, parapen = c(0,0)
  )

  plot(pnn)
  summary(pnn) #Approx inference


%  ##Predicting for the test set
%  pr <- predict(pnn, newX = X[te,], fe.newX = id[te], new.param = P[te,])
%  plot(y[te], pr)
%  mean((y[te] - pr)^2)
%  abline(0,1)


%  #Comparing to a vanilla neural net
%  library(nnet)
%  nn <- nnet(y = y[tr], x = cbind(id, X, P)[tr,], size = 20, linout = TRUE, maxit = 10000)
%  pr <- predict(nn, newdata = cbind(id, X, P)[te,])
%  plot(y[te], pr)
%  abline(0,1)
%  mean((y[te] - pr)^2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
