\name{evolvability_all}
\alias{evolvability_all}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to estimate Evolvability Conditoinal evolvability, autonomy and Integration
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
evolvability_all(name, nb_pop, names_pop, cpus = 1, simu = TRUE, analytic = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{
prefix of the Gmatrices
}
  \item{names_pop}{
names of the populations used for the comparison
}
  \item{cpus}{
nb_cpus used, let 1 for analytic results
}
  \item{simu}{
%%     ~~Describe \code{simu} here~~
}
  \item{analytic}{
%%     ~~Describe \code{analytic} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

###################  EXEMPLE for 3 populations

data(animal_results)

Gmat_1=ani_1$VCV #extract variance part of the model to the first MCMCglmm result
Gmat_1=Gmat_1[,grep(".animal",colnames(Gmat_1))] # extract the genetic variances and covariances
Gmat_2=ani_2$VCV 
Gmat_2=Gmat_2[,grep(".animal",colnames(Gmat_1))]
Gmat_3=ani_3$VCV #extract variance part of the model
Gmat_3=Gmat_3[,grep(".animal",colnames(Gmat_1))]


evolvability_all(name = "Gmat_",names_pop = c("pop1","pop2","pop3"))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
