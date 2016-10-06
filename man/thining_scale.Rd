\name{thining_scale}
\alias{thining_scale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to merge results from multiple MCMCglmmm outputs, perform scaling and increase thining intarval
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
thining_scale(name = NULL, nb_gmat = NULL, scale = NULL, traits = NULL, div_thin = NULL, mean_std = FALSE, var_std = FALSE, data_ext = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{name}{
name corresonping toi the first part of the Rdata
}
  \item{nb_gmat}{
number of file ti merge
}
  \item{scale}{
set scale for each traits (ie c(1,1,10) if 3 traits
}
  \item{traits}{
names of traits which have to be kept
}
  \item{div_thin}{
number to increase the thining interval (2 = intial thining*2)
}
  \item{mean_std}{
%%     ~~Describe \code{mean_std} here~~
}
  \item{var_std}{
%%     ~~Describe \code{var_std} here~~
}
  \item{data_ext}{
data set to perform mean standardisation
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
