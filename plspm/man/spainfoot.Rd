\name{spainfoot}
\alias{spainfoot}
\docType{data}
\title{Spanish football dataset}
\description{
  This data set contains the results of the teams in the Spanish football league 2008-2009.
}
\usage{spainfoot}
\format{
  A data frame with 20 observations on 14 variables. The variables may be used to construct four latent concepts: 1) \code{ATTACK}=Attack, 2) \code{DEFENSE}=Defense, 3) \code{SUCCESS}=Success, 4) \code{INDIS}=Indiscipline. 

  \code{ATTACK}: Recommended variables GSH, GSA, SSH, SSA. \cr
  \code{DEFENSE}: Recommended variables GCH, GCA, CSH, CSA. \cr
  \code{SUCCESS}: Recommended variables WMH, WMA, LWR, LRWL. \cr
  \code{INDIS}:  Recommended variables YC and RC. \cr
  \describe{
    \item{\code{GSH}}{Goals Scored Home: total number of goals scored at home}
    \item{\code{GSA}}{Goals Scored Away: total number of goals scored away}
    \item{\code{SSH}}{Success to Score Home: Percentage of matches with scores goals at home}
    \item{\code{SSA}}{Success to Score Away: Percentage of matches with scores goals away}
    \item{\code{GCH}}{Goals Conceded Home: total number of goals conceded at home}
    \item{\code{GCA}}{Goals Conceded Away: total number of goals conceded away}
    \item{\code{CSH}}{Clean Sheets Home: percentage of matches with no conceded goals at home}
    \item{\code{CSA}}{Clean Sheets Away: percentage of matches with no conceded goals away}
    \item{\code{WMH}}{Won Matches Home: total number of won matches at home }
    \item{\code{WMA}}{Won matches Away: total number of won matches away}
    \item{\code{LWR}}{Longest Winning Run: longest run of won matches}
    \item{\code{LRWL}}{Longest Run Without Loss: longest run of matches without losing}
    \item{\code{YC}}{Yellow Cards: total number of yellow cards}
    \item{\code{RC}}{Red Cards: total number of red cards}
  }
}
\source{
  League Day. \url{http://www.leagueday.com} \cr
  Statto. \url{http://www.statto.com} \cr
  BDFutbol. \url{http://www.statto.com} \cr
  Cero a cero. \url{http://www.ceroacero.es/} \cr
}
\references{
  Applying plspm Beginners Guide. \url{http://www.scribd.com/plspm}
}
\examples{
  data(spainfoot)
  spainfoot
}
\keyword{datasets}
