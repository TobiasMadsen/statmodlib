% Generated by roxygen2 (4.1.0.9001): do not edit by hand
% Please edit documentation in R/statmodlib.r
\name{estimates2dSSD}
\alias{estimates2dSSD}
\title{Estimer 2D normal fordeling}
\usage{
estimates2dSSD(S_x, S_y, SSD_x, SSD_y, SPD_xy, n)
}
\arguments{
\item{S_x}{summen af x'er}

\item{S_y}{summen af y'er}

\item{SSD_x}{summen af (x-mean(x))^2}

\item{SSD_y}{summen af (y-mean(y))^2}

\item{SPD_xy}{summen af (x-mean(x))*(y-mean(y))}

\item{n}{stikprøvestørrelsen}
}
\description{
(x,y) er en observationsrække fra den 2dimensionale
normalfordeling.
}

