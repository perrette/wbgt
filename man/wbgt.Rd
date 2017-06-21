\name{wbgt}
\alias{wbgt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute the web bulb globe temperature (wbgt)}
\description{
    An R wrapper for Liljegren's C implemenation of his method for computing the wbgt.
}
\usage{
wbgt(year, month, day, hour, minute, gmt, avg, 
     lat, lon, solar, pres, Tair, relhum, speed, 
     zspeed, dT, urban)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{year}{
    An integer vector of four-digit years.
}
	\item{month}{
		An integer vector of months (1-12), or 0, indicating that `day` is the day of the year.
	}
	\item{day}{
		An integer vector of days, either of month of year depending upon value of `month`.
	}
	\item{hour}{
		An integer vector of hours in local standard time (LST).
	}
	\item{minute}{
		An integer vector of minutes.
	}
	\item{gmt}{
		An integer vector of LST-GMT in hours.
	}
	\item{avg}{
		An integer vector minutes, giving averaging time of meteorological inputs.
	}
	\item{lat}{
		A numeric vector of lattitude values.
	}
	\item{lon}{
		A numeric vector of longitude values.
	}
	\item{solar}{
		A numeric vector solar irradiance in W/m^2.
	}
	\item{pres}{
		A numeric vector barometric pressures in mb.
	}
	\item{Tair}{
		A numeric vector of air temperatures in degrees C.
	}
	\item{relhum}{
		A numeric vector of relative humidities in percent.
	}
	\item{speed}{
		A numeric vector of wind speeds in m/s.
	}
	\item{zspeed}{
		A numeric vector of heights at which wind speeds were measured in m.
	}
	\item{dT}{
		A numeric vector of temperature differential between `zspeed` and reference height.
	}
	\item{urban}{
		An integer vector indicating "urban" (1) or "rural" (0) measurement location.
	}

}
\value{
    A list containing all of the input variables (described above) in addition to
	`Twbg` (the wet bulb globe temperature in degrees C), `Tpsy` (psychrometric wet bulb temperature in degrees C), 
	`Tnwb` (natural wet bulb temperature in degrees C), `Tg` (the globe temperature in degrees C), and `status`
	(with 0 indicating convergence and -1 indicating non-convergence of the Liljegren algorithm).
}
\author{
    Max Lieblich, June Spector
}

\examples{
wbgt(2015, 1, 2, 4, 4, -8, 15, 46.2, -120.0, 1000, 1017.8, 12, 45, 1.24, 2.1, -0.1, 0)$Twbg
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{manip}
