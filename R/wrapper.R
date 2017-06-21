wbgt <- function (year, month, day, hour, minute,
                  gmt, avg, lat, lon, solar, pres, 
                  Tair, relhum, speed, zspeed, dT, urban)
{
  num_obs <- length(year)
  Tg=rep(0.0, num_obs)
  Tnwb=rep(0.0, num_obs)
  Tpsy=rep(0.0, num_obs)
  Twbg=rep(0.0, num_obs)
  status=rep(0, num_obs)
  
out <- .C("wbgt", num_obs=as.integer(num_obs), year=as.integer(year), month=as.integer(month), day=as.integer(day),
            hour=as.integer(hour), minute=as.integer(minute), gmt=as.integer(gmt), avg=as.integer(avg),
            lat=as.double(lat), lon=as.double(lon), solar=as.double(solar), pres=as.double(pres), 
            Tair=as.double(Tair), relhum=as.double(relhum), speed=as.double(speed), 
            zspeed=as.double(zspeed), 
            dT=as.double(dT), urban=as.integer(urban),
            Tg=as.double(Tg), 
            Tnwb=as.double(Tnwb), Tpsy=as.double(Tpsy), Twbg=as.double(Twbg), status=as.integer(status), PACKAGE="wbgt")

# any calculation based on pressure outside of normal range yields NA
out[pres > 2000 | pres < 800] <- NA
out
}

# Apply the wbgt function above to a data frame, assumed to contain columns with names
# equal to the names of the inputs to wbgt
# Returns copy of the original data frame with the wbgt appended as a new column 
wbgt_df <- function(data) {
  if (!requireNamespace("dplyr", quietly=TRUE)) {
    stop("Please install dplyr in order to use this function.")
  }
  with(data, dplyr::mutate(data, wbgt=wbgt(year, month, day, hour, minute,
             gmt, avg, lat, lon, solar, pres, 
             Tair, relhum, speed, zspeed, dT, urban)$Twbg))
}
