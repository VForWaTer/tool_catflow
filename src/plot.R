###
#
plot_profile_geo <- function()
  {
plot(hill$short_rep_hill$short_dist, hill$short_rep_hill$short_elev, col=as.factor(hill$short_rep_hill$geo))
  }


###
# plot representiv hillslope 
plot_rep_hill <- function(hillslope_dist, hillslope_width, rect_col=NA){


  delta_d <- (diff(hillslope_width, lag=1)*-1)/2
  plot(1, type="n", axes=T, xlab="", ylab="", xlim=c(max(hillslope_width)*-0.5,
                                                     max(hillslope_width)),
       ylim=c(0, max(hillslope_dist)), bty='n')

  rect(delta_d[1], 0, hillslope_width[1] - delta_d[1], hillslope_dist[1], col=rect_col)

  for(x in 2:length(delta_d))
  { 
    rect(sum(delta_d[1:x-1]) + delta_d[x], hillslope_dist[x-1],
         hillslope_width[1] - (sum(delta_d[1:x-1]) + delta_d[x]),
         hillslope_dist[x], col=rect_col)  
  }
}


###
#
plot_2d_profile <- function(rep_hill, sp_elev, short_dist, short_elev, hillslope_nr)
  {
  plot(rep_hill$mean_dist, rep_hill$mean_elev, xlab='distance to river [m]', ylab='elevation above river [m]', main=hillslope_nr)
  lines(sp_elev, lwd=2, col='blue')
  points(short_dist, short_elev, pch=20, col='green')
  abline(v=max(short_dist), col='red', lwd=2)
  }


###
#
plot_hillslope_width <- function(rep_hill, short_dist, short_width, short_width_area_corr, sm_spline)
  {
  plot(rep_hill$mean_dist, rep_hill$width, pch=1, , xlab='distance to river [m]', ylab='hillslope width [m]')
  lines(sm_spline, lwd=2, col='blue')
  abline(h=quantile(rep_hill$width, 0.1), lwd=2, col='black', pch=2)
  points(short_dist, short_width, pch=20, col='darkgreen')
  points(short_dist, short_width_area_corr, pch=5, col='red')
  }