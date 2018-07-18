require(rgl)
require(mvtnorm)
require(RColorBrewer)
require(mgcv)
require(scam)

Directions <- rbind(
  c(-1,0,0),
  c(1,0,0),
  c(0,-1,0),
  c(0,1,0),
  c(0,0,-1),
  c(0,0,1))




MakePlot <- function(WinRange,PlotType,zmax = WinRange){
	Ticks <- seq(0,WinRange,by=0.25)
	TickLen <- 0.05
	Sign <- 2*((1:6-1)%%2)-1
	
	AllTicks <- vector("list",6)

	for (i in 1:6){
		for (j in 2:(length(Ticks)-1)){
			for (k in 1:6){
				AllTicks[[i]] <- rbind(AllTicks[[i]],
															Ticks[j]*Directions[i,],Ticks[j]*Directions[i,]+TickLen*Directions[k,])
			}
		}
	}

  plot3d(0,0,0, xlab="", ylab="", zlab="", xlim=c(-WinRange,WinRange), ylim=c(-WinRange,WinRange), zlim=c(-WinRange,WinRange),axes=FALSE)
  AxesDraw <- 1:6
	if (PlotType == 3) AxesDraw <- c(1:4,6)
	for (i in AxesDraw){
    arrow3d(c(0,0,0), WinRange*Directions[i,], type = "lines",     col = "black",s=.1)
    segments3d(AllTicks[[i]])
    for (j in 2:(length(Ticks)-1)){
      if (Ticks[j]%%1 == 0 & i<5){
        text3d(Ticks[j]*Directions[i,]+0.2*Directions[6,], texts=Sign[i]*Ticks[j])
      } else if (Ticks[j]%%1 == 0){
        text3d(Ticks[j]*Directions[i,]+0.2*Directions[2,], texts=Sign[i]*Ticks[j]*zmax/WinRange)
      }
    }
  }
	if (PlotType == 1){
		Labels <- c("PCA 1","PCA 2","PCA 3")
	} else if (PlotType == 2){
		Labels <- c("H/A z","W/A z","W/H z")
	} else {
		Labels <- c("PCA 1","PCA 2","U5M rate")
	}
  
  for (i in 1:3){
    text3d(WinRange*Directions[2*i,]+.3*Directions[6,], texts = Labels[i],cex=1.2)
  }
}

