legend.col <- function(col, lev){
 	opar <- par
 	n <- length(col)
 	
 	
 	bx <- par("usr")
 	box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
	bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
	box.cy <- c(bx[3], bx[3])
	box.sy <- (bx[4] - bx[3]) / n

	xwd<-box.cx[2]-box.cx[1]
 	box.cx[2]<-box.cx[2] - 0.5*xwd

	ywd<-box.cy[2]-box.cy[1]
 	box.cy[2]<-box.cy[2] - 0.5*xwd
## 	box.cy[1]<-box.cy[1] + 0.25*xwd
	box.sy<-box.sy*0.5
 	
	xx <- rep(box.cx, each = 2)
 
	par(xpd = TRUE)

	text(box.cx[1]-xwd*0.1,1+box.cy[1],"0",adj=c(1,0),cex=1.5)
	for(i in 1:n){
		yy <- c(box.cy[1] + (box.sy * (i - 1)),
		box.cy[1] + (box.sy * (i)),
		box.cy[1] + (box.sy * (i)),
		box.cy[1] + (box.sy * (i - 1)))
		polygon(xx, yy+1, col = col[i], border = col[i])
		if(i%%25==0){
			text(box.cx[1]-xwd*0.1,1+box.sy*(i)+box.cy[1],i/100,adj=c(1,1),cex=1.5)
		}
	}
	text(box.cx[1],1+box.sy*102+box.cy[1],expression(paste("D","'")),adj=c(0,0),cex=1.5)

	par(new = TRUE)
	plot(0, 0, type = "n",
	ylim = c(min(lev), max(lev)),
	yaxt = "n", ylab = "",
	xaxt = "n", xlab = "",
	frame.plot = FALSE)

	par <- opar
}

