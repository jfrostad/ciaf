if (!require("pacman")) install.packages("pacman")
pacman::p_load(RColorBrewer)


DropboxFolder <- "IHME/CGF_CBH"
if (Sys.info()["nodename"] == "DT-30200298"){
  Prefix <- "E:/Dropbox/"
} else {
  Prefix <- "C:/Users/bcreiner/Dropbox/"
}

folder <- paste0(Prefix,DropboxFolder)

setwd(folder)

cgf <- read.csv("data/ciaf_collapsed_2017.12.27.csv")
u5m <- read.csv("data/u5m_for_bobby.csv")
cgf$admin_1 <- as.character(cgf$admin_1)
u5m$admin_1 <- as.character(u5m$admin_1)




head(u5m)
head(cgf)

NData <- length(cgf[,1])

year <- cgf$year
nid <- cgf$nid
admin_1 <- cgf$admin_1
meanHAZ <- cgf$avg_HAZ
meanWAZ <- cgf$avg_WAZ
meanWHZ <- cgf$avg_WHZ
cgfN <- cgf$total_N
u5mN <- u5mD <- u5mq <- rep(NA, NData)
error <- {}
good <- {}
for (i in 1:NData){
	tmploc <- which(u5m$nid == nid[i] & u5m$admin_1 == admin_1[i])
	if (length(tmploc) != 1){
		error <- c(error,i)
	} else {
		good <- c(good,i)
		u5mN[i] <- u5m$childyrs[tmploc]
		u5mq[i] <- u5m$q[tmploc]
		u5mD[i] <- u5m$died[tmploc]
	}
}

data <- data.frame(nid, year, admin_1, meanHAZ, meanWAZ, meanWHZ, cgfN, u5mN, u5mD, u5mq)
data <- data[-error,]

baddata <- which(data$meanHAZ == 0)
data <- data[-baddata,]

