###############################################################################
# R functions for analysing European gypsy moth, L. dispar, 
# range expansion in heterogenous environments 
###############################################################################
# Audrey Lustig
# Bio-Protection Research Centre, Lincoln University, NZ
# March 2016

require(ggplot2)
require(reshape2)
require(gridExtra)
library(plotly)


path_to_mdig_casetudy_folder <- 'path_to_mdig-casetudy_folder'

working_directory <- paste(path_to_mdig_casetudy_folder, 'L.Dispar/Results/analysis/',sep='') # working directory

###############################################################################
# Read Gypsy moth density and rate of spread data for the following parameter

alleeN=2 # alle threshold
introN=0 # number of introduced individual
REP=seq(1,10) # landscape replication
H=c(0.3,0.5,0.7) # spatial autocorrelation, H
P=c(35,55,75) # Percentage of suitable habitat cover, P
R=c(0.815,1.223) # intrinsic rate of increase
DIST=c(3,5) # dispersal distance

data.den <- NULL # density data frame
data.ROS <- NULL # rate of spread data frame
data.OA <- NULL
# Loop over the data set
for (h in H){
for (p in P){
for (rep in REP){
for (dist in DIST){
for (r in R){

dt <- read.table(paste(working_directory,'Frac.E.1280.R.10.H.', h ,'.P.',p,'.rawASCIIgrid-rep',rep,'.asc/allee.',alleeN,'.intro.',introN,'.nbp.5.dist.',dist,'.repro.',r,'.K.50.prob_disp.0.05.gz',sep=''))

data.OA <- rbind(data.OA,c(h, p, r, dist, rep, dt[,2]))
data.den <- rbind(data.den,c(h, p, r, dist, rep, dt[,1]))
data.ROS <- rbind(data.ROS,c(h, p, r, dist, rep, dt[,4]))
}}}}}


###############################################################################
# Plot density as a function of time
data.den<-as.data.frame(data.den) 
mymelt <- NULL
for (line in seq(1, nrow(data.den))){
	tts <- data.den[line,7:ncol(data.den)]/(16384*data.OA[line,7:ncol(data.OA)]) 
	r = data.den[line,3]
	dist = data.den[line,4]
	if(r==0.815 & dist ==3){c <- 1}
	if(r==0.815 & dist ==5){c <- 2}
	if(r==1.223 & dist ==3){c <- 3}
	if(r==1.223 & dist ==5){c <- 4}
	mymelt <- rbind(mymelt,data.frame(t=seq(1,(length(tts))),density=as.numeric(tts),ts=rep(paste('t',line,sep=''),length(tts)),r=rep(data.den[line,3],length(tts)),dist=rep(data.den[line,4],length(tts)),parameters=rep(c,length(tts))))}
mymelt$parameters<-as.factor(mymelt$parameters)
pTime_density<-ggplot(mymelt, aes(x=t, y=density)) + geom_line(aes(fill=ts,colour=parameters))+xlab(' ')+ylab('Population density (d)')+ guides(color=guide_legend(title="Parameters"))+scale_color_manual(labels = c("r=0.815, dist=3", "r=0.815, dist=5","r=1.223, dist=3","r=1.223, dist=5"), values = c(1,2,3,4)) +theme_bw() + theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), legend.position=c(0.3,0.7), legend.text = element_text(size = 11),legend.title = element_text(size=11, face="bold"), axis.title.y = element_text(size = 12),axis.title.x = element_text(size = 12)) + annotate("text", x = min(mymelt$t), y = max(mymelt$density),label="(A)")
pTime_density


# Plot ROS as a function of time
data.ROS<-as.data.frame(data.ROS)
mymelt <- NULL
for (line in seq(1, nrow(data.ROS))){
	tts <- data.ROS[line,7:ncol(data.ROS)]#ncol(data.ROS)]
	r = data.ROS[line,3]
	dist = data.ROS[line,4]
	if(r==0.815 & dist ==3){c <- 1}
	if(r==0.815 & dist ==5){c <- 2}
	if(r==1.223 & dist ==3){c <- 3}
	if(r==1.223 & dist ==5){c <- 4}
	mymelt <- rbind(mymelt,data.frame(t=seq(1,(length(tts))),ROS=as.numeric(tts),ts=rep(paste('t',line,sep=''),length(tts)),r=rep(data.ROS[line,3],length(tts)),dist=rep(data.den[line,4],length(tts)),parameters=rep(c,length(tts))))
}
mymelt$parameters<-as.factor(mymelt$parameters)
pTime_ROS<-ggplot(mymelt, aes(x=t, y=ROS)) + geom_line(aes(fill=ts,colour=parameters))+xlab('Time (years)')+ylab('Rate of spread (ROS)')+ guides(color=guide_legend(title="Parameters"))+scale_color_manual(labels = c("r=0.815, dist=3", "r=0.815, dist=5","r=1.223, dist=3","r=1.223, dist=5"), values = c(1,2,3,4)) +theme_bw() + theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), legend.position=c(0.3,0.8), legend.text = element_text(size = 12),legend.title = element_text(size=14, face="bold"), axis.title.y = element_text(size = 12),axis.title.x = element_text(size = 12))+ theme(legend.position="none") + annotate("text", x = min(mymelt$t), y = max(mymelt$ROS),label="(B)")
pTime_ROS



###############################################################################
# Plot density as a function of PLAND
dat.denb <- NULL
for (line in seq(1, nrow(data.den))){
	tts <- data.den[line,7:ncol(data.den)]/(16384*data.OA[line,7:ncol(data.OA)]) 
	dat.denb <- c(dat.denb,mean(as.numeric(tts)))
} 
dat.denb <- cbind(data.den[,1:5], dat.denb)
dat.denb <-as.data.frame(dat.denb)
names(dat.denb) <- c('H','P','r','dist','rep','d')

df.m <- melt(dat.denb[,c(2,3,4,6)], id.var = c('P','r','dist'))
Parameters <- rep(1,nrow(df.m))
for (j in seq(1,nrow(df.m))){
r <- df.m[j,2]
dist <- df.m[j,3]
if(r==0.815 & dist ==3)
Parameters[j] = 1
if(r==0.815 & dist ==5)
Parameters[j] = 2
if(r==1.223 & dist ==3)
Parameters[j]=  3	
if(r==1.223 & dist ==5)
Parameters[j]= 4
} 
df.t<-cbind(df.m[,c(1,2,3,5)],Parameters)
df.t$P <- as.factor(df.t$P)
df.t$Parameters <- as.factor(df.t$Parameters)
pPLAND_density<-ggplot(data = df.t, aes(x=P, y=value)) + geom_boxplot(aes(fill=Parameters)) + xlab(" ") + ylab(" ") + guides(color=guide_legend(title="Parameters"))+scale_fill_manual(labels = c("r=0.815, dist=3", "r=0.815, dist=5","r=1.223, dist=3","r=1.223, dist=5"), values = c(1,2,3,4))+theme(legend.position=c(0.1,0.9))+theme_bw() + theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), legend.position=c(0.3,0.8), legend.text = element_text(size = 12),legend.title = element_text(size=14, face="bold"), axis.title.y = element_text(size = 12),axis.title.x = element_text(size = 12))+ theme(legend.position="none") + annotate("text", x = 0.6, y = max(df.t$value),label="(C)")
pPLAND_density

# Plot ROS as a function of PLAND
dat.ROSb <- cbind(data.ROS[,1:5], apply(data.ROS, 1, function(x) mean(x[7:ncol(data.ROS)]))) # Mean rate of spread over landscape replication
dat.ROSb <-as.data.frame(dat.ROSb)
names(dat.ROSb) <- c('H','P','r','dist','rep','d')

df.m <- melt(dat.ROSb[,c(2,3,4,6)], id.var = c('P','r','dist'))
Parameters <- rep(1,nrow(df.m))
for (j in seq(1,nrow(df.m))){
r <- df.m[j,2]
dist <- df.m[j,3]
if(r==0.815 & dist ==3)
Parameters[j] = 1
if(r==0.815 & dist ==5)
Parameters[j] = 2
if(r==1.223 & dist ==3)
Parameters[j]=  3	
if(r==1.223 & dist ==5)
Parameters[j]= 4
} 
df.t<-cbind(df.m[,c(1,2,3,5)],Parameters)
df.t$P <- as.factor(df.t$P)
df.t$Parameters <- as.factor(df.t$Parameters)
pPLAND_ROS<-ggplot(data = df.t, aes(x=P, y=value)) + geom_boxplot(aes(fill=Parameters)) + xlab("Percentage of suitable habitat cover (PLAND)") + ylab(" ") + guides(color=guide_legend(title="Parameters"))+scale_fill_manual(labels = c("r=0.815, dist=3", "r=0.815, dist=5","r=1.223, dist=3","r=1.223, dist=5"), values = c(1,2,3,4))+theme_bw() + theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), legend.position=c(0.3,0.8), legend.text = element_text(size = 12),legend.title = element_text(size=14, face="bold"), axis.title.y = element_text(size = 12),axis.title.x = element_text(size = 12))+ theme(legend.position="none") + annotate("text", x = 0.6, y = max(df.t$value),label="(D)")
pPLAND_ROS




###############################################################################
### Load Fragstat landscape metrics and select connectivity index
load(paste(path_to_mdig_casetudy_folder, 'L.Dispar/Data/LandscapeMetrics/Landscape_metrics.RData',sep=''))
metrics.names <- listdata[1][[1]] # name of Fragstat landscape metrics
metrics.type <- listdata[2][[1]] # metrics class: in this case the metrics are always calculated for the suitable area
landscape.name <- listdata[3][[1]] # name of the ASCII file
landscape.param <- listdata[4][[1]] # Parameter of the landscape (extent, resolution, H, P, replication, class metrics, and landscape name)
metrics.value <- listdata[5][[1]] # Metrics value

CONNECT <- metrics.value[,metrics.names=='CONNECT'] # Select metrics CONNECT
connect.mat <- cbind(landscape.param[,-c(6,7)],CONNECT) # combine landscape parameter and metrics value
connect.mat <- connect.mat[order(connect.mat$H,connect.mat$P,connect.mat$rep),] # Order the table as a function of H, P and replication


### Plot density as a function of CONNECT
dat.denb <- NULL
for (line in seq(1, nrow(data.den))){
	tts <- data.den[line,7:ncol(data.den)]/(16384*data.OA[line,7:ncol(data.OA)]) 
	dat.denb <- c(dat.denb,mean(as.numeric(tts)))
} 
dat.denb <- cbind(data.den[,1:5], dat.denb)
dat.denb <-as.data.frame(dat.denb)
names(dat.denb) <- c('H','P','r','dist','rep','d')
dat.denb <- dat.denb[order(dat.denb$H,dat.denb$P,dat.denb$rep),]
connect.index<-as.numeric(apply(connect.mat, 1, function(x) rep(x[6],4)))
db <- cbind(connect.index,dat.denb)
db$Parameters <- rep(1,nrow(db))
for (line in seq(1, nrow(db))){
	r = db[line,4]
	dist = db[line,5]
	if(r==0.815 & dist ==3){db$Parameters[line] = 1}
	if(r==0.815 & dist ==5){db$Parameters[line] = 2}
	if(r==1.223 & dist ==3){db$Parameters[line] = 3}
	if(r==1.223 & dist ==5){db$Parameters[line] = 4}
}
db$connect.index <- as.numeric(db$connect.index)
db$Parameters <- as.factor(db$Parameters)
pCONNECT_density<-ggplot(db) + geom_jitter(aes(y=d,x=connect.index, colour=Parameters),) + geom_smooth(aes(y=d,x=connect.index, colour=Parameters), method=lm, se=FALSE) + xlab(" ") + ylab(" ") + guides(color=guide_legend(title="Parameters"))+scale_colour_manual(labels = c("r=0.815, dist=3", "r=0.815, dist=5","r=1.223, dist=3","r=1.223, dist=5"), values = c(1,2,3,4))+theme(legend.position=c(0.9,0.9))+ theme(legend.text = element_text(size = 12))+theme_bw() + theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), legend.position=c(0.3,0.8), legend.text = element_text(size = 12),legend.title = element_text(size=14, face="bold"), axis.title.y = element_text(size = 12),axis.title.x = element_text(size = 12))+ theme(legend.position="none")+ annotate("text", x = min(db$connect.index), y = max(db$d),label="(E)")
pCONNECT_density

### Plot ROS as a function of CONNECT
dat.ROSb <- cbind(data.ROS[,1:5], apply(data.ROS, 1, function(x) mean(x[7:ncol(data.ROS)])))
dat.dROSb <-as.data.frame(dat.ROSb)
names(dat.ROSb) <- c('H','P','r','dist','rep','ROS')
dat.ROSb <- dat.ROSb[order(dat.ROSb$H,dat.ROSb$P,dat.ROSb$rep),]
connect.index<-as.numeric(apply(connect.mat, 1, function(x) rep(x[6],4)))
db <- cbind(connect.index,dat.ROSb)
db$Parameters <- rep(1,nrow(db))
for (line in seq(1, nrow(db))){
	r = db[line,4]
	dist = db[line,5]
	if(r==0.815 & dist ==3){db$Parameters[line] = 1}
	if(r==0.815 & dist ==5){db$Parameters[line] = 2}
	if(r==1.223 & dist ==3){db$Parameters[line] = 3}
	if(r==1.223 & dist ==5){db$Parameters[line] = 4}
}
db$connect.index <- as.numeric(db$connect.index)
db$Parameters <- as.factor(db$Parameters)
pCONNECT_ROS<-ggplot(db) + geom_jitter(aes(y=ROS,x=connect.index, colour=Parameters),) + geom_smooth(aes(y=ROS,x=connect.index, colour=Parameters), method=lm, se=FALSE) + xlab("Connectance (CONNECT)") + ylab(" ") + guides(color=guide_legend(title="Parameters"))+scale_colour_manual(labels = c("r=0.815, dist=3", "r=0.815, dist=5","r=1.223, dist=3","r=1.223, dist=5"), values = c(1,2,3,4))+theme(legend.position=c(0.9,0.9))+theme_bw() + theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), legend.position=c(0.3,0.8), legend.text = element_text(size = 12),legend.title = element_text(size=14, face="bold"), axis.title.y = element_text(size = 12),axis.title.x = element_text(size = 12))+ theme(legend.position="none")+ annotate("text", x = min(db$connect.index), y = max(db$ROS),label="(F)")
pCONNECT_ROS




#####################################################################################################################
### Figure 4
grid.arrange(pTime_density, pPLAND_density, pCONNECT_density, pTime_ROS, pPLAND_ROS, pCONNECT_ROS, nrow = 2, ncol = 3)




