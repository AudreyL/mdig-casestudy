###############################################################################
# R functions for analysing Great White Butterfly, P. Brassicae, 
# range expansion in heterogenous environments (Nelson area, New Zealand)
###############################################################################
# Dr. Audrey Lustig
# Bio-Protection Research Centre, Lincoln University, NZ
# May 2016

require(ggplot2)
require(gridExtra)

working_directory <- 'path_to_folder_mdig-casestudy'

###############################################################################
# Read P. Brassicae spread data 
dt_noE <- read.table(paste(working_directory,'P.Brassicae/Results/P.Brassicae_spread_noElevation_LS1/analysis/average_area.dat',sep='')) # load average number of sites occupied in  survival layer 1 (LS_1)
dt_withE <- read.table(paste(working_directory,'P.Brassicae/Results/P.Brassicae_spread_withElevation_LS2/analysis/average_area.dat',sep='')) # load average number of sites occupied in  survival layer 2 (LS_2)

# redifine time step
dt_noE[,1]<-seq(2010,2030,0.5)
dt_withE[,1] <- seq(2010,2030,0.5)
###############################################################################
# Plot average occupied area and confidence interval as a function of time

dt_nb_sites_noE<-cbind(dt_noE[,1],apply(dt_noE[,2:dim(dt_noE)[2]], 1, mean),apply(dt_noE[,2:dim(dt_noE)[2]], 1, sd),rep('LS_1',dim(dt_noE)[1])) # number of sites occupied as a function of time for LS1
dt_nb_sites_noE<-as.data.frame(dt_nb_sites_noE)
names(dt_nb_sites_noE)<-c('t','Nb.sites','CI','ID')
dt_nb_sites_noE

dt_nb_sites_withE<-cbind(dt_withE[,1],apply(dt_withE[,2:dim(dt_withE)[2]], 1, mean),apply(dt_withE[,2:dim(dt_withE)[2]], 1, sd),rep('LS_2',dim(dt_withE)[1])) # number of sites occupied as a function of time for LS2
dt_nb_sites_withE<-as.data.frame(dt_nb_sites_withE)
names(dt_nb_sites_withE)<-c('t','Nb.sites','CI','ID')
dt_nb_sites_withE

# transform data for plotting
mymelt<-rbind(dt_nb_sites_noE,dt_nb_sites_withE)
mymelt[,1]<-as.numeric(as.character(mymelt[,1]))
mymelt[,2]<-as.numeric(as.character(mymelt[,2]))*0.0025 #0.0025 is the raster cell size area in km2
mymelt[,3]<-as.numeric(as.character(mymelt[,3]))*0.0025

cbPalette <- c("#E69F00", "#0072B2") # define color
# Figure 2-A
pTime_area<-ggplot(mymelt, aes(x=t, y=Nb.sites, ymin=Nb.sites-CI,ymax=Nb.sites+CI, color=factor(ID)))  + geom_line() + geom_pointrange() +xlab('')+ylab('Invaded Area (km sq.)') + guides(color=guide_legend(title="Survival layer"))+theme_bw() + theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), legend.position=c(0.15,0.8), legend.text = element_text(size = 12),legend.title = element_text(size=14,face="bold"), axis.title.y = element_text(size = 12),axis.title.x = element_text(size = 12)) + scale_colour_manual(values=cbPalette) + annotate("text", x = 2010, y = 370,label="(A)")
pTime_area 


###############################################################################
# Plot rate of spread and confidence interval as a function of time
# Calculate ROS as the number of new cell occupied at each time step
# Plot average occupied area and confidence interval as a function of time
dt_noE=cbind(dt_noE[,1],dt_noE[,-1]/200712*100) # number of highly suitable sites occupied as a function of time for LS1
dt_nb_sites_noE<-cbind(dt_noE[,1],apply(dt_noE[,2:dim(dt_noE)[2]], 1, mean),apply(dt_noE[,2:dim(dt_noE)[2]], 1, sd),rep('LS_1',dim(dt_noE)[1]))
dt_nb_sites_noE<-as.data.frame(dt_nb_sites_noE)
names(dt_nb_sites_noE)<-c('t','Nb.sites','CI','ID')
dt_nb_sites_noE



dt_withE=cbind(dt_withE[,1],dt_withE[,-1]/52970.39*100) # number of highly suitable sites occupied as a function of time for LS1
dt_nb_sites_withE<-cbind(dt_withE[,1],apply(dt_withE[,2:dim(dt_withE)[2]], 1, mean),apply(dt_withE[,2:dim(dt_withE)[2]], 1, sd),rep('LS_2',dim(dt_withE)[1]))
dt_nb_sites_withE<-as.data.frame(dt_nb_sites_withE)
names(dt_nb_sites_withE)<-c('t','Nb.sites','CI','ID')
dt_nb_sites_withE

# transform data for plotting
mymelt<-rbind(dt_nb_sites_noE,dt_nb_sites_withE)
mymelt[,1]<-as.numeric(as.character(mymelt[,1]))
mymelt[,2]<-as.numeric(as.character(mymelt[,2]))
mymelt[,3]<-as.numeric(as.character(mymelt[,3]))

# Figure 2-B
pTime_lativeArea<-ggplot(mymelt, aes(x=t, y=Nb.sites, ymin=Nb.sites-CI,ymax=Nb.sites+CI, color=factor(ID)))  + geom_line() + geom_pointrange() +xlab('Time')+ylab('Proportion of highly suitable area invaded (%)') + guides(color=guide_legend(title="Survival layer"))+theme_bw() + theme(axis.line = element_line(colour = "black"), panel.background = element_blank(), legend.position='none', legend.text = element_text(size = 12),legend.title = element_text(size=14,face="bold"), axis.title.y = element_text(size = 12),axis.title.x = element_text(size = 12)) + scale_colour_manual(values=cbPalette) + ylim(0,100) + annotate("text", x = 2010, y = 100,label="(B)")
pTime_lativeArea 



###############################################################################
# Figure 3-C
dt_metrics <- read.table(paste(working_directory,'P.Brassicae/Data/LandscapeMetrics/Landscape_metrics.csv',sep=''),sep=',',header=TRUE) # read table of landscape metrics extracted from FRAGSTAT http://www.umass.edu/landeco/research/fragstats/fragstats.html
dt_metrics_scale <- cbind(dt_metrics[,c(1,2)], apply(dt_metrics[,c(3,4,5,6)], 2, function(x) scale(x, scale = TRUE))) # scale the landscape metrics

# transform data for plotting
mymelt_m <- NULL
for (i in seq(1,dim(dt_metrics_scale)[1])){
	for (j in seq(4,dim(dt_metrics_scale)[2])){
		mymelt_m <- rbind(mymelt_m, cbind(dt_metrics_scale[i,c(1,2)], names(dt_metrics_scale[j]), dt_metrics_scale[i,j]))
	}
}
names(mymelt_m) <- c('Survival_layer', 'Habitat_characteristics', 'Metrics', 'Value')


cbPalette <- c("#E69F00", "#0072B2")# define color
mymelt_m$x <- factor(mymelt_m$Habitat_characteristics, levels=c("Overall habitat", "highly suitable","less suitable","marginal areas")) # define factor order
mymelt_m$x 

ann_text <- data.frame(Survival_layer=factor("LS_1", levels=c("LS_1", "LS_2")), Habitat_characteristics=factor("Overall habitat", levels=c("Overall habitat", "highly suitable","less suitable","marginal areas")), Metrics=factor("PLAND", levels=c("PLAND", "NP","CONNECT")), lab='(C)', Value=1.2) # define annotation text
# Figure 3-C
t<-ggplot(mymelt_m, aes(x=x, y=Value, fill=Survival_layer)) +  geom_bar(stat = "identity", position = position_dodge(), colour="black") + scale_alpha_manual(values=c(0.3, 1)) +xlab('Habitat type')+ylab('Landscape metrics (normalised value)') + facet_grid(Metrics~., scale='free_y')+theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+ guides(color=guide_legend(title="Survival layer"))+theme_bw() + theme(legend.position='none', legend.text = element_text(size = 12),legend.title = element_text(size=14,face="bold"), axis.title.y = element_text(size = 12),axis.title.x = element_text(size = 12))+ scale_fill_manual(values=cbPalette) + geom_text(data = ann_text,aes(label =lab, x=0.6, y=1.5))
t 

grid.arrange(grid.arrange(pTime_area,pTime_lativeArea, nrow=2), t, nrow = 1, ncol =2)
