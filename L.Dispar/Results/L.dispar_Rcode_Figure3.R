###############################################################################
# R functions for analysing European gypsy moth, L. dispar, 
# range expansion in heterogenous environments 
###############################################################################
# Audrey Lustig
# Bio-Protection Research Centre, Lincoln University, NZ
# March 2016


path_to_mdig_casetudy_folder <- 'path_to_mdig_casetudy_folder/'
working_directory_path <- paste(path_to_mdig_casetudy_folder, 'L.Dispar/Data/RawLandscapeComponents',sep='') # working directory

# Plot map (Figure 3)
# read one of the replication
mat1<- read.table(paste(working_directory_path,'Frac.E.1280.R.10.H.0.3.P.35.rawASCIIgrid-rep1.asc',sep='/'), header=FALSE) 
mat2<- read.table(paste(working_directory_path,'Frac.E.1280.R.10.H.0.3.P.55.rawASCIIgrid-rep1.asc',sep='/'), header=FALSE) 
mat3<- read.table(paste(working_directory_path,'Frac.E.1280.R.10.H.0.3.P.75.rawASCIIgrid-rep1.asc',sep='/'), header=FALSE) 

mat4<- read.table(paste(working_directory_path,'Frac.E.1280.R.10.H.0.5.P.35.rawASCIIgrid-rep1.asc',sep='/'), header=FALSE) 
mat5<- read.table(paste(working_directory_path,'Frac.E.1280.R.10.H.0.5.P.55.rawASCIIgrid-rep1.asc',sep='/'), header=FALSE) 
mat6<- read.table(paste(working_directory_path,'Frac.E.1280.R.10.H.0.5.P.75.rawASCIIgrid-rep1.asc',sep='/'), header=FALSE) 

mat7<- read.table(paste(working_directory_path,'Frac.E.1280.R.10.H.0.7.P.35.rawASCIIgrid-rep1.asc',sep='/'), header=FALSE) 
mat8<- read.table(paste(working_directory_path,'Frac.E.1280.R.10.H.0.7.P.55.rawASCIIgrid-rep1.asc',sep='/'), header=FALSE) 
mat9<- read.table(paste(working_directory_path,'Frac.E.1280.R.10.H.0.7.P.75.rawASCIIgrid-rep1.asc',sep='/'), header=FALSE) 


# create layout and plot
layout(matrix(seq(9),nrow=3,ncol=3,byrow=TRUE),widths=c(1,1,1),heights=c(1,1,1))
par(mar=c(1,1,1,1), oma=c(3,3,3,3))
plot_ly(z = as.matrix(mat1), x = xx, y = yy, ,col=c('black','white'), xaxt='H=0.3',yaxt='P = 35',type = "heatmap")
plot_ly(z = as.matrix(mat2), x = xx, y = yy, ,col=c('black','white'), xaxt='H=0.3',yaxt='P = 55',type = "heatmap")

layout(matrix(seq(9),nrow=3,ncol=3,byrow=TRUE),widths=c(1,1,1),heights=c(1,1,1))
image(as.matrix(mat1),col=c('black','white'),  xaxt='n',yaxt='n', ylab='', useRaster=TRUE, cex=2)
mtext('P = 35', side=2, cex=1.5)
image(as.matrix(mat2),col=c('black','white'), xaxt='n',yaxt='n', useRaster=TRUE)
image(as.matrix(mat3),col=c('black','white'), xaxt='n',yaxt='n', useRaster=TRUE)

image(as.matrix(mat4),col=c('black','white'), xaxt='n',yaxt='n', ylab='', useRaster=TRUE)
mtext('P = 55', side=2, cex=1.5)
image(as.matrix(mat5),col=c('black','white'), xaxt='n',yaxt='n', useRaster=TRUE)
image(as.matrix(mat6),col=c('black','white'), xaxt='n',yaxt='n', useRaster=TRUE)

image(as.matrix(mat7),col=c('black','white'), xaxt='n',yaxt='n', ylab='', useRaster=TRUE)
mtext('P = 75', side=2, cex=1.5)
mtext('H = 0.1', side=1, cex=1.5, line=1)
image(as.matrix(mat8),col=c('black','white'), xaxt='n',yaxt='n', useRaster=TRUE)
mtext('H = 0.5', side=1, cex=1.5, line=1)
image(as.matrix(mat9),col=c('black','white'), xaxt='n',yaxt='n', useRaster=TRUE)
mtext('H = 0.9', side=1, cex=1.5, line=1)

