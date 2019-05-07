# ****************************************************************************
# MODEL CODE FOR:
#   "General rules for environmental management to prioritise social-ecological
#   systems research based on a value of information approach"
# AUTHORS:
#   Katrina J Davis, Iadine Chadès, Jonathan R Rhodes & Michael Bode
# DESCRIPTION:
#   Code recreates Figure 4 in the manuscript (but at lower discretisation)
# NOTE: 
#   This code should be run after 'Data_Aggregation'
# ****************************************************************************

# Update with personal file path
workDir <- "C:/Users/trins/Dropbox/0_Work/1_Workshops_shared/SNA Workshop/Projects/Value_of_information/Submission/6_JAE/Proofs/Code"

# Set working directory
setwd(workDir) 

# libraries
library(dplyr) # general data wrangling
library(ggplot2) # graphing
library(grDevices) # colours

# clear environment
rm(list=ls())

# Color settings ---------------------------------------------------------

# Colours
DBlue = "#1D7BE5" # Social
LBlue = "#00B0F0" # Social-ecological
LGreen = "#64C39B" # Ecological
White = "White"

# Order for figure: Soc, Soc, Ecol, Ecol, Soc-Ecol, space, Soc, Ecol, Soc-Ecol
color <- c(DBlue,DBlue,LGreen,LGreen,LBlue,White,DBlue,LGreen,LBlue)


# Aggregate systems data ------------------------------------------------------------------


# Read in systems data
S1 <- read.csv("S1.csv", header = F) # System 1
S2 <- read.csv("S2.csv", header = F) # System 2
S3 <- read.csv("S3.csv", header = F) # System 3
S4 <- read.csv("S4.csv", header = F) # System 4

# Combine data into single dataframe
d <- rbind(S1, S2, S3, S4)

# Specify column names
colnames(d) <- c("I","q","r","C","H","Soc","Eco","Soc-Eco")
# Add an empty column
d$No <- NA
# Re-order
d <- d[,c("I","q","r","C","H","No","Soc","Eco","Soc-Eco")]

# Data to percentages
d=d*100

# Create Manuscript Figure 4 as PDF ---------------------------------------------------------------------

# Individual plot titles
mains <- c("1. Territorial use rights fishery", 
           "2. Salmon fishery", 
           "3. Agricultural system", 
           "4. Non-timber forest products")

pdf(file = "EVPXI.pdf", width = 3.5, height = 8) # call pdf (units in inches)
par(mfrow=c(4,1), # multipane plot
    mar=c(0,5,2,0.5), # inner plot margins (space btw plots)
    oma=c(5,0,0,0)) # outer margins (space outside plots)

j <- 0 # initialise counter j
for (i in c(2,5,8,11)){
  j=j+1
  barCenters <- barplot(as.numeric(d[i,]), col=color,
                        names.arg = colnames(d),
                        xaxt="n",
                        axes=T,
                        ylim = c(0, 100),
                        main = mains[j],
                        ylab = "", cex.lab=1, cex.main=1.25, cex.names=0.5, cex.axis=1)
  segments(6.7, 0, 6.7, 85, lwd=1.5, lty=2, col="gray")
  
  if (j==4) {
    par(xpd=NA)
    text(x=barCenters[1:5],y = par("usr")[3] - 9,labels=colnames(d[1:5]), 
         cex = 1, pos = 1, offset = 0.01)
    text(x=barCenters[7:9],y = par("usr")[3] - 9,labels=colnames(d[7:9]), 
         cex = 1, pos = 1, offset = 0.3, srt=45)
  }
  
  segments(barCenters,as.matrix(d[i-1,]), barCenters, as.matrix(d[i+1,]), lwd = 1)
  arrows(barCenters, as.matrix(d[i-1,]), barCenters, as.matrix(d[i+1,]), lwd = 1, 
         angle = 90,code = 3, length = 0.05)
}
par(xpd=NA)

# y axis label
text(x = -2, y = 240, labels = "Expected value of partial information (%)", srt=90, cex=1.25)

dev.off() # close graph



# Sensitivity Analysis ----------------------------------------------------

#Read data
voiPath <- "C:/Users/trins/Dropbox/0_Work/1_Workshops_shared/SNA Workshop/Projects/Value_of_information/Code/Code_080617/SavedResults"

# M
mBase <- read.csv(paste(voiPath,"/SA_Base.csv", sep = ""),header=F)
mLower <- read.csv(paste(voiPath,"/MLower_0_225.csv", sep = ""),header=F)
mUpper <- read.csv(paste(voiPath,"/MUpper_0_275.csv", sep = ""),header=F)

# rbind all rows
mSA <- rbind(mLower, mBase, mUpper)

# Specify column names
colnames(mSA) <- c("I","q","r","C","H","Soc","Eco","Soc-Eco")

# Add an empty column
mSA$No <- NA
mSA <-  mSA[,c("I","q","r","C","H","No","Soc","Eco","Soc-Eco")]

# Data to percentages
mSA <- mSA * 100

# D
dBase <- read.csv(paste(voiPath,"/SA_Base.csv", sep = ""),header=F)
dLower <- read.csv(paste(voiPath,"/DLower_0_225.csv", sep = ""),header=F)
dUpper <- read.csv(paste(voiPath,"/DUpper_0_275.csv", sep = ""),header=F)

# rbind all rows
dSA <- rbind(dLower, dBase, dUpper)

# Specify column names
colnames(dSA) <- c("I","q","r","C","H","Soc","Eco","Soc-Eco")

# Add an empty column
dSA$No <- NA
dSA <-  dSA[,c("I","q","r","C","H","No","Soc","Eco","Soc-Eco")]

# Data to percentages
dSA <- dSA * 100



# Plot SA - M --------------------------------------------------------------------

# Bar plot of Low, Base and High
# pdf("VOI M_SA.pdf", w=8, h=20) #as pdf
jpeg(filename = "VOI_M_SA.jpeg", width = 550, height = 900) #as jpeg
# quartz(w=6,h=24)
par(mfrow=c(3,1), oma=c(10,2,1,2), mar=c(2,5,2,2)) #pdf
mains <- c("Low M (-10%)", "Base M", "High M (+10%)")

j=0
# NB, have 3 tables, dim: 9,10 so mid rows are 2,5,8. 
for (i in c(2,5,8)){
  j=j+1
  barCenters_M <- barplot(as.numeric(mSA[i,]), col=color,
                          names.arg = colnames(mSA),
                          xaxt="n",
                          #border="black",
                          axes=T,
                          ylim = c(0, 100),
                          main = mains[j],
                          ylab = "",cex.lab=2, cex.main=2.5, cex.names=1.5, cex.axis=2)
  segments(6.7, 0, 6.7, 85, lwd=2, lty=2, col="gray")
  
  if (j==3) {#text(x = barCenters, y = par("usr")[3] - 1.5, srt = 45, adj = 1, labels = c(colnames(d)[1:5],"",colnames(d)[7:9]), xpd = TRUE, cex=2.4)
    par(xpd=NA)
    text(x=barCenters_M[1:5],y = par("usr")[3] - 9,labels=colnames(mSA[1:5]),cex=2.4)
    text(x=barCenters_M[7:9],y = par("usr")[3] - 9,labels=colnames(mSA[7:9]),srt=45,cex=2.4)
  }
  # Error (quartile) lines
  segments(barCenters_M,as.matrix(mSA[i-1,]), barCenters_M, as.matrix(mSA[i+1,]), lwd = 1.5)
  # Error (quartile) bars
  arrows(barCenters_M, as.matrix(mSA[i-1,]), barCenters_M, as.matrix(mSA[i+1,]), lwd = 1.5, angle = 90,code = 3, length = 0.05)
}
par(xpd=NA)
text(-1.3,200,"Expected value of partial information (%)",srt=90,cex=3)

dev.off()


# Table -  SA M -----------------------------------------------------------


# TABLE
mSA_table <- mSA %>% 
  slice(c(2,5,8))

# Label columns
mSA_table$M <- c("lowM","baseM","highM")
# Reorder columns
mSA_table <- mSA_table %>% 
  select(M, everything(), -No)
# Export as table - to paste in excel or word
s2M_export <- format.data.frame(mSA_table) # Remove formatting so I can copy
write.table(s2M_export, 'clipboard', sep='\t',row.names=F)  # **For copy and paste to table in manuscript

# Difference relative to base case

mDiff <- mSA_table %>% 
  select(-M)

(mDiffL <- (mDiff[1,] - mDiff[2,])/mDiff[2,] * 100)
(mDiffH <- (mDiff[3,] - mDiff[2,])/mDiff[2,] * 100)

mDiff[1,1]
mDiff[2,1]
mDiff[3,1]

# Plot SA - D --------------------------------------------------------------------

# Bar plot of Low, Base and High
# pdf("VOI D_SA.pdf", w=8, h=20) #as pdf
tiff(filename = "VOI_D_SA.tiff", width = 550, height = 900) #as tiff
# quartz(w=6,h=24)
par(mfrow=c(3,1), oma=c(10,2,1,2), mar=c(2,5,2,2)) #pdf
mains <- c("Low D (-10%)", "Base D", "High D (+10%)")

j=0
# NB, have 3 tables, dim: 9,10 so mid rows are 2,5,8. 
for (i in c(2,5,8)){
  j=j+1
  barCenters_D <- barplot(as.numeric(dSA[i,]), col=color,
                          names.arg = colnames(dSA),
                          xaxt="n",
                          #border="black",
                          axes=T,
                          ylim = c(0, 100),
                          main = mains[j],
                          ylab = "",cex.lab=2, cex.main=2.5, cex.names=1.5, cex.axis=2)
  segments(6.7, 0, 6.7, 85, lwd=2, lty=2, col="gray")
  
  if (j==3) {
    par(xpd=NA)
    text(x=barCenters_D[1:5],y = par("usr")[3] - 9,labels=colnames(dSA[1:5]),cex=2.4)
    text(x=barCenters_D[7:9],y = par("usr")[3] - 9,labels=colnames(dSA[7:9]),srt=45,cex=2.4)
  }
  # Error (quartile) lines
  segments(barCenters_D,as.matrix(dSA[i-1,]), barCenters_D, as.matrix(dSA[i+1,]), lwd = 1.5)
  # Error (quartile) bars
  arrows(barCenters_D, as.matrix(dSA[i-1,]), barCenters_D, as.matrix(dSA[i+1,]), lwd = 1.5, angle = 90,code = 3, length = 0.05)
}
par(xpd=NA)
text(-1.3,200,"Expected value of partial information (%)",srt=90,cex=3)

dev.off()


# Table -  SA D -----------------------------------------------------------


# TABLE
dSA_table <- dSA %>% 
  slice(c(2,5,8))

# Label columns
dSA_table$D <- c("lowD","baseD","highD")
# Reorder columns
dSA_table <- dSA_table %>% 
  select(D, everything(), -No)
# Export as table - to paste in excel or word
s2D_export <- format.data.frame(dSA_table) # Remove formatting so I can copy
write.table(s2D_export, 'clipboard', sep='\t',row.names=F)  # **For copy and paste to table in manuscript

# Difference relative to base case

mDiff <- mSA_table %>% 
  select(-M)

(mDiffL <- (mDiff[1,] - mDiff[2,])/mDiff[2,] * 100)
(mDiffH <- (mDiff[3,] - mDiff[2,])/mDiff[2,] * 100)

mDiff[1,1]
mDiff[2,1]
mDiff[3,1]

# Next --------------------------------------------------------------------




# subset analysis just for system 2
dim(d) # 12 9
dim(d)[1] / 4 # rows per system = 3

# I want system 2, so rows 4-6
s2 <-  d[4:6,]

# The plot for just system 2

# Pull out the center observation: 2 (mean VOI)
barCenters_s2 <- barplot(as.numeric(s2[2,]), col=color,
                         names.arg = colnames(d),
                         xaxt="n",
                         #border="black",
                         axes=T,
                         ylim = c(0, 100),
                         main = mains[2],
                         ylab = "",cex.lab=2, cex.main=2.5, cex.names=1.5, cex.axis=2)
segments(6.7, 0, 6.7, 85, lwd=2, lty=2, col="gray") # Include line between single parameters and pairs
# Labels for parameters
par(xpd=NA)
text(x=barCenters_s2[1:5],y = par("usr")[3] - 9,labels=colnames(s2[1:5]),cex=2.4)
text(x=barCenters_s2[7:9],y = par("usr")[3] - 9,labels=colnames(s2[7:9]),srt=45,cex=2.4)
# Error bars
arrows(barCenters_s2, as.matrix(s2[2-1,]), barCenters_s2, as.matrix(s2[2+1,]), lwd = 1.5, angle = 90,code = 3, length = 0.05)
par(xpd=NA)
# Vertical axis label
text(-0.8,50,"Expected value of partial information (%)",srt=90,cex=2)



# pdf("VOI Fig 2.pdf", w=8, h=20) #as pdf
jpeg(filename = "VOI SA_M.jpeg", width = 550, height = 1200) #as jpeg
# quartz(w=6,h=24)
par(mfrow=c(4,1), oma=c(10,2,1,2), mar=c(2,5,2,2)) #pdf

barCenters_s2 <- barplot(as.numeric(s2[2,]), col=color,
                         names.arg = colnames(d),
                         xaxt="n",
                         #border="black",
                         axes=T,
                         ylim = c(0, 100),
                         main = mains[2],
                         ylab = "",cex.lab=2, cex.main=2.5, cex.names=1.5, cex.axis=2)
segments(6.7, 0, 6.7, 85, lwd=2, lty=2, col="gray") # Include line between single parameters and pairs
# Labels for parameters
par(xpd=NA)
text(x=barCenters_s2[1:5],y = par("usr")[3] - 9,labels=colnames(s2[1:5]),cex=2.4)
text(x=barCenters_s2[7:9],y = par("usr")[3] - 9,labels=colnames(s2[7:9]),srt=45,cex=2.4)
# Error bars
arrows(barCenters_s2, as.matrix(s2[2-1,]), barCenters_s2, as.matrix(s2[2+1,]), lwd = 1.5, angle = 90,code = 3, length = 0.05)
par(xpd=NA)
# Vertical axis label
text(-0.8,50,"Expected value of partial information (%)",srt=90,cex=2)
dev.off()
