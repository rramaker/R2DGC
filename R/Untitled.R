#Add colored rectangles
plot.new()
rect(0,0.6,1/3,0.867, col="indianred1")
rect(1/3,0.6,2/3,0.867, col="steelblue2")
rect(2/3,0.6,3/3,0.867, col="springgreen3")
#Label boxes
text(0.5,0.9, "Metabolomics Data Processing", cex=1.5)
text(1/6,0.85, "Raw Signal Filtering\nand Peak Calling\n(Chromatof)", pos=1,cex=0.8,)
text(3/6,0.85, "Multi-sample Peak\nAlignment and\nMetabolite Identification\n(R2DGC)", cex=0.8, pos=1)
text(5/6,0.85, "Data Normalization\nand Statistical\nInference\n(Metaboanalyst)", pos=1,cex=0.8,)
arrows(0.1,0.65,0.9,0.65,length = 0.1 ,lwd=3)
#Add box1 peaks
points(seq(0,0.2,(0.2-0)/8),(c(0,0.08,0.12,0.15,0.16,0.15,0.13,0.09,0)*0.5)+0.45,type="l", lwd=3)
points(seq(0.1,0.3,(0.3-0.1)/8),(c(0,0.08,0.12,0.15,0.16,0.15,0.13,0.09,0)*0.5)+0.45,type="l", lwd=3)
points(seq(0,0.15,(0.15-0)/8),(c(0,0.08,0.12,0.15,0.16,0.15,0.13,0.09,0)*0.5)+0.22,type="l", lwd=3)
points(seq(0.16,0.3,(0.3-0.16)/8),(c(0,0.08,0.12,0.15,0.16,0.15,0.13,0.09,0)*0.5)+0.22,type="l", lwd=3)
#Add box2 peaks
points(seq(0,0.13,(0.13-0)/8)+(0.36),(c(0.0,0.08,0.12,0.15,0.16,0.15,0.13,0.09,0)/3)+0.26,type="l", lwd=3)
points(seq(0.17,0.3,(0.3-0.17)/8)+(0.34),(c(0,0.08,0.12,0.15,0.16,0.15,0.13,0.09,0)/3)+0.26,type="l", lwd=3)
points(seq(0,0.13,(0.13-0)/8)+(0.36),(c(0,0.08,0.12,0.15,0.16,0.15,0.13,0.09,0)/3)+0.32,type="l", lwd=3)
points(seq(0.17,0.3,(0.3-0.17)/8)+(0.34),(c(0,0.08,0.12,0.15,0.16,0.15,0.13,0.09,0)/3)+0.32,type="l", lwd=3)
points(seq(0,0.13,(0.13-0)/8)+(0.36),(c(0,0.08,0.12,0.15,0.16,0.15,0.13,0.09,0)/3)+0.38,type="l", lwd=3)
points(seq(0.17,0.3,(0.3-0.17)/8)+(0.34),(c(0,0.08,0.12,0.15,0.16,0.15,0.13,0.09,0)/3)+0.38,type="l", lwd=3)
points(seq(0,0.13,(0.13-0)/8)+(0.36),(c(0,0.08,0.12,0.15,0.16,0.15,0.13,0.09,0)/3)+0.44,type="l", lwd=3)
points(seq(0.17,0.3,(0.3-0.17)/8)+(0.34),(c(0,0.08,0.12,0.15,0.16,0.15,0.13,0.09,0)/3)+0.44,type="l", lwd=3)
#Add box3 peaks
points(seq(0,0.13,(0.13-0)/8)+(0.7),(c(0.0,0.08,0.12,0.15,0.16,0.15,0.13,0.09,0)/3)+0.2,type="l", lwd=3)
points(seq(0.17,0.3,(0.3-0.17)/8)+(0.69),(c(0,0.08,0.12,0.15,0.16,0.15,0.13,0.09,0)/3)+0.2,type="l", lwd=3)
points(seq(0,0.13,(0.13-0)/8)+(0.7),(c(0,0.08,0.12,0.15,0.16,0.15,0.13,0.09,0)/3)+0.26,type="l", lwd=3)
points(seq(0.17,0.3,(0.3-0.17)/8)+(0.69),(c(0,0.08,0.12,0.15,0.16,0.15,0.13,0.09,0)/3)+0.26,type="l", lwd=3)
points(seq(0,0.13,(0.13-0)/8)+(0.7),(c(0,0.08,0.12,0.15,0.16,0.15,0.13,0.09,0)/3)+0.42,type="l", lwd=3)
points(seq(0.17,0.3,(0.3-0.17)/8)+(0.69),(c(0,0.08,0.12,0.15,0.16,0.15,0.13,0.09,0)/3)+0.42,type="l", lwd=3)
points(seq(0,0.13,(0.13-0)/8)+(0.7),(c(0,0.08,0.12,0.15,0.16,0.15,0.13,0.09,0)/3)+0.48,type="l", lwd=3)
points(seq(0.17,0.3,(0.3-0.17)/8)+(0.69),(c(0,0.08,0.12,0.15,0.16,0.15,0.13,0.09,0)/3)+0.48,type="l", lwd=3)
rect(xleft = 0.69,ybottom = 0.19,xright = 1,0.325, lwd = 3)
rect(xleft = 0.69,ybottom = 0.41,xright = 1,0.546, lwd = 3)
text(0.85, 0.565, "Group A")
text(0.85, 0.345, "Group B")
#Add dashed lines between peaks
arrows(1/3,0.2,1/3,0.6,code = 0,lty=2)
arrows(2/3,0.2,2/3,0.6,code = 0,lty=2)
arrows(1/6.5,0.28,1/6.5,0.45,code = 1, length=0.1,lwd=2)
text(0.5, 0, "Figure 1. Overview of metabolomics data processing",cex=0.8)

########R2DGC basic workflow

plot.new()
text(0.5,1,"Input Sample\nChromatof File Paths", cex=0.8)
arrows(0.5,0.96,0.5,0.905,length=0.1, lwd=2)
rect(0.3,0.75,0.7,0.9, col="springgreen3")
text(0.5,0.83,"Ion Filtering\n(FindProblemIons)")
arrows(0.5,0.745,0.5,0.705,length=0.1, lwd=2)
rect(0.3,0.55,0.7,0.7, col="steelblue2")
text(0.5,0.63,"Intra-Sample\nPeak Compression\n(PrecompressFiles)")
arrows(0.5,0.545,0.5,0.505,length=0.1, lwd=2)
rect(0.3,0.3,0.7,0.5, col="indianred1")
text(0.5,0.4,"Multi-Sample Peak\nAlignment and\nMetabolite Identification\n(ConsensusAlign)")
rect(0.8,0.5,1,0.65, col="plum3")
text(0.9,0.58,"Standard Library\nCreation\n(MakeReference)", cex=0.7)
text(0.9,0.77,"Input Metabolite\nStandard\nChromatof File Paths", cex=0.8)
arrows(0.9,0.71,0.9,0.655,length=0.1, lwd=2)
arrows(0.9,0.495,0.705,0.4,length=0.1, lwd=2)
arrows(0.5,0.295,0.5,0.255,length=0.1, lwd=2)
rect(-0.02,0.03,1.02,0.25,col="snow2")
text(0.5,0.225,"Output List")
rect(0,0.05,0.3,0.2, col="wheat1")
text(0.15,0.13,"Peak Alignment Table", cex=0.8)
rect(0.35,0.05,0.65,0.2, col="khaki2")
text(0.5,0.13,"Peak Info Table", cex=0.8)
rect(0.7,0.05,1,0.2, col="lightgoldenrod")
text(0.85,0.13,"Incongruent Quant\nMass List", cex=0.8)
text(0.5, 0, "Figure 2. Basic Workflow. Functions shown in parenthesis",cex=0.8)


a<-c(1,10,300)
b<-c(300,10,1)
a%*%b/(sqrt(sum(a*a))*sqrt(sum(b*b)))

sqrt(sum(a*a))
sqrt(sum(b*b))

sqrt(sum(a*a))*sqrt(sum(b*b))

plot.new()
points(c(.1,.35,.3,.8),rep(0.8,4), pch=16, col=c("black","red","black","black"))
text(c(.1,.35,.5,.8),rep(0.85,5), labels = c("Std1","Analyte","Std2","Std3"))
arrows(x0 = 0, y0 = 0.7,x1 = 1,y1 = 0.7, code=0, lwd=3)
arrows(x0 = c(0,0.2,0.4,0.6,0.8,1), y0 = rep(0.7,6),x1 = c(0,0.2,0.4,0.6,0.8,1),y1 = rep(0.65,6), code=0, lwd=3)
text(c(0,0.2,0.4,0.6,0.8,1),rep(0.62,6), labels = c(seq(0,10,2)))
text(0.5,0.72, labels = "Retention Time",cex=0.75)

arrows(x0 = .1, y0 = 0.5,x1 = .8,y1 = 0.5, code=0, lwd=2)
text(0.45,0.52, labels = substitute(paste(italic("T"),"=7")))
arrows(x0 = c(.1,.8), y0 = rep(0.49,2),x1 = c(.1,.8),y1 = rep(0.51,2), code=0, lwd=2)

arrows(x0 = .1, y0 = 0.4,x1 = .35,y1 = 0.4, code=0, lwd=2)
text(.225,0.43, labels = expression(italic("d"[1])*"= -2.5"))
arrows(x0 = c(.1,.35), y0 = rep(0.39,2),x1 = c(.1,.35),y1 = rep(0.41,2), code=0, lwd=2)

arrows(x0 = .5, y0 = 0.3,x1 = .35,y1 = 0.3, code=0, lwd=2)
text(.425,0.33, labels = expression(italic("d"[2])*"= 1.5"))
arrows(x0 = c(.5,.35), y0 = rep(0.29,2),x1 = c(.5,.35),y1 = rep(0.31,2), code=0, lwd=2)

arrows(x0 = .8, y0 = 0.2,x1 = .35,y1 = 0.2, code=0, lwd=2)
text(.575,0.23, labels = expression(italic("d"[3])*"= 4.5"))
arrows(x0 = c(.8,.35), y0 = rep(0.19,2),x1 = c(.8,.35),y1 = rep(0.21,2), code=0, lwd=2)
text(.5,0, labels = "Figure 4. Computing retention index values", cex=0.75)


plot.new()
points(c(.1,.35,.3,.5),rep(0.8,4), pch=16, col=c("black","red","black","black"))
text(c(.1,.35,.3,.5),c(0.85,.75,.85,0.85), labels = c("Std1","","Std2","Std3"), cex=0.8)
rect(xleft = 0.05,ybottom = 0.75,xright = 0.55,ytop = 0.88)
text(0.3,0.9,"Column1")
points(c(.2,.4,.6),rep(0.6,3), pch=16, col=c("black","black","black"))
text(c(.2,.4,.6),rep(0.65,5), labels = c("Std1","Std2","Std3"),cex=0.8)
rect(xleft = 0.15,ybottom = 0.55,xright = 0.65,ytop = 0.68)
text(0.4,0.7,"Column2")
points(c(.1,.5,1),rep(0.4,3), pch=16, col=c("black","black","black"))
text(c(.1,.5,1),rep(0.45,5), labels = c("Std1","Std2","Std3"),cex=0.8)
rect(xleft = 0.05,ybottom = 0.35,xright = 1.03,ytop = 0.48)
text(0.55,0.5,"Column3")
arrows(x0 = 0, y0 = 0.3,x1 = 1,y1 = 0.3, code=0, lwd=3)
arrows(x0 = c(0,0.2,0.4,0.6,0.8,1), y0 = rep(0.3,6),x1 = c(0,0.2,0.4,0.6,0.8,1),y1 = rep(0.25,6), code=0, lwd=3)
text(c(0,0.2,0.4,0.6,0.8,1),rep(0.22,6), labels = c(seq(0,10,2)))
text(0.5,0.32, labels = "Retention Time",cex=0.75)
text(.5,.05, labels = "Figure 5. Testing retention index values", cex=0.75)


FAME_Standards<-c(1,3,5)
Total_Standard_Length<-(max(FAME_Standards)-min(FAME_Standards))
Analyte1_RT<-3.5
Analyte1_RetentionIndex<- (FAME_Standards-Analyte1_RT)/Total_Standard_Length

IndexDifferences<-c()
for(Analyte2_RT in seq(0,10,0.1)){
  Analyte2_RetentionIndex<- (FAME_Standards-Analyte2_RT)/Total_Standard_Length
  IndexDifferences[as.character(Analyte2_RT)]<-sum(abs(Analyte2_RetentionIndex-Analyte1_RetentionIndex))
}
plot(seq(0,10,0.1),IndexDifferences, pch=16, xlab="Analyte 2 Retention Time", ylab="Difference")

FAME_Standards1<-c(1,3,5)
Total_Standard_Length1<-(max(FAME_Standards1)-min(FAME_Standards1))
Analyte1_RT<-3.5
Analyte1_RetentionIndex<- (FAME_Standards1-Analyte1_RT)/Total_Standard_Length
FAME_Standards2<-c(2,4,6)
Total_Standard_Length2<-(max(FAME_Standards2)-min(FAME_Standards2))
IndexDifferences<-c()
for(Analyte2_RT in seq(0,10,0.1)){
  Analyte2_RetentionIndex<- (FAME_Standards2-Analyte2_RT)/Total_Standard_Length2
  IndexDifferences[as.character(Analyte2_RT)]<-sum(abs(Analyte2_RetentionIndex-Analyte1_RetentionIndex))
}
plot(seq(0,10,0.1),IndexDifferences, pch=16, xlab="Analyte 2 Retention Time", ylab="Retention Index Difference")


FAME_Standards1<-c(1,3,5)
Total_Standard_Length1<-(max(FAME_Standards1)-min(FAME_Standards1))
Analyte1_RT<-3.5
Analyte1_RetentionIndex<- (FAME_Standards1-Analyte1_RT)/Total_Standard_Length
FAME_Standards2<-c(1,5,10)
Total_Standard_Length2<-(max(FAME_Standards2)-min(FAME_Standards2))
IndexDifferences<-c()
for(Analyte2_RT in seq(0,10,0.1)){
  Analyte2_RetentionIndex<- (FAME_Standards2-Analyte2_RT)/Total_Standard_Length2
  IndexDifferences[as.character(Analyte2_RT)]<-sum(abs(Analyte2_RetentionIndex-Analyte1_RetentionIndex))
}
plot(seq(0,10,0.1),IndexDifferences, pch=16, xlab="Analyte 2 Retention Time", ylab="Retention Index Difference")


plot.new()
rect(xleft = 0.25,ybottom = 0.9,xright = 0.75,ytop = 1, col="springgreen3")
text(0.5,0.95,"Read in sample files", cex=0.75)
arrows(x0 = 0.5,y0 = 0.9,x1 = 0.5,y1 = 0.875,length=0.05, lwd=2)

rect(0.25,0.775,0.75,0.875, col="wheat1")
text(0.5,0.825,"Optional: Compute retention indices", cex=0.75)
arrows(x0 = 0.5,y0 = 0.775,x1 = 0.5,y1 = 0.75,length=0.05, lwd=2)

rect(0.25,0.65,0.75,.75, col="indianred1")
text(0.5,0.7,"Compute pairwise sample-seed\npeak similarity scores", cex=0.75)
arrows(x0 = 0.5,y0 = 0.65,x1 = 0.5,y1 = 0.625,length=0.05, lwd=2)

rect(0.25,0.525,0.75,.625, col="khaki2")
text(0.5,0.575,"Optional: compute optimal\npeak similarity threshold", cex=0.75)
arrows(x0 = 0.5,y0 = 0.525,x1 = 0.5,y1 = 0.5,length=0.05, lwd=2)

rect(0.25,0.4,0.75,.5, col="steelblue2")
text(0.5,0.45,"Find best peak pairs above\npeak similarity threshold", cex=0.75)
arrows(x0 = 0.5,y0 = 0.4,x1 = 0.5,y1 = 0.375,length=0.05, lwd=2)

rect(0.25,0.275,0.75,.375, col="plum3")
text(0.5,0.325,"Optional: Relaxed threshold search\nfor high likelihood missing peaks", cex=0.75)
arrows(x0 = 0.5,y0 = 0.275,x1 = 0.5,y1 = 0.25,length=0.05, lwd=2)

rect(0.25,0.15,0.75,.25, col="lightgoldenrod")
text(0.5,0.2,"Optional: Identify aligned\npeaks with reference library", cex=0.75)
arrows(x0 = 0.5,y0 = 0.15,x1 = 0.5,y1 = 0.125,length=0.05, lwd=2)

arrows(x0 = 0.1,y0 = 0.325,x1 = 0.1,y1 = 0.7, lwd=2, code=0)
arrows(x0 = 0.25,y0 = 0.325,x1 = 0.1,y1 = 0.325, lwd=2, code=0)
arrows(x0 = 0.25,y0 = 0.7,x1 = 0.1,y1 = 0.7, lwd=2, length=0.05, code=1)
rect(0,0.45,0.2,0.575, col="coral")
text(0.1,0.515, "Optional: repeat\nalignment with\nmultiple seeds", cex=0.7)

rect(0,0,0.3,0.09, col="cadetblue2")
text(0.15,0.045,"Peak Alignment Table", cex=0.8)
rect(0.35,0.0,0.65,0.09, col="navajowhite")
text(0.5,0.045,"Peak Info Table", cex=0.8)
rect(0.7,0.0,1,0.09, col="pink2")
text(0.85,0.045,"Incongruent Quant\nMass List", cex=0.8)
text(0.5,.105, "Outputs", cex=0.8)
rect(0,0,1,0.09)
text(0.5,-0.02, "Figure 9. Overview of ConsensusAlign function.", cex=0.7)
