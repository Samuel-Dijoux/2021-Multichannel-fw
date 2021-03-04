###################################################################################### 
#                -> Bifurcation analyses over K1 and Beta gradients for varying K2 
#
# Codimension bifurcation tracking selected BP
# as a function of two parameters : K2 and K1+Beta
# Representation of the coexistence phase of both consumers: P-C1C2-R1R2 ~ K1+Beta            
######################################################################################
library(PSPManalysis)
if (!exists("PSPMsrc.fullpath")) source("../PSPMworkshop/PSPManalysis/PSPManalysis.r")

################## IA. K2=1e-5 ========================
### EQ ========================
### with Beta = 0.5
cmd = output8IAXOB <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[7,2], 1E-8, output2$curvepoints[7,c(5,7)]),
                               0.1, c(11, 0, 4E-3), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = output8IAXOBb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IAXOB$bifpoints[1,c(1:3,6:8)]),
                                0.1, c(11, 0, 4E-3), NULL, c("envZE", "2","envZE", "4"), clean=FALSE)
cmd = output8IAXOBc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IAXOBb$bifpoints[,c(1:8)]),
                                -0.5, c(11, 0, 4E-3), NULL, NULL, clean=FALSE)
### with Beta = 2.5, detection of LP2
cmd = output8IAXOBd <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IAXOBc$bifpoints[1,c(1:4,6,8)]),
                                -0.5, c(11, 0, 4E-3), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
### No detection of LP3 with Beta = 3.0 & 3.5

Eq_K1_IAXOB_BP1 <- output8IAXOB$bifpoints[1,]; #inv of C1 in C2-R2+R1
Eq_K1_IAXOB_BP2 <- output8IAXOBc$bifpoints[2,];
Eq_K1_IAXOB_BPE <- output8IAXOBb$bifpoints[1,];
Eq_K1_IAXOB_LP <- output8IAXOBc$bifpoints[1,];
Eq_K1_IAXOB_LPb <- output8IAXOBd$bifpoints[2,];
### Codim ========================
### BP 1 
cmd = output14IAXOBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IAXOB$bifpoints[1,c(1:3,7)], 0.5), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "2", "envZE", "3", "envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output14IAXOBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IAXOB$bifpoints[,c(1:3,7)], 0.5), 0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "2", "envZE", "3", "envZE", "4", "popBP", "1"), clean=FALSE)
output14IAXOB<-rbind(output14IAXOBa$curvepoints[dim(output14IAXOBa$curvepoints)[1]:1,],
                     output14IAXOBb$curvepoints[2:dim(output14IAXOBb$curvepoints)[1],]);
### BP 2
cmd = output15IAXOBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IAXOBc$bifpoints[2,c(1:6,8)], 0.5), -0.1, c(11, 0, 9.5E-1, 20, 0.1, 4.5), NULL, c("popBP", "0"), clean=FALSE)
cmd = output15IAXOBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IAXOBc$bifpoints[2,c(1:6,8)], 0.5), 0.1, c(11, 0, 9.5E-1, 20, 0.1, 4.5), NULL, c("popBP", "0"), clean=FALSE)
output15IAXOB<-rbind(output15IAXOBa$curvepoints[dim(output15IAXOBa$curvepoints)[1]:1,],
                     output15IAXOBb$curvepoints[2:dim(output15IAXOBb$curvepoints)[1],]);
### BPE
cmd = output16IAXOBa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output8IAXOBb$bifpoints[1,c(1:3,7,8)], 0.5), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)
cmd = output16IAXOBb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output8IAXOBb$bifpoints[1,c(1:3,7,8)], 0.5), 0.1, c(11, 0, 9.5E-01, 20, -0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)
output16IAXOB<-rbind(output16IAXOBa$curvepoints[dim(output16IAXOBa$curvepoints)[1]:1,],
                     output16IAXOBb$curvepoints[2:dim(output16IAXOBb$curvepoints)[1],]);
### LP
cmd = output17IAXOBa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IAXOBc$bifpoints[1,c(1:8)], 0.5), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, NULL, clean=FALSE)
cmd = output17IAXOBb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IAXOBc$bifpoints[1,c(1:8)], 0.5), 0.1, c(11, 0, 9.5E-01, 20, -0.1, 4.5), NULL, NULL, clean=FALSE)
output17IAXOB<-rbind(output17IAXOBa$curvepoints[dim(output17IAXOBa$curvepoints)[1]:1,],
                     output17IAXOBb$curvepoints[2:dim(output17IAXOBb$curvepoints)[1],]);

cmd = output17IAXOBc <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IAXOBd$bifpoints[2,c(1:4,6,8)], 2.5), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output17IAXOBd <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IAXOBd$bifpoints[2,c(1:4,6,8)], 2.5), 1.2, c(11, 0, 9.5E-01, 20, -0.1, 4.5), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
output17IAXOB2<-rbind(output17IAXOBc$curvepoints[dim(output17IAXOBc$curvepoints)[1]:1,],
                      output17IAXOBd$curvepoints[2:dim(output17IAXOBd$curvepoints)[1],]);

################## IB. K2=2e-5 ========================
### EQ ========================
cmd = output8IAXIA <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[21,2], 1E-8, output2$curvepoints[21,c(5,7)]),
                               0.1, c(11, 0, 4E-4), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = output8IAXIAb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IAXIA$bifpoints[1,c(1:3,6:8)]),
                                0.1, c(11, 0, 4E-4), NULL, c("envZE", "2","envZE", "4"), clean=FALSE)
cmd = output8IAXIAc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IAXIAb$bifpoints[,c(1:8)]),
                                -0.5, c(11, 0, 4E-4), NULL, NULL, clean=FALSE)
### with Beta = 2.5, detection of LP2
cmd = output8IAXIAd <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output3.2$bifpoints[2,c(1:8)]),
                                0.8, c(11, 0, 4E-3), NULL, NULL, clean=FALSE)
### with Beta = 30, detection of LP3
cmd = output8IAXIAe <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[21,2], 1E-8, output2$curvepoints[21,c(5,7)]),
                                0.1, c(11, 0, 4E-4), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = output8IAXIAf <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8ICXIBe$bifpoints[1,c(1:3,6:8)]),
                                0.1, c(11, 0, 5E-2), NULL, c("envZE", "2","envZE", "4"), clean=FALSE)
cmd = output8IAXIAg <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8ICXIBf$bifpoints[,c(1:8)]),
                                -0.5, c(11, 0, 5E-2), NULL, NULL, clean=FALSE)

Eq_K1_IAXIA_BP1 <- output8IAXIA$bifpoints[1,]; #inv of C1 in C2-R2+R1
Eq_K1_IAXIA_BP2 <- output8IAXIAc$bifpoints[3,];
Eq_K1_IAXIA_BPE <- output8IAXIAb$bifpoints[1,];
Eq_K1_IAXIA_LP <- output8IAXIAc$bifpoints[2,];
Eq_K1_IAXIA_LPb <- output3.2$bifpoints[3,];
Eq_K1_IAXIA_LPc <- output8IAXIAg$bifpoints[2,];
### Codim ========================
### BP 1 
cmd = output14IBXIAa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IAXIA$bifpoints[,c(1:3,7)], 1.0), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "2", "envZE", "3", "envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output14IBXIAb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IAXIA$bifpoints[,c(1:3,7)], 1.0), 0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "2", "envZE", "3", "envZE", "4", "popBP", "1"), clean=FALSE)
output14IBXIA<-rbind(output14IBXIAa$curvepoints[dim(output14IBXIAa$curvepoints)[1]:1,],
                     output14IBXIAb$curvepoints[2:dim(output14IBXIAb$curvepoints)[1],]);
### BP 2
cmd = output15IBXIAa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IAXIAc$bifpoints[3,c(1:6,8)], 1.0), -0.1, c(11, 0, 9.5E-1, 20, 0.1, 4.5), NULL, c("popBP", "0"), clean=FALSE)
cmd = output15IBXIAb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IAXIAc$bifpoints[3,c(1:6,8)], 1.0), 0.1, c(11, 0, 9.5E-1, 20, 0.1, 4.5), NULL, c("popBP", "0"), clean=FALSE)
output15IBXIA<-rbind(output15IBXIAa$curvepoints[dim(output15IBXIAa$curvepoints)[1]:1,],
                     output15IBXIAb$curvepoints[2:dim(output15IBXIAb$curvepoints)[1],]);
### BPE
cmd = output16IBXIAa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output8IAXIAb$bifpoints[1,c(1:3,7,8)], 1.0), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)
cmd = output16IBXIAb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output8IAXIAb$bifpoints[1,c(1:3,7,8)], 1.0), 0.1, c(11, 0, 9.5E-01, 20, -0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)
output16IBXIA<-rbind(output16IBXIAa$curvepoints[dim(output16IBXIAa$curvepoints)[1]:1,],
                     output16IBXIAb$curvepoints[2:dim(output16IBXIAb$curvepoints)[1],]);
### LP
cmd = output17IBXIAa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IAXIAc$bifpoints[2,c(1:8)], 1.0), -0.1, c(11, 0, 9.5E-01, 20, -0.1, 4.5), NULL, NULL, clean=FALSE)
cmd = output17IBXIAb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IAXIAc$bifpoints[2,c(1:8)], 1.0), 0.1, c(11, 0, 9.5E-01, 20, -0.1, 4.5), NULL, NULL, clean=FALSE)
output17IBXIA<-rbind(output17IBXIAa$curvepoints[dim(output17IBXIAa$curvepoints)[1]:1,],
                     output17IBXIAb$curvepoints[2:dim(output17IBXIAb$curvepoints)[1],]);

cmd = output17IBXIAc <- PSPMequi("Tritrophicmod5_fin", "LP", c(output3.2$bifpoints[3,c(1:4,6,8)], 2.5), -0.1, c(11, 0, 9.5E-01, 20, -0.1, 4.5), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output17IBXIAd <- PSPMequi("Tritrophicmod5_fin", "LP", c(output3.2$bifpoints[3,c(1:4,6,8)], 2.5), 1.2, c(11, 0, 0.5E-01, 20, -0.1, 4.5), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
output17IBXIA2<-rbind(output17IBXIAc$curvepoints[dim(output17IBXIAc$curvepoints)[1]:1,],
                      output17IBXIAd$curvepoints[2:dim(output17IBXIAd$curvepoints)[1],]);

cmd = output17IBXIAe <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IAXIAg$bifpoints[2,c(1:8)], 3.5), -0.1, c(11, 0, 9.5E-01, 20, -0.1, 4.5), NULL, NULL, clean=FALSE)
cmd = output17IBXIAf <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IAXIAg$bifpoints[2,c(1:8)], 3.5), 1.2, c(11, 0, 0.5E-01, 20, -0.1, 4.5), NULL, NULL, clean=FALSE)
output17IBXIA3<-rbind(output17IBXIAe$curvepoints[dim(output17IBXIAe$curvepoints)[1]:1,],
                      output17IBXIAf$curvepoints[2:dim(output17IBXIAf$curvepoints)[1],]);
################## IC. K2=3e-5 ========================
### EQ ========================
cmd = output8ICXIB <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[27,2], 1E-8, output2$curvepoints[27,c(5,7)]),
                               0.1, c(11, 0, 4E-4), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = output8ICXIBb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8ICXIB$bifpoints[1,c(1:3,6:8)]),
                                0.1, c(11, 0, 4E-4), NULL, c("envZE", "2","envZE", "4"), clean=FALSE)
cmd = output8ICXIBc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8ICXIBb$bifpoints[,c(1:8)]),
                                -0.5, c(11, 0, 4E-4), NULL, NULL, clean=FALSE)
### with Beta = 2.5, detection of LP2
cmd = output8ICXIBd <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output3.2$bifpoints[1,c(1:8)]),
                                0.8, c(11, 0, 4E-3), NULL, NULL, clean=FALSE)
### with Beta = 3.5, detection of LP3
cmd = output8ICXIBe <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[27,2], 1E-8, output2$curvepoints[27,c(5,7)]),
                                0.1, c(11, 0, 4E-4), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = output8ICXIBf <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8ICXIBe$bifpoints[1,c(1:3,6:8)]),
                                0.1, c(11, 0, 5E-2), NULL, c("envZE", "2","envZE", "4"), clean=FALSE)
cmd = output8ICXIBg <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8ICXIBf$bifpoints[,c(1:8)]),
                                -0.5, c(11, 0, 5E-2), NULL, NULL, clean=FALSE)

Eq_K1_ICXIB_BP1 <- output8ICXIB$bifpoints[1,]; #inv of C1 in C2-R2+R1
Eq_K1_ICXIB_BP2 <- output8ICXIBc$bifpoints[3,];
Eq_K1_ICXIB_BPE <- output8ICXIBb$bifpoints[1,];
Eq_K1_ICXIB_LP <- output8ICXIBc$bifpoints[2,];
Eq_K1_ICXIB_LPb <- output3.2$bifpoints[2,];
Eq_K1_ICXIB_LPc <- output8ICXIBg$bifpoints[2,];
### Codim ========================
### BP 1 
cmd = output14ICXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8ICXIB$bifpoints[,c(1:3,7)], 1.0), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "2", "envZE", "3", "envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output14ICXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8ICXIB$bifpoints[,c(1:3,7)], 1.0), 0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "2", "envZE", "3", "envZE", "4", "popBP", "1"), clean=FALSE)
output14ICXIB<-rbind(output14ICXIBa$curvepoints[dim(output14ICXIBa$curvepoints)[1]:1,],
                     output14ICXIBb$curvepoints[2:dim(output14ICXIBb$curvepoints)[1],]);
### BP 2
cmd = output15ICXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8ICXIBc$bifpoints[3,c(1:6,8)], 1.0), -0.1, c(11, 0, 9.5E-1, 20, 0.1, 4.5), NULL, c("popBP", "0"), clean=FALSE)
cmd = output15ICXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8ICXIBc$bifpoints[3,c(1:6,8)], 1.0), 0.1, c(11, 0, 9.5E-1, 20, 0.1, 4.5), NULL, c("popBP", "0"), clean=FALSE)
output15ICXIB<-rbind(output15ICXIBa$curvepoints[dim(output15ICXIBa$curvepoints)[1]:1,],
                     output15ICXIBb$curvepoints[2:dim(output15ICXIBb$curvepoints)[1],]);
### BPE
cmd = output16ICXIBa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output8ICXIBb$bifpoints[1,c(1:3,7,8)], 1.0), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)
cmd = output16ICXIBb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output8ICXIBb$bifpoints[1,c(1:3,7,8)], 1.0), 0.1, c(11, 0, 9.5E-01, 20, -0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)
output16ICXIB<-rbind(output16ICXIBa$curvepoints[dim(output16ICXIBa$curvepoints)[1]:1,],
                     output16ICXIBb$curvepoints[2:dim(output16ICXIBb$curvepoints)[1],]);
### LP
cmd = output17ICXIBa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8ICXIBc$bifpoints[2,c(1:8)], 1.0), -0.1, c(11, 0, 9.5E-01, 20, -0.1, 4.5), NULL, NULL, clean=FALSE)
cmd = output17ICXIBb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8ICXIBc$bifpoints[2,c(1:8)], 1.0), 0.1, c(11, 0, 9.5E-01, 20, -0.1, 4.5), NULL, NULL, clean=FALSE)
output17ICXIB<-rbind(output17ICXIBa$curvepoints[dim(output17ICXIBa$curvepoints)[1]:1,],
                     output17ICXIBb$curvepoints[2:dim(output17ICXIBb$curvepoints)[1],]);
cmd = output17ICXIBc <- PSPMequi("Tritrophicmod5_fin", "LP", c(output3.2$bifpoints[2,c(1:4,6,8)], 2.5), -0.1, c(11, 0, 9.5E-01, 20, -0.1, 4.5), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output17ICXIBd <- PSPMequi("Tritrophicmod5_fin", "LP", c(output3.2$bifpoints[2,c(1:4,6,8)], 2.5), 1.2, c(11, 0, 0.5E-01, 20, -0.1, 4.5), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
output17ICXIB2<-rbind(output17ICXIBc$curvepoints[dim(output17ICXIBc$curvepoints)[1]:1,],
                      output17ICXIBd$curvepoints[2:dim(output17ICXIBd$curvepoints)[1],]);
cmd = output17ICXIBe <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8ICXIBg$bifpoints[2,c(1:8)], 3.5), -0.1, c(11, 0, 9.5E-01, 20, -0.1, 4.5), NULL, NULL, clean=FALSE)
cmd = output17ICXIBf <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8ICXIBg$bifpoints[2,c(1:8)], 3.5), 1.2, c(11, 0, 0.5E-01, 20, -0.1, 4.5), NULL, NULL, clean=FALSE)
output17ICXIB3<-rbind(output17ICXIBe$curvepoints[dim(output17ICXIBe$curvepoints)[1]:1,],
                      output17ICXIBf$curvepoints[2:dim(output17ICXIBf$curvepoints)[1],]);
################## IIA. K2=5e-5 ========================
### EQ ========================
cmd = output8IIAXIB <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[37,2], 1E-8, output2$curvepoints[37,c(5,7)]),
                                0.1, c(11, 0, 4E-4), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = output8IIAXIBb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IIAXIB$bifpoints[1,c(1:3,6:8)]),
                                 0.1, c(11, 0, 4E-4), NULL, c("envZE", "2","envZE", "4"), clean=FALSE)
cmd = output8IIAXIBc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IIAXIBb$bifpoints[,c(1:8)]),
                                 -0.5, c(11, 0, 4E-4), NULL, NULL, clean=FALSE)
### with Beta = 2.5, detection of LP2
cmd = output8IIAXIBd <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output3.2$bifpoints[1,c(1:8)]),
                                 0.8, c(11, 0, 4E-3), NULL, NULL, clean=FALSE)
### with Beta = 3.5, detection of LP3
cmd = output8IIAXIBe <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[37,2], 1E-8, output2$curvepoints[37,c(5,7)]),
                                 0.1, c(11, 0, 4E-4), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = output8IIAXIBf <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IIAXIBe$bifpoints[1,c(1:3,6:8)]),
                                 0.1, c(11, 0, 5E-2), NULL, c("envZE", "2","envZE", "4"), clean=FALSE)
cmd = output8IIAXIBg <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IIAXIBf$bifpoints[,c(1:8)]),
                                 -0.5, c(11, 0, 5E-2), NULL, NULL, clean=FALSE)

Eq_K1_IIAXIB_BP1 <- output8IIAXIB$bifpoints[1,]; #inv of C1 in C2-R2+R1
Eq_K1_IIAXIB_BP2 <- output8IIAXIBc$bifpoints[2,];
Eq_K1_IIAXIB_BPE <- output8IIAXIBb$bifpoints[1,];
Eq_K1_IIAXIB_LP <- output8IIAXIBc$bifpoints[1,];
Eq_K1_IIAXIB_LPb <- output3.2$bifpoints[2,];
Eq_K1_IIAXIB_LPc <- output8IIAXIBg$bifpoints[2,];

### Codim ========================
### BP 1
cmd = output14IIAXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IIAXIB$bifpoints[,c(1:3,7)], 1.2), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "2", "envZE", "3", "envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output14IIAXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IIAXIB$bifpoints[,c(1:3,7)], 1.2), 0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "2", "envZE", "3", "envZE", "4", "popBP", "1"), clean=FALSE)
output14IIAXIB<-rbind(output14IIAXIBa$curvepoints[dim(output14IIAXIBa$curvepoints)[1]:1,],
                      output14IIAXIBb$curvepoints[2:dim(output14IIAXIBb$curvepoints)[1],]);
### BP 2
cmd = output15IIAXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IIAXIBc$bifpoints[2,c(1:6,8)], 1.2), -0.1, c(11, 0, 9.5E-1, 20, 0.1, 4.5), NULL, c("popBP", "0"), clean=FALSE)
cmd = output15IIAXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IIAXIBc$bifpoints[2,c(1:6,8)], 1.2), 0.1, c(11, 0, 9.5E-1, 20, 0.1, 4.5), NULL, c("popBP", "0"), clean=FALSE)
output15IIAXIB<-rbind(output15IIAXIBa$curvepoints[dim(output15IIAXIBa$curvepoints)[1]:1,],
                      output15IIAXIBb$curvepoints[2:dim(output15IIAXIBb$curvepoints)[1],]);
### BPE
cmd = output16IIAXIBa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output8IIAXIBb$bifpoints[1,c(1:3,7,8)], 1.2), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)
cmd = output16IIAXIBb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output8IIAXIBb$bifpoints[1,c(1:3,7,8)], 1.2), 0.1, c(11, 0, 9.5E-01, 20, -0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)
output16IIAXIB<-rbind(output16IIAXIBa$curvepoints[dim(output16IIAXIBa$curvepoints)[1]:1,],
                      output16IIAXIBb$curvepoints[2:dim(output16IIAXIBb$curvepoints)[1],]);
### LP
cmd = output17IIAXIBa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IIAXIBc$bifpoints[1,c(1:8)],1.2), -0.1, c(11, 0, 9.5E-01, 20, -0.1, 4.5), NULL, NULL, clean=FALSE)
cmd = output17IIAXIBb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IIAXIBc$bifpoints[1,c(1:8)],1.2), 0.1, c(11, 0, 9.5E-01, 20, -0.1, 4.5), NULL, NULL, clean=FALSE)
output17IIAXIB<-rbind(output17IIAXIBa$curvepoints[dim(output17IIAXIBa$curvepoints)[1]:1,],
                      output17IIAXIBb$curvepoints[2:dim(output17IIAXIBb$curvepoints)[1],]);
cmd = output17IIAXIBc <- PSPMequi("Tritrophicmod5_fin", "LP", c(output3.2$bifpoints[2,c(1:4,6,8)], 2.5), -0.1, c(11, 0, 9.5E-01, 20, -0.1, 4.5), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output17IIAXIBd <- PSPMequi("Tritrophicmod5_fin", "LP", c(output3.2$bifpoints[2,c(1:4,6,8)], 2.5), 1.2, c(11, 0, 0.5E-01, 20, -0.1, 4.5), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
output17IIAXIB2<-rbind(output17IIAXIBc$curvepoints[dim(output17IIAXIBc$curvepoints)[1]:1,],
                       output17IIAXIBd$curvepoints[2:dim(output17IIAXIBd$curvepoints)[1],]);
cmd = output17IIAXIBe <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8ICXIBg$bifpoints[2,c(1:8)], 3.5), -0.5, c(11, 0, 9.5E-01, 20, -0.1, 4.5), NULL, NULL, clean=FALSE)
cmd = output17IIAXIBf <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8ICXIBg$bifpoints[2,c(1:8)], 3.5), 1.2, c(11, 0, 0.5E-01, 20, -0.1, 4.5), NULL, NULL, clean=FALSE)
output17IIAXIB3<-rbind(output17IIAXIBe$curvepoints[dim(output17IIAXIBe$curvepoints)[1]:1,],
                       output17IIAXIBf$curvepoints[2:dim(output17IIAXIBf$curvepoints)[1],]);
################## IIB. K2=8e-5 ========================
### EQ ========================
cmd = output8IIBXIB <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[52,2], 1E-8, output2$curvepoints[52,c(5,7)]),
                                0.1, c(11, 0, 4E-4), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = output8IIBXIBb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IIBXIB$bifpoints[1,c(1:3,6:8)]),
                                 0.1, c(11, 0, 4E-4), NULL, c("envZE", "2","envZE", "4"), clean=FALSE)
cmd = output8IIBXIBc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IIBXIBb$bifpoints[,c(1:8)]),
                                 -0.5, c(11, 0, 4E-4), NULL, NULL, clean=FALSE)
### with Beta = 2.5, detection of LP2
cmd = output8IIBXICd <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output3.2$bifpoints[1,c(1:8)]),
                                 -2, c(11, 0, 4E-3), NULL, NULL, clean=FALSE)

Eq_K1_IIBXIB_BP1 <- output8IIBXIB$bifpoints[1,]; #inv of C1 in C2-R2+R1
Eq_K1_IIBXIB_BP2 <- output8IIBXIBc$bifpoints[2,];
Eq_K1_IIBXIB_BPE <- output8IIBXIBb$bifpoints[1,];
Eq_K1_IIBXIB_LP <- output8IIBXIBc$bifpoints[1,];
Eq_K1_IIBXIB_LPb <- output3.2$bifpoints[2,];
### Codim ========================
### BP 1 
cmd = output14IIBXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IIBXIB$bifpoints[,c(1:3,7)], 1.2), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "2", "envZE", "3", "envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output14IIBXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IIBXIB$bifpoints[,c(1:3,7)], 1.2), 0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "2", "envZE", "3", "envZE", "4", "popBP", "1"), clean=FALSE)
output14IIBXIB<-rbind(output14IIBXIBa$curvepoints[dim(output14IIBXIBa$curvepoints)[1]:1,],
                      output14IIBXIBb$curvepoints[2:dim(output14IIBXIBb$curvepoints)[1],]);
### BP 2
cmd = output15IIBXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IIBXIBc$bifpoints[2,c(1:6,8)], 1.2), -0.1, c(11, 0, 9.5E-1, 20, 0.1, 4.5), NULL, c("popBP", "0"), clean=FALSE)
cmd = output15IIBXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IIBXIBc$bifpoints[2,c(1:6,8)], 1.2), 0.1, c(11, 0, 9.5E-1, 20, 0.1, 4.5), NULL, c("popBP", "0"), clean=FALSE)
output15IIBXIB<-rbind(output15IIBXIBa$curvepoints[dim(output15IIBXIBa$curvepoints)[1]:1,],
                      output15IIBXIBb$curvepoints[2:dim(output15IIBXIBb$curvepoints)[1],]);
### BPE
cmd = output16IIBXIBa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output8IIBXIBb$bifpoints[1,c(1:3,7,8)], 1.2), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)
cmd = output16IIBXIBb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output8IIBXIBb$bifpoints[1,c(1:3,7,8)], 1.2), 0.1, c(11, 0, 9.5E-01, 20, -0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)
output16IIBXIB<-rbind(output16IIBXIBa$curvepoints[dim(output16IIBXIBa$curvepoints)[1]:1,],
                      output16IIBXIBb$curvepoints[2:dim(output16IIBXIBb$curvepoints)[1],]);
### LP
cmd = output17IIBXIBa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IIBXIBc$bifpoints[1,c(1:8)],1.2), -0.1, c(11, 0, 9.5E-01, 20, -0.1, 4.5), NULL, NULL, clean=FALSE)
cmd = output17IIBXIBb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IIBXIBc$bifpoints[1,c(1:8)],1.2), 0.1, c(11, 0, 9.5E-01, 20, -0.1, 4.5), NULL, NULL, clean=FALSE)
output17IIBXIB<-rbind(output17IIBXIBa$curvepoints[dim(output17IIBXIBa$curvepoints)[1]:1,],
                      output17IIBXIBb$curvepoints[2:dim(output17IIBXIBb$curvepoints)[1],]);
cmd = output17IIBXIBc <- PSPMequi("Tritrophicmod5_fin", "LP", c(output3.2$bifpoints[2,c(1:4,6,8)], 2.5), -0.1, c(11, 0, 9.5E-01, 20, -0.1, 4.5), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output17IIBXIBd <- PSPMequi("Tritrophicmod5_fin", "LP", c(output3.2$bifpoints[2,c(1:4,6,8)], 2.5), 1.2, c(11, 0, 0.5E-01, 20, -0.1, 4.5), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
output17IIBXIB2<-rbind(output17IIBXIBc$curvepoints[dim(output17IIBXIBc$curvepoints)[1]:1,],
                       output17IIBXIBd$curvepoints[2:dim(output17IIBXIBd$curvepoints)[1],]);
################## IIC. K2=1e-4 ========================
### EQ ========================
cmd = output8IICXIB <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[62,2], 1E-8, output2$curvepoints[62,c(5,7)]),
                                0.1, c(11, 0, 4E-4), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = output8IICXIBb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IICXIB$bifpoints[1,c(1:3,5:8)]),
                                 0.1, c(11, 0, 4E-4), NULL, c("envZE", "2"), clean=FALSE)
### with Beta = 1.0 to locate Bp1#3
cmd = output8IICXIBc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1.3e-04,output3.2$bifpoints[3,c(2:6,7,8)]), -0.1, c(11, 0, 9.5E-01), NULL, NULL, clean=FALSE)
### with Beta = 0.5 to locate LP & Bp1#2
cmd = output8IICXIBd <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IICXIBb$bifpoints[2,c(1:8)]),
                                 -0.5, c(11, 0, 4E-4), NULL, NULL, clean=FALSE)

### with Beta = 2 to locate LP#2 to define the tri-stability 2/4/5
cmd = output8IICXIBe <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output3$curvepoints[277,2], 1E-8, output3$curvepoints[277,c(4,5,7)]),
                                 0.1, c(11, 0, 2E-3), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)
cmd = output8IICXIBf <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IICXIBe$bifpoints[1,c(1:8)]),
                                 0.1, c(11, 0, 2E-3), NULL, NULL, clean=FALSE)
cmd = output8IICXIBg <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IICXIBf$bifpoints[1,c(1:4,6,8)]),
                                 -0.1, c(11, 0, 2E-3), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

### with Beta = 0.3 to locate LP#3
cmd = output8IICXIBh <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IICXIBb$bifpoints[,c(1:8)]),
                                 0.1, c(11, 0, 4E-4), NULL, NULL, clean=FALSE)
cmd = output8IICXIBi <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IICXIBh$bifpoints[2,c(1:5,7)]),
                                 -0.1, c(11, 0, 8E-4), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

Eq_K1_IICXIB_BP1 <- output8IICXIB$bifpoints[1,]; #inv of C1 in C2-R2+R1
Eq_K1_IICXIB_BP1b <- output8IICXIBb$bifpoints[2,]; #inv of C1 on unstable branch
Eq_K1_IICXIB_BP1c <- output8IICXIBc$bifpoints[1,]; #inv of C1 in P-C2-R2+R1
Eq_K1_IICXIB_BP2 <- output3.2$bifpoints[3,]; #with beta=1.0
Eq_K1_IICXIB_BPE <- output8IICXIBb$bifpoints[1,];
Eq_K1_IICXIB_LP <- output8IICXIBd$bifpoints[1,];
Eq_K1_IICXIB_LPb <- output8IICXIBg$bifpoints[1,];
Eq_K1_IICXIB_LPc <- output8IICXIBh$bifpoints[1,];
### Codim ========================
### BP 1 
cmd = output14IICXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IICXIB$bifpoints[,c(1:3,7)], 1.0), -0.3, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "2", "envZE", "3", "envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output14IICXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IICXIB$bifpoints[,c(1:3,7)], 1.0), 0.3, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "2", "envZE", "3", "envZE", "4", "popBP", "1"), clean=FALSE)
output14IICXIB<-rbind(output14IICXIBa$curvepoints[dim(output14IICXIBa$curvepoints)[1]:1,],
                      output14IICXIBb$curvepoints[2:dim(output14IICXIBb$curvepoints)[1],]);

cmd = output14IICXIBc <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IICXIBd$bifpoints[2,c(1:5,7)], 0.5), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output14IICXIBd <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IICXIBd$bifpoints[2,c(1:5,7)], 0.5), 0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
output14IICXIB2<-rbind(output14IICXIBc$curvepoints[dim(output14IICXIBc$curvepoints)[1]:1,],
                       output14IICXIBd$curvepoints[2:dim(output14IICXIBd$curvepoints)[1],]);

cmd = output14IICXIBe <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IICXIBc$bifpoints[1,c(1:5,7)], 1), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output14IICXIBf <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IICXIBc$bifpoints[1,c(1:5,7)], 1), 0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
output14IICXIB3<-rbind(output14IICXIBe$curvepoints[dim(output14IICXIBe$curvepoints)[1]:1,],
                       output14IICXIBf$curvepoints[2:dim(output14IICXIBf$curvepoints)[1],]);
### BP 2
cmd = output15IICXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output3.2$bifpoints[3,c(1:6,8)], 1.0), -0.1, c(11, 0, 9.5E-1, 20, 0.1, 4.5), NULL, c("popBP", "0"), clean=FALSE)
cmd = output15IICXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output3.2$bifpoints[3,c(1:6,8)], 1.0), 0.1, c(11, 0, 9.5E-1, 20, 0.1, 4.5), NULL, c("popBP", "0"), clean=FALSE)
output15IICXIB<-rbind(output15IICXIBa$curvepoints[dim(output15IICXIBa$curvepoints)[1]:1,],
                      output15IICXIBb$curvepoints[2:dim(output15IICXIBb$curvepoints)[1],]);
### BPE
cmd = output16IICXIBa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output8IICXIBb$bifpoints[2,c(1:3,7,8)], 1.0), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)
cmd = output16IICXIBb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output8IICXIBb$bifpoints[2,c(1:3,7,8)], 1.0), 0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)
output16IICXIB<-rbind(output16IICXIBa$curvepoints[dim(output16IICXIBa$curvepoints)[1]:1,],
                      output16IICXIBb$curvepoints[2:dim(output16IICXIBb$curvepoints)[1],]);
### LP 
cmd = output17IICXIBa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IICXIBd$bifpoints[1,c(1:8)], 0.5), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, NULL, clean=FALSE)
cmd = output17IICXIBb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IICXIBd$bifpoints[1,c(1:8)], 0.5), 0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, NULL, clean=FALSE)
output17IICXIB<-rbind(output17IICXIBa$curvepoints[dim(output17IICXIBa$curvepoints)[1]:1,],
                      output17IICXIBb$curvepoints[2:dim(output17IICXIBb$curvepoints)[1],]);

cmd = output17IICXIBc <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IICXIBg$bifpoints[1,c(1:4,6,8)], 2), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output17IICXIBd <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IICXIBg$bifpoints[1,c(1:4,6,8)], 2), 0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
output17IICXIB2<-rbind(output17IICXIBc$curvepoints[dim(output17IICXIBc$curvepoints)[1]:1,],
                       output17IICXIBd$curvepoints[2:dim(output17IICXIBd$curvepoints)[1],]);

cmd = output17IICXIBe <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IICXIBh$bifpoints[1,c(1:8)], 0.3), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, NULL, clean=FALSE)
cmd = output17IICXIBf <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IICXIBh$bifpoints[1,c(1:8)], 0.3), 1.3, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, NULL, clean=FALSE)
output17IICXIB3<-rbind(output17IICXIBe$curvepoints[dim(output17IICXIBe$curvepoints)[1]:1,],
                       output17IICXIBf$curvepoints[2:dim(output17IICXIBf$curvepoints)[1],]);

################## IIIA. K2=2e-4 ========================
### EQ ========================
cmd = output8IIIAXIB <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[78,2], 1E-8, output2$curvepoints[78,c(5,7)]),
                                 0.1, c(11, 0, 4E-4), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = output8IIIAXIBb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IIIAXIB$bifpoints[1,c(1:3,6:8)]),
                                  0.1, c(11, 0, 4E-4), NULL, c("envZE", "2","envZE", "4"), clean=FALSE)

### with Beta = 1.0 to locate Bp1#3
cmd = output8IIIAXIBc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1.3e-04,output3.2$bifpoints[3,c(2:6,7,8)]), -0.1, c(11, 0, 9.5E-01), NULL, NULL, clean=FALSE)

### with Beta = 0.3 to locate LP & Bp1#2
cmd = output8IIIAXIBd <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IIIAXIBb$bifpoints[2,c(1:8)]),
                                  -0.5, c(11, 0, 4E-4), NULL, NULL, clean=FALSE)
### with Beta = 0.5 to locate Bp1#2
cmd = output8IIIAXIBe <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IIIAXIBb$bifpoints[1,c(1:8)]),
                                  -0.5, c(11, 0, 4E-4), NULL, NULL, clean=FALSE)
### with Beta = 3 to locate LP#2
cmd = output8IIIAXIBf <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output3$curvepoints[371,2], 1E-8, output3$curvepoints[371,c(4,5,7)]),
                                  0.1, c(11, 0, 2E-3), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)
cmd = output8IIIAXIBg <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IIIAXIBf$bifpoints[1,c(1:8)]),
                                  0.1, c(11, 0, 2E-3), NULL, NULL, clean=FALSE)
cmd = output8IIIAXIBh <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IIIAXIBg$bifpoints[1,c(1:4,6,8)]),
                                  -0.1, c(11, 0, 2E-3), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K1_IIIAXIB_BP1 <- output8IIIAXIB$bifpoints[1,]; #inv of C1 in C2-R2+R1
Eq_K1_IIIAXIB_BP1b <- output8IIIAXIBe$bifpoints[1,]; #inv of C1 on unstable branch
Eq_K1_IIIAXIB_BP1c <- output8IIIAXIBc$bifpoints[1,]; #inv of C1 in P-C2-R2+R1
Eq_K1_IIIAXIB_BP2 <- output8IIIAXIBg$biftypes; #with beta=1.0
Eq_K1_IIIAXIB_BPE <- output8IIIAXIBb$bifpoints[1,];
Eq_K1_IIIAXIB_LP <- output8IIIAXIBd$bifpoints[1,];
Eq_K1_IIIAXIB_LPb <- output8IIIAXIBh$bifpoints[1,];
### Codim ========================
### BP 1 
cmd = output14IIIAXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IIIAXIB$bifpoints[,c(1:3,7)], 1.0), -0.6, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "2", "envZE", "3", "envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output14IIIAXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IIIAXIB$bifpoints[,c(1:3,7)], 1.0), 0.7, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "2", "envZE", "3", "envZE", "4", "popBP", "1"), clean=FALSE)
output14IIIAXIB<-rbind(output14IIIAXIBa$curvepoints[dim(output14IIIAXIBa$curvepoints)[1]:1,],
                       output14IIIAXIBb$curvepoints[2:dim(output14IIIAXIBb$curvepoints)[1],]);

cmd = output14IIIAXIBc <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IIIAXIBe$bifpoints[1,c(1:5,7)], 0.5), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output14IIIAXIBd <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IIIAXIBe$bifpoints[1,c(1:5,7)], 0.5), 0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
output14IIIAXIB2<-rbind(output14IIIAXIBc$curvepoints[dim(output14IIIAXIBc$curvepoints)[1]:1,],
                        output14IIIAXIBd$curvepoints[2:dim(output14IIIAXIBd$curvepoints)[1],]);

cmd = output14IIIAXIBe <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IIIAXIBc$bifpoints[1,c(1:5,7)], 1), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output14IIIAXIBf <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IIIAXIBc$bifpoints[1,c(1:5,7)], 1), 0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
output14IIIAXIB3<-rbind(output14IIIAXIBe$curvepoints[dim(output14IIIAXIBe$curvepoints)[1]:1,],
                        output14IIIAXIBf$curvepoints[2:dim(output14IIIAXIBf$curvepoints)[1],]);

### BP 2
cmd = output15IIIAXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IIIAXIBg$bifpoints[1,c(1:6,8)], 3.0), -0.1, c(11, 0, 9.5E-1, 20, 0.1, 4.5), NULL, c("popBP", "0"), clean=FALSE)
cmd = output15IIIAXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IIIAXIBg$bifpoints[1,c(1:6,8)], 3.0), 0.1, c(11, 0, 9.5E-1, 20, 0.1, 4.5), NULL, c("popBP", "0"), clean=FALSE)
output15IIIAXIB<-rbind(output15IIIAXIBa$curvepoints[dim(output15IIIAXIBa$curvepoints)[1]:1,],
                       output15IIIAXIBb$curvepoints[2:dim(output15IIIAXIBb$curvepoints)[1],]);
### BPE
cmd = output16IICXIBa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output8IICXIBb$bifpoints[2,c(1:3,7,8)], 1.0), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)

cmd = output16IIIAXIBa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output8IIIAXIBb$bifpoints[1,c(1:3,7,8)], 1.0), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)
cmd = output16IIIAXIBb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output8IIIAXIBb$bifpoints[1,c(1:3,7,8)], 1.0), 0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)
output16IIIAXIB<-rbind(output16IIIAXIBa$curvepoints[dim(output16IIIAXIBa$curvepoints)[1]:1,],
                       output16IIIAXIBb$curvepoints[2:dim(output16IIIAXIBb$curvepoints)[1],]);

### LP with maybe beta=0.3
cmd = output17IIIAXIBa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IIIAXIBd$bifpoints[1,c(1:8)], 0.3), -0.1, c(11, 0, 9.5E-01, 20, -0.1, 4.5), NULL, NULL, clean=FALSE)
cmd = output17IIIAXIBb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IIIAXIBd$bifpoints[1,c(1:8)], 0.3), 0.1, c(11, 0, 9.5E-01, 20, -0.1, 4.5), NULL, NULL, clean=FALSE)
output17IIIAXIB<-rbind(output17IIIAXIBa$curvepoints[dim(output17IIIAXIBa$curvepoints)[1]:1,],
                       output17IIIAXIBb$curvepoints[2:dim(output17IIIAXIBb$curvepoints)[1],]);
cmd = output17IIIAXIBc <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IIIAXIBh$bifpoints[1,c(1:4,6,8)], 3), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output17IIIAXIBd <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IIIAXIBh$bifpoints[1,c(1:4,6,8)], 3), 0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
output17IIIAXIB2<-rbind(output17IIIAXIBc$curvepoints[dim(output17IIIAXIBc$curvepoints)[1]:1,],
                        output17IIIAXIBd$curvepoints[2:dim(output17IIIAXIBd$curvepoints)[1],]);

################## IIIB. K2=3e-4 ========================
### EQ ========================
cmd = output7IIIBXIB <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1e-8, output3$curvepoints[426,2], 1e-8,output3$curvepoints[426,c(4,5,7)]), 0.1, c(11, 0, 9.5E-01), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)
cmd = output7IIIBXIBb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIIBXIB$bifpoints[1,c(1:8)]), 0.1, c(11, 0, 9.5E-01), NULL, NULL, clean=FALSE)
cmd = output7IIIBXIBc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIIBXIBb$bifpoints[1,c(1:4,6,8)]), -0.1, c(11, 0, 9.5E-01), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K1_IIIBXIB_BP1 <- output7IIIBXIB$bifpoints[1,];
Eq_K1_IIIBXIB_BP2 <- output7IIIBXIBb$bifpoints[1,];
Eq_K1_IIIBXIB_LP <- output7IIIBXIBc$bifpoints[2,];
### Codim ========================
### BP 1 
cmd = output14IIIBXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIIBXIB$bifpoints[2,c(1:5,7)], 1.5), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output14IIIBXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIIBXIB$bifpoints[2,c(1:5,7)], 1.5), 0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
output14IIIBXIB<-rbind(output14IIIBXIBa$curvepoints[dim(output14IIIBXIBa$curvepoints)[1]:1,],
                       output14IIIBXIBb$curvepoints[2:dim(output14IIIBXIBb$curvepoints)[1],]);
### BP 2 
cmd = output15IIIBXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIIBXIB$bifpoints[1,c(1:5,8)], 1.5), -0.1, c(11, 0, 9.5E-1, 20, 0.1, 4.5), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
cmd = output15IIIBXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIIBXIB$bifpoints[1,c(1:5,8)], 1.5), 0.1, c(11, 0, 9.5E-1, 20, 0.1, 4.5), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
output15IIIBXIB<-rbind(output15IIIBXIBa$curvepoints[dim(output15IIIBXIBa$curvepoints)[1]:1,],
                       output15IIIBXIBb$curvepoints[2:dim(output15IIIBXIBb$curvepoints)[1],]);
### LP
cmd = output17IIIBXIBa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIIBXIBc$bifpoints[2,c(1:4,6,8)], 3.5), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output17IIIBXIBb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIIBXIBc$bifpoints[2,c(1:4,6,8)], 3.5), 1.3, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
output17IIIBXIB<-rbind(output17IIIBXIBa$curvepoints[dim(output17IIIBXIBa$curvepoints)[1]:1,],
                       output17IIIBXIBb$curvepoints[2:dim(output17IIIBXIBb$curvepoints)[1],]);
################## IIIC. K2=4e-4 ========================
### EQ ========================
cmd = output7IIICXIB <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1e-8,output3$curvepoints[426,2],1e-8, output3$curvepoints[426,c(4,5,7)]), 0.1, c(11, 0, 9.5E-01), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)
cmd = output7IIICXIBb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIICXIB$bifpoints[1,c(1:8)]), 0.1, c(11, 0, 9.5E-01), NULL, NULL, clean=FALSE)
cmd = output7IIICXIBc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIICXIBb$bifpoints[1,c(1:4,6,8)]), -0.1, c(11, 0, 9.5E-01), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K1_IIICXIB_BP1 <- output7IIICXIB$bifpoints[1,];
Eq_K1_IIICXIB_BP2 <- output7IIICXIB$bifpoints[1,];
Eq_K1_IIICXIB_LP <- output7IIICXIBc$bifpoints[2,];
### Codim ========================
### BP 1 
cmd = output14IIICXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIICXIB$bifpoints[2,c(1:5,7)], 1.5), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output14IIICXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIICXIB$bifpoints[2,c(1:5,7)], 1.5), 0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
output14IIICXIB<-rbind(output14IIICXIBa$curvepoints[dim(output14IIICXIBa$curvepoints)[1]:1,],
                       output14IIICXIBb$curvepoints[2:dim(output14IIICXIBb$curvepoints)[1],]);
### BP 2
cmd = output15IIICXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIICXIB$bifpoints[1,c(1:5,8)], 1.5), -0.1, c(11, 0, 9.5E-1, 20, 0.1, 4.5), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
cmd = output15IIICXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIICXIB$bifpoints[1,c(1:5,8)], 1.5), 0.1, c(11, 0, 9.5E-1, 20, 0.1, 4.5), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
output15IIICXIB<-rbind(output15IIICXIBa$curvepoints[dim(output15IIICXIBa$curvepoints)[1]:1,],
                       output15IIICXIBb$curvepoints[2:dim(output15IIICXIBb$curvepoints)[1],]);
### LP
cmd = output17IIICXIBa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIICXIBc$bifpoints[2,c(1:4,6,8)], 3.5), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output17IIICXIBb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIICXIBc$bifpoints[2,c(1:4,6,8)], 3.5), 1.3, c(11, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
output17IIICXIB<-rbind(output17IIICXIBa$curvepoints[dim(output17IIICXIBa$curvepoints)[1]:1,],
                       output17IIICXIBb$curvepoints[2:dim(output17IIICXIBb$curvepoints)[1],]);