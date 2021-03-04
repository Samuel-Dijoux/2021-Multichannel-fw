###################################################################################### 
#                -> Bifurcation analyses over K2 and Beta gradients for varying K1 
#
# Codimension bifurcation tracking selected BP
# as a function of two parameters : K1 and K2+Beta
# Representation of the coexistence phase of both consumers: P-C1C2-R1R2 ~ K2+Beta            
######################################################################################
library(PSPManalysis)
if (!exists("PSPMsrc.fullpath")) source("../PSPMworkshop/PSPManalysis/PSPManalysis.r")

################## IA. K1=1e-5 ========================
### EQ ========================
cmd = output7IAXIC <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output2.2$curvepoints[8,c(3,6,8)]),
                               0.5, c(10, 0, 4E-4), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7IAXICb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IAXIC$bifpoints[1,c(1:3,6:8)]),
                                0.1, c(10, 0, 9E-4), NULL, c("envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7IAXICc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IAXICb$bifpoints[,c(1:8)]),
                                -0.1, c(10, 0, 6E-4), NULL, NULL, clean=FALSE)
# with Beta = 0.5 to detect LPb, After extinction of C2
cmd = output7IAXICd <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IAXICc$bifpoints[2,c(1:5,7)]), -0.1, c(10, 0, 2E-3), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

Eq_K2_IAXIC_BP1 <- output7IAXICc$bifpoints[2,]; 
Eq_K2_IAXIC_BP2 <- output7IAXIC$bifpoints[1,]; #inv of C2 in C1-R1+R2
Eq_K2_IAXIC_BPE <- output7IAXICb$bifpoints[1,];
Eq_K2_IAXIC_LP <- output7IAXICc$bifpoints[1,];
Eq_K2_IAXIC_LPb <- output7IAXICd$bifpoints[2,];
### Codim ========================
### BP 1
cmd = output18IAXICa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IAXICc$bifpoints[2,c(1:6,8)], 1.5), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popBP", "1"), clean=FALSE);
cmd = output18IAXICb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IAXICc$bifpoints[2,c(1:6,8)], 1.5), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popBP", "1"), clean=FALSE);
output18IAXIC<-rbind(output18IAXICa$curvepoints[dim(output18IAXICa$curvepoints)[1]:1,],
                     output18IAXICb$curvepoints[2:dim(output18IAXICb$curvepoints)[1],]);

### BPE
cmd = output20IAXICa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IAXICb$bifpoints[1,c(1:3,7,8)], 1.5), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE);
cmd = output20IAXICb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IAXICb$bifpoints[1,c(1:3,7,8)], 1.5), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE);
output20IAXIC<-rbind(output20IAXICa$curvepoints[dim(output20IAXICa$curvepoints)[1]:1,],
                     output20IAXICb$curvepoints[2:dim(output20IAXICb$curvepoints)[1],]);

### LP
cmd = output21IAXICa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IAXICc$bifpoints[1,c(1:8)], 1.5), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, NULL, clean=FALSE);
cmd = output21IAXICb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IAXICc$bifpoints[1,c(1:8)], 1.5), 0.1, c(10, 0, 9.5E-04, 20, 0.5, 4.5), NULL, NULL, clean=FALSE);
output21IAXIC<-rbind(output21IAXICb$curvepoints[dim(output21IAXICb$curvepoints)[1]:1,],
                     output21IAXICa$curvepoints[2:dim(output21IAXICa$curvepoints)[1],]);

cmd = output21IAXICc <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IAXICd$bifpoints[2,c(1:5,7)], 0.5), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE);
cmd = output21IAXICd <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IAXICd$bifpoints[2,c(1:5,7)], 0.5), 0.1, c(10, 0, 9.5E-04, 20, 0.5, 4.5), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE);
output21IAXIC2<-rbind(output21IAXICc$curvepoints[dim(output21IAXICc$curvepoints)[1]:1,],
                      output21IAXICd$curvepoints[2:dim(output21IAXICd$curvepoints)[1],]);

################## IB. K1=2e-5 ========================
### EQ ========================
cmd = output7IBXIB <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output2.2$curvepoints[22,c(3,6,8)]),
                               0.5, c(10, 0, 4E-4), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7IBXIBb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IBXIB$bifpoints[1,c(1:3,6:8)]),
                                0.1, c(10, 0, 9E-4), NULL, c("envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7IBXIBc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IBXIBb$bifpoints[,c(1:8)]),
                                -0.1, c(10, 0, 6E-4), NULL, NULL, clean=FALSE)
cmd = output7IBXIBd <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IBXIBc$bifpoints[3,c(1:5,7)]),
                                -0.1, c(10, 0, 6E-4), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

Eq_K2_IBXIB_BP1 <- output7IBXIBc$bifpoints[3,]; 
Eq_K2_IBXIB_BP2 <- output7IBXIB$bifpoints[1,]; #inv of C2 in C1-R1+R2
Eq_K2_IBXIB_BPE <- output7IBXIBb$bifpoints[1,];
Eq_K2_IBXIB_LP <- output7IBXIBc$bifpoints[2,];
Eq_K2_IBXIB_LPb <- output7IBXIBd$bifpoints[1,];
### Codim ========================
### BP 1
cmd = output18IBXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IBXIBc$bifpoints[3,c(1:6,8)], 1), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popBP", "1"), clean=FALSE);
cmd = output18IBXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IBXIBc$bifpoints[3,c(1:6,8)], 1), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popBP", "1"), clean=FALSE);
output18IBXIB<-rbind(output18IBXIBa$curvepoints[dim(output18IBXIBa$curvepoints)[1]:1,],
                     output18IBXIBb$curvepoints[2:dim(output18IBXIBb$curvepoints)[1],]);

### BPE
cmd = output20IBXIBa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IBXIBb$bifpoints[1,c(1:3,7,8)], 1), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE);
cmd = output20IBXIBb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IBXIBb$bifpoints[1,c(1:3,7,8)], 1), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE);
output20IBXIB<-rbind(output20IBXIBa$curvepoints[dim(output20IBXIBa$curvepoints)[1]:1,],
                     output20IBXIBb$curvepoints[2:dim(output20IBXIBb$curvepoints)[1],]);

### LP
cmd = output21IBXIBa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IBXIBc$bifpoints[2,c(1:8)], 0.5), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, NULL, clean=FALSE);
cmd = output21IBXIBb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IBXIBc$bifpoints[2,c(1:8)], 0.5), 1.3, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, NULL, clean=FALSE);
output21IBXIB<-rbind(output21IBXIBa$curvepoints[dim(output21IBXIBa$curvepoints)[1]:1,],
                     output21IBXIBb$curvepoints[2:dim(output21IBXIBb$curvepoints)[1],]);

cmd = output21IBXIBc <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IBXIBd$bifpoints[1,c(1:5,7)], 0.5), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE);
cmd = output21IBXIBd <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IBXIBd$bifpoints[1,c(1:5,7)], 0.5), 1.3, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE);
output21IBXIB2<-rbind(output21IBXIBc$curvepoints[dim(output21IBXIBc$curvepoints)[1]:1,],
                      output21IBXIBd$curvepoints[2:dim(output21IBXIBd$curvepoints)[1],]);
################## IC. K1=3e-5 ========================
### EQ ========================
cmd = output7ICXIB <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output2.2$curvepoints[27,c(3,6,8)]),
                               0.5, c(10, 0, 4E-4), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7ICXIBb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7ICXIB$bifpoints[1,c(1:3,6:8)]),
                                0.1, c(10, 0, 9E-4), NULL, c("envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7ICXIBc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7ICXIBb$bifpoints[,c(1:8)]),
                                -0.1, c(10, 0, 6E-4), NULL, NULL, clean=FALSE)

### Starting from P-C1-R1+R2 with Beta=0.5 to detect LP#2 
cmd = output7ICXIBd <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output3.2$curvepoints[87,c(3,4,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output7ICXIBe <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7ICXIBd$bifpoints[1,c(1:8)]), 0.1, c(10, 0, 2E-3), NULL, NULL, clean=FALSE)
cmd = output7ICXIBf <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7ICXIBe$bifpoints[1,c(1:5,7)]), -0.1, c(10, 0, 2E-3), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

Eq_K2_ICXIB_BP1 <- output7ICXIBc$bifpoints[3,]; 
Eq_K2_ICXIB_BP2 <- output7ICXIB$bifpoints[1,]; #inv of C2 in C1-R1+R2
Eq_K2_ICXIB_BPE <- output7ICXIBb$bifpoints[1,];
Eq_K2_ICXIB_LP <- output7ICXIBc$bifpoints[2,];
Eq_K2_ICXIB_LPb <- output7ICXIBf$bifpoints[2,];
### Codim ========================
### BP 1
cmd = output18ICXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7ICXIBc$bifpoints[3,c(1:6,8)], 1), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popBP", "1"), clean=FALSE);
cmd = output18ICXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7ICXIBc$bifpoints[3,c(1:6,8)], 1), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popBP", "1"), clean=FALSE);
output18ICXIB<-rbind(output18ICXIBa$curvepoints[dim(output18ICXIBa$curvepoints)[1]:1,],
                     output18ICXIBb$curvepoints[2:dim(output18ICXIBb$curvepoints)[1],]);

### BPE
cmd = output20ICXIBa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7ICXIBb$bifpoints[1,c(1:3,7,8)], 1), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE);
cmd = output20ICXIBb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7ICXIBb$bifpoints[1,c(1:3,7,8)], 1), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE);
output20ICXIB<-rbind(output20ICXIBa$curvepoints[dim(output20ICXIBa$curvepoints)[1]:1,],
                     output20ICXIBb$curvepoints[2:dim(output20ICXIBb$curvepoints)[1],]);

### LP
cmd = output21ICXIBa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7ICXIBc$bifpoints[2,c(1:8)], 1), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, NULL, clean=FALSE);
cmd = output21ICXIBb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7ICXIBc$bifpoints[2,c(1:8)], 1), 0.1, c(10, 0, 9.5E-04, 20, 0.1, 4.5), NULL, NULL, clean=FALSE);
output21ICXIB<-rbind(output21ICXIBb$curvepoints[dim(output21ICXIBb$curvepoints)[1]:1,],
                     output21ICXIBa$curvepoints[2:dim(output21ICXIBa$curvepoints)[1],]);

cmd = output21ICXIBc <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7ICXIBf$bifpoints[2,c(1:5,7)], 0.5), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE);
cmd = output21ICXIBd <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7ICXIBf$bifpoints[2,c(1:5,7)], 0.5), 0.1, c(10, 0, 9.5E-04, 20, 0.1, 4.5), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE);
output21ICXIB2<-rbind(output21ICXIBc$curvepoints[dim(output21ICXIBc$curvepoints)[1]:1,],
                      output21ICXIBd$curvepoints[2:dim(output21ICXIBd$curvepoints)[1],]);
################## IIA. K1=5e-5 ========================
### EQ ========================
cmd = output7IIAXIB <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output2.2$curvepoints[37,c(3,6,8)]),
                                0.5, c(10, 0, 4E-4), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7IIAXIBb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIAXIB$bifpoints[1,c(1:3,6:8)]),
                                 0.1, c(10, 0, 9E-4), NULL, c("envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7IIAXIBc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIAXIBb$bifpoints[,c(1:8)]),
                                 -0.1, c(10, 0, 6E-4), NULL, NULL, clean=FALSE)

### Starting from P-C1-R1+R2 with Beta=0.5 to detect LP#2 
cmd = output7IIAXIBd <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output3.2$curvepoints[113,c(3,4,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output7IIAXIBe <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIAXIBd$bifpoints[1,c(1:8)]), 0.1, c(10, 0, 2E-3), NULL, NULL, clean=FALSE)
cmd = output7IIAXIBf <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIAXIBe$bifpoints[1,c(1:5,7)]), -0.1, c(10, 0, 2E-3), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

Eq_K2_IIAXIB_BP1 <- output7IIAXIBc$bifpoints[3,]; 
Eq_K2_IIAXIB_BP2 <- output7IIAXIB$bifpoints[1,]; #inv of C2 in C1-R1+R2
Eq_K2_IIAXIB_BPE <- output7IIAXIBb$bifpoints[1,];
Eq_K2_IIAXIB_LP <- output7IIAXIBc$bifpoints[2,];
Eq_K2_IIAXIB_LP2 <- output7IIAXIBf$bifpoints[1,];
### Codim ========================
### BP 1
cmd = output18IIAXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIAXIBc$bifpoints[3,c(1:6,8)], 1.0), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popBP", "1"), clean=FALSE);
cmd = output18IIAXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIAXIBc$bifpoints[3,c(1:6,8)], 1.0), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popBP", "1"), clean=FALSE);
output18IIAXIB<-rbind(output18IIAXIBa$curvepoints[dim(output18IIAXIBa$curvepoints)[1]:1,],
                      output18IIAXIBb$curvepoints[2:dim(output18IIAXIBb$curvepoints)[1],]);

### BPE
cmd = output20IIAXIBa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IIAXIBb$bifpoints[,c(1:3,7,8)], 1.0), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE);
cmd = output20IIAXIBb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IIAXIBb$bifpoints[,c(1:3,7,8)], 1.0), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE);
output20IIAXIB<-rbind(output20IIAXIBa$curvepoints[dim(output20IIAXIBa$curvepoints)[1]:1,],
                      output20IIAXIBb$curvepoints[2:dim(output20IIAXIBb$curvepoints)[1],]);
### LP
cmd = output21IIAXIBa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIAXIBc$bifpoints[2,c(1:8)], 1.0), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, NULL, clean=FALSE);
cmd = output21IIAXIBb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIAXIBc$bifpoints[2,c(1:8)], 1.0), 0.1, c(10, 0, 9.5E-04, 20, 0.1, 4.5), NULL, NULL, clean=FALSE);
output21IIAXIB<-rbind(output21IIAXIBa$curvepoints[dim(output21IIAXIBa$curvepoints)[1]:1,],
                      output21IIAXIBb$curvepoints[2:dim(output21IIAXIBb$curvepoints)[1],]);

cmd = output21IIAXIBc <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIAXIBf$bifpoints[1,c(1:5,7)], 0.5), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE);
cmd = output21IIAXIBd <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIAXIBf$bifpoints[1,c(1:5,7)], 0.5), 0.1, c(10, 0, 9.5E-04, 20, 0.1, 4.5), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE);
output21IIAXIB2<-rbind(output21IIAXIBc$curvepoints[dim(output21IIAXIBc$curvepoints)[1]:1,],
                       output21IIAXIBd$curvepoints[2:dim(output21IIAXIBd$curvepoints)[1],]);

################## IIB. K1=8e-5 ========================
### EQ ========================
# with Beta = 2.0 To detect BPE, LP and unstable BP C2
cmd = output7IIBXIB <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output2.2$curvepoints[53,c(3,6,8)]),
                                0.5, c(10, 0, 4E-4), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7IIBXIBb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIBXIB$bifpoints[1,c(1:3,5:8)]),
                                 0.1, c(10, 0, 9E-4), NULL, c("envZE", "2"), clean=FALSE)
cmd = output7IIBXIBc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIBXIBb$bifpoints[,c(1:8)]),
                                 -0.1, c(10, 0, 6E-4), NULL, NULL, clean=FALSE)

# with Beta = 1.0 Detection of BP C1 and C2 in presence of P
cmd = output7IIBXIBd <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1.2e-04,output3$bifpoints[3,c(2:6,7,8)]),
                                 -0.1, c(10, 0, 9E-1), NULL, NULL, clean=FALSE)
cmd = output7IIBXIBe <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1.2e-04,output3$bifpoints[3,c(2:6,7,8)]),
                                 0.1, c(10, 0, 9E-1), NULL, NULL, clean=FALSE)

### Starting from P-C1-R1+R2 with Beta=0.5 to detect LP#2 
cmd = output7IIBXIBf <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output3.2$curvepoints[141,c(3,4,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output7IIBXIBg <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIBXIBf$bifpoints[1,c(1:8)]), 0.1, c(10, 0, 2E-3), NULL, NULL, clean=FALSE)
cmd = output7IIBXIBh <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIBXIBg$bifpoints[1,c(1:5,7)]), -0.1, c(10, 0, 2E-3), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

# with Beta = 2 Continuing output7IIBXIBc > fixing K2=1E-3 > Beta gradient to detect LP#3
cmd = output7IIBXIBi <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIBXIBc$bifpoints[2,c(1:4,6,8)]), 0.1, c(10, 0, 1.1E-3), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output7IIBXIBj <- PSPMequi("Tritrophicmod5_fin", "EQ",c(2, output7IIBXIBi$curvepoints[120,c(2:4,6,8)]), 0.1, c(20, 0, 4.1), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output7IIBXIBk <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIBXIBj$bifpoints[5,c(1:4,6,8)]), 0.6, c(20, 0, 3.8), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K2_IIBXIB_BP1 <- output7IIBXIBe$bifpoints[1,]; 
Eq_K2_IIBXIB_BP2 <- output7IIBXIB$bifpoints[1,]; #inv of C2 in C1-R1+R2
Eq_K2_IIBXIB_BP2b <- output7IIBXIBc$bifpoints[2,]; #inv of C2 on unstable branch
Eq_K2_IIBXIB_BP2c <- output7IIBXIBd$bifpoints[1,]; #inv of C2 in P-C1-R1+R2
Eq_K2_IIBXIB_BPE <- output7IIBXIBb$bifpoints[1,];
Eq_K2_IIBXIB_LP <- output7IIBXIBc$bifpoints[1,];
Eq_K2_IIBXIB_LPb <- output7IIBXIBh$bifpoints[2,];

Eq_K2_IIBXIB_LPc <- output7IIBXIBk$bifpoints[2,];
### Codim ========================
### BP 1
cmd = output18IIBXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIBXIBe$bifpoints[1,c(1:4,6,8)], 1), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE);
cmd = output18IIBXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIBXIBe$bifpoints[1,c(1:4,6,8)], 1), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE);
output18IIBXIB<-rbind(output18IIBXIBa$curvepoints[dim(output18IIBXIBa$curvepoints)[1]:1,],
                      output18IIBXIBb$curvepoints[2:dim(output18IIBXIBb$curvepoints)[1],]);
### BP 2
cmd = output19IIBXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIBXIBd$bifpoints[1,c(1:4,5,7)], 1), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE);
cmd = output19IIBXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIBXIBd$bifpoints[1,c(1:4,5,7)], 1), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE);
output19IIBXIB<-rbind(output19IIBXIBa$curvepoints[dim(output19IIBXIBa$curvepoints)[1]:1,],
                      output19IIBXIBb$curvepoints[2:dim(output19IIBXIBb$curvepoints)[1],]);

cmd = output19IIBXIBc <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIBXIBc$bifpoints[2,c(1:5,8)], 2), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE);
cmd = output19IIBXIBd <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIBXIBc$bifpoints[2,c(1:5,8)], 2), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE);
output19IIBXIB2<-rbind(output19IIBXIBc$curvepoints[dim(output19IIBXIBc$curvepoints)[1]:1,],
                       output19IIBXIBd$curvepoints[2:dim(output19IIBXIBd$curvepoints)[1],]);

### BPE
cmd = output20IIBXIBa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IIBXIBb$bifpoints[1,c(1:3,7,8)], 2), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)
cmd = output20IIBXIBb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IIBXIBb$bifpoints[1,c(1:3,7,8)], 2), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)
output20IIBXIB<-rbind(output20IIBXIBa$curvepoints[dim(output20IIBXIBa$curvepoints)[1]:1,],
                      output20IIBXIBb$curvepoints[2:dim(output20IIBXIBb$curvepoints)[1],]);
### LP
cmd = output21IIBXIBa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIBXIBc$bifpoints[1,c(1:8)], 2), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, NULL, clean=FALSE);
cmd = output21IIBXIBb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIBXIBc$bifpoints[1,c(1:8)], 2), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, NULL, clean=FALSE);
output21IIBXIB<-rbind(output21IIBXIBa$curvepoints[dim(output21IIBXIBa$curvepoints)[1]:1,],
                      output21IIBXIBb$curvepoints[2:dim(output21IIBXIBb$curvepoints)[1],]);

cmd = output21IIBXIBc <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIBXIBh$bifpoints[2,c(1:5,7)], 0.5), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE);
cmd = output21IIBXIBd <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIBXIBh$bifpoints[2,c(1:5,7)], 0.5), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE);
output21IIBXIB2<-rbind(output21IIBXIBc$curvepoints[dim(output21IIBXIBc$curvepoints)[1]:1,],
                       output21IIBXIBd$curvepoints[2:dim(output21IIBXIBd$curvepoints)[1],]);
################## IIC. K1=1e-4 ========================
### EQ ========================
# with Beta = 2.0 To detect BPE, LP and unstable BP C2
cmd = output7IICXIB <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output2.2$curvepoints[63,c(3,6,8)]),
                                0.5, c(10, 0, 4E-4), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7IICXIBb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IICXIB$bifpoints[1,c(1:3,6:8)]),
                                 0.1, c(10, 0, 9E-4), NULL, c("envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7IICXIBc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IICXIBb$bifpoints[,c(1:8)]),
                                 -0.1, c(10, 0, 6E-4), NULL, NULL, clean=FALSE)

# with Beta = 1.0 Detection of BP C1 and C2 in presence of P
cmd = output7IICXIBd <- PSPMequi("Tritrophicmod5_fin", "EQ",c(2.3e-04,output3$bifpoints[2,c(2:6,7,8)]),
                                 -0.1, c(10, 0, 9E-1), NULL, NULL, clean=FALSE)
cmd = output7IICXIBe <- PSPMequi("Tritrophicmod5_fin", "EQ",c(2.3e-04,output3$bifpoints[2,c(2:6,7,8)]),
                                 0.1, c(10, 0, 9E-1), NULL, NULL, clean=FALSE)

### Starting from P-C1-R1+R2 with Beta=0.5 to detect LP#2 
cmd = output7IICXIBf <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output3.2$curvepoints[157,c(3,4,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output7IICXIBg <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IICXIBf$bifpoints[1,c(1:8)]), 0.1, c(10, 0, 2E-3), NULL, NULL, clean=FALSE)
cmd = output7IICXIBh <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IICXIBg$bifpoints[1,c(1:5,7)]), -0.1, c(10, 0, 2E-3), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

# with Beta = 2 Continuing output7IIBXIBc > fixing K2=1E-3 > Beta gradient to detect LP#3
cmd = output7IICXIBi <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IICXIBc$bifpoints[5,c(1:4,6,8)]), 0.1, c(10, 0, 1.1E-3), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output7IICXIBj <- PSPMequi("Tritrophicmod5_fin", "EQ",c(2, output7IICXIBi$curvepoints[121,c(2:4,6,8)]), 0.1, c(20, 0, 4.1), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K2_IICXIB_BP1 <- output7IICXIBe$bifpoints[1,]; 
Eq_K2_IICXIB_BP2 <- output7IICXIB$bifpoints[1,]; #inv of C2 in C1-R1+R2
Eq_K2_IICXIB_BP2b <- output7IICXIBc$bifpoints[5,]; #inv of C2 on unstable branch
Eq_K2_IICXIB_BP2c <- output7IICXIBd$bifpoints[1,]; #inv of C2 in P-C1-R1+R2
Eq_K2_IICXIB_BPE <- output7IICXIBb$bifpoints[1,];
Eq_K2_IICXIB_LP <- output7IICXIBc$bifpoints[2,];
Eq_K2_IICXIB_LPb <- output7IICXIBh$bifpoints[1,];

Eq_K2_IICXIB_LPc <- output7IICXIBj$bifpoints[2,];
### Codim ========================
### BP 1
cmd = output18IICXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IICXIBe$bifpoints[,c(1:4,6,7)], 1), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE);
cmd = output18IICXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IICXIBe$bifpoints[,c(1:4,6,7)], 1), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE);
output18IICXIB<-rbind(output18IICXIBa$curvepoints[dim(output18IICXIBa$curvepoints)[1]:1,],
                      output18IICXIBb$curvepoints[2:dim(output18IICXIBb$curvepoints)[1],]);
### BP 2
cmd = output19IICXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IICXIBd$bifpoints[,c(1:5,8)], 1), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE);
cmd = output19IICXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IICXIBd$bifpoints[,c(1:5,8)], 1), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE);
output19IICXIB<-rbind(output19IICXIBa$curvepoints[dim(output19IICXIBa$curvepoints)[1]:1,],
                      output19IICXIBb$curvepoints[2:dim(output19IICXIBb$curvepoints)[1],]);

cmd = output19IICXIBc <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IICXIBc$bifpoints[5,c(1:5,8)], 2), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE);
cmd = output19IICXIBd <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IICXIBc$bifpoints[5,c(1:5,8)], 2), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE);
output19IICXIB2<-rbind(output19IICXIBc$curvepoints[dim(output19IICXIBc$curvepoints)[1]:1,],
                       output19IICXIBd$curvepoints[2:dim(output19IICXIBd$curvepoints)[1],]);
### BPE
cmd = output20IICXIBa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IICXIBb$bifpoints[1,c(1:3,7,8)], 2), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)
cmd = output20IICXIBb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IICXIBb$bifpoints[1,c(1:3,7,8)], 2), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)
output20IICXIB<-rbind(output20IICXIBa$curvepoints[dim(output20IICXIBa$curvepoints)[1]:1,],
                      output20IICXIBb$curvepoints[2:dim(output20IICXIBb$curvepoints)[1],]);
### LP
cmd = output21IICXIBa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IICXIBc$bifpoints[2,c(1:8)], 2), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, NULL, clean=FALSE);
cmd = output21IICXIBb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IICXIBc$bifpoints[2,c(1:8)], 2), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, NULL, clean=FALSE);
output21IICXIB<-rbind(output21IICXIBa$curvepoints[dim(output21IICXIBa$curvepoints)[1]:1,],
                      output21IICXIBb$curvepoints[2:dim(output21IICXIBb$curvepoints)[1],]);

cmd = output21IICXIBc <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IICXIBh$bifpoints[1,c(1:5,7)], 0.5), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE);
cmd = output21IICXIBd <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IICXIBh$bifpoints[1,c(1:5,7)], 0.5), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE);
output21IICXIB2<-rbind(output21IICXIBc$curvepoints[dim(output21IICXIBc$curvepoints)[1]:1,],
                       output21IICXIBd$curvepoints[2:dim(output21IICXIBd$curvepoints)[1],]);
################## IIIA. K1=2e-4 ========================
### EQ ========================
cmd = output7IIIAXIB <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output2.2$curvepoints[78,c(3,6,8)]),
                                 0.5, c(10, 0, 4E-4), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7IIIAXIBb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIIAXIB$bifpoints[1,c(1:3,6:8)]),
                                  0.1, c(10, 0, 9E-4), NULL, c("envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7IIIAXIBc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIIAXIBb$bifpoints[,c(1:8)]),
                                  -0.1, c(10, 0, 6E-4), NULL, NULL, clean=FALSE)

# with Beta = 1.0 Detection of BP C1 and C2 in presence of P
cmd = output7IIIAXIBd <- PSPMequi("Tritrophicmod5_fin", "EQ",c(3.4e-04,output3$bifpoints[3,c(2:6,7,8)]),
                                  -0.1, c(10, 0, 9E-1), NULL, NULL, clean=FALSE)
cmd = output7IIIAXIBe <- PSPMequi("Tritrophicmod5_fin", "EQ",c(3.4e-04,output3$bifpoints[3,c(2:6,7,8)]),
                                  0.1, c(10, 0, 9E-1), NULL, NULL, clean=FALSE)

### Starting from P-C1-R1+R2 with Beta=0.5 to detect LP#2 
cmd = output7IIIAXIBf <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output3.2$curvepoints[195,c(3,4,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output7IIIAXIBg <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIIAXIBf$bifpoints[1,c(1:8)]), 0.1, c(10, 0, 2E-3), NULL, NULL, clean=FALSE)
cmd = output7IIIAXIBh <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIIAXIBg$bifpoints[1,c(1:5,7)]), 0.1, c(10, 0, 2E-3), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

# with Beta = 2.3 Continuing output7IIBXIBc > fixing K2=1E-3 > Beta gradient to detect LP#3
cmd = output7IIIAXIBi <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIIAXIBc$bifpoints[4,c(1:4,6,8)]), 0.1, c(10, 0, 1.1E-3), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output7IIIAXIBj <- PSPMequi("Tritrophicmod5_fin", "EQ",c(2.3, output7IIIAXIBi$curvepoints[121,c(2:4,6,8)]), 0.1, c(20, 0, 4.1), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K2_IIIAXIB_BP1 <- output7IIIAXIBe$bifpoints[1,]; 
Eq_K2_IIIAXIB_BP2 <- output7IIIAXIB$bifpoints[1,]; #inv of C2 in C1-R1+R2
Eq_K2_IIIAXIB_BP2b <- output7IIIAXIBc$bifpoints[4,]; #inv of C2 on unstable branch
Eq_K2_IIIAXIB_BP2c <- output7IIIAXIBd$bifpoints[1,]; #inv of C2 in P-C1-R1+R2
Eq_K2_IIIAXIB_BPE <- output7IIIAXIBb$bifpoints[1,];
Eq_K2_IIIAXIB_LP <- output7IIIAXIBc$bifpoints[1,];
Eq_K2_IIIAXIB_LPb <- output7IIIAXIBh$bifpoints[2,];

Eq_K2_IIIAXIB_LPc <- output7IIIAXIBj$bifpoints[2,];
### Codim ========================
### BP 1
cmd = output18IIIAXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIIAXIBe$bifpoints[,c(1:4,6,7)], 1), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE);
cmd = output18IIIAXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIIAXIBe$bifpoints[,c(1:4,6,7)], 1), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE);
output18IIIAXIB<-rbind(output18IIIAXIBa$curvepoints[dim(output18IIIAXIBa$curvepoints)[1]:1,],
                       output18IIIAXIBb$curvepoints[2:dim(output18IIIAXIBb$curvepoints)[1],]);
### BP 2
cmd = output19IIIAXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIIAXIBd$bifpoints[,c(1:5,8)], 1), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE);
cmd = output19IIIAXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIIAXIBd$bifpoints[,c(1:5,8)], 1), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE);
output19IIIAXIB<-rbind(output19IIIAXIBa$curvepoints[dim(output19IIIAXIBa$curvepoints)[1]:1,],
                       output19IIIAXIBb$curvepoints[2:dim(output19IIIAXIBb$curvepoints)[1],]);

cmd = output19IIIAXIBc <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIIAXIBc$bifpoints[4,c(1:5,8)], 2.3), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE);
cmd = output19IIIAXIBd <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIIAXIBc$bifpoints[4,c(1:5,8)], 2.3), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE);
output19IIIAXIB2<-rbind(output19IIIAXIBc$curvepoints[dim(output19IIIAXIBc$curvepoints)[1]:1,],
                        output19IIIAXIBd$curvepoints[2:dim(output19IIIAXIBd$curvepoints)[1],]);
### BPE
cmd = output20IIIAXIBa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IIIAXIBb$bifpoints[1,c(1:3,7,8)], 2.3), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)
cmd = output20IIIAXIBb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IIIAXIBb$bifpoints[1,c(1:3,7,8)], 2.3), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)

output20IIIAXIB<-rbind(output20IIIAXIBa$curvepoints[dim(output20IIIAXIBa$curvepoints)[1]:1,],
                       output20IIIAXIBb$curvepoints[2:dim(output20IIIAXIBb$curvepoints)[1],]);
### LP
cmd = output21IIIAXIBa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIIAXIBc$bifpoints[1,c(1:8)], 2.3), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, NULL, clean=FALSE);
cmd = output21IIIAXIBb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIIAXIBc$bifpoints[1,c(1:8)], 2.3), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, NULL, clean=FALSE);
output21IIIAXIB<-rbind(output21IIIAXIBa$curvepoints[dim(output21IIIAXIBa$curvepoints)[1]:1,],
                       output21IIIAXIBb$curvepoints[2:dim(output21IIIAXIBb$curvepoints)[1],]);

cmd = output21IIIAXIBc <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIIAXIBh$bifpoints[2,c(1:5,7)], 0.5), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE);
cmd = output21IIIAXIBd <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIIAXIBh$bifpoints[2,c(1:5,7)], 0.5), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE);
output21IIIAXIB2<-rbind(output21IIIAXIBc$curvepoints[dim(output21IIIAXIBc$curvepoints)[1]:1,],
                        output21IIIAXIBd$curvepoints[2:dim(output21IIIAXIBd$curvepoints)[1],]);

cmd = output21IIIAXIBe <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIIAXIBl$bifpoints[1,c(1:8)], 1E-3), -0.1, c(20, 0, 4.1, 10, 0.1, 9.5E-01), NULL, NULL, clean=FALSE);
cmd = output21IIIAXIBf <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIIAXIBl$bifpoints[1,c(1:5,7)], 0.5), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE);
output21IIIAXIB3<-rbind(output21IIIAXIBe$curvepoints[dim(output21IIIAXIBe$curvepoints)[1]:1,],
                        output21IIIAXIBf$curvepoints[2:dim(output21IIIAXIBf$curvepoints)[1],]);

################## IIIB. K1=3e-4 ========================
### EQ ========================
cmd = output7IIIBXIB <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output2.2$curvepoints[82,c(3,6,8)]),
                                 0.5, c(10, 0, 4E-4), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7IIIBXIBb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIIBXIB$bifpoints[1,c(1:3,6:8)]),
                                  0.1, c(10, 0, 9E-4), NULL, c("envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7IIIBXIBc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIIBXIBb$bifpoints[,c(1:8)]),
                                  -0.1, c(10, 0, 6E-4), NULL, NULL, clean=FALSE)

# with Beta = 1.0 Detection of BP C1 and C2 in presence of P
cmd = output7IIIBXIBd <- PSPMequi("Tritrophicmod5_fin", "EQ",c(4.5e-04,output3$bifpoints[2,c(2:6,7,8)]),
                                  -0.1, c(10, 0, 9E-1), NULL, NULL, clean=FALSE)
cmd = output7IIIBXIBe <- PSPMequi("Tritrophicmod5_fin", "EQ",c(4.5e-04,output3$bifpoints[2,c(2:6,7,8)]),
                                  0.1, c(10, 0, 9E-1), NULL, NULL, clean=FALSE)

# Starting from P-C1-R1+R2 with Beta=0.5 to detect LP#2 
cmd = output7IIIBXIBf <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output3.2$curvepoints[202,c(3,4,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output7IIIBXIBg <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIIBXIBf$bifpoints[1,c(1:8)]), 0.1, c(10, 0, 2E-3), NULL, NULL, clean=FALSE)
cmd = output7IIIBXIBh <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIIBXIBg$bifpoints[1,c(1:5,7)]), -0.1, c(10, 0, 2E-3), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

# with Beta = 2.5 Continuing output7IIBXIBc > fixing K2=1E-3 > Beta gradient to detect LP#3
cmd = output7IIIBXIBi <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIIBXIBc$bifpoints[2,c(1:4,6,8)]), 0.1, c(10, 0, 1.1E-3), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output7IIIBXIBj <- PSPMequi("Tritrophicmod5_fin", "EQ",c(2.5, output7IIIBXIBi$curvepoints[121,c(2:4,6,8)]), 0.1, c(20, 0, 4.1), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K2_IIIBXIB_BP1 <- output7IIIBXIBe$bifpoints[1,]; 
Eq_K2_IIIBXIB_BP2 <- output7IIIBXIB$bifpoints[1,]; #inv of C2 in C1-R1+R2
Eq_K2_IIIBXIB_BP2b <- output7IIIBXIBc$bifpoints[2,]; #inv of C2 on unstable branch
Eq_K2_IIIBXIB_BP2c <- output7IIIBXIBd$bifpoints[1,]; #inv of C2 in P-C1-R1+R2
Eq_K2_IIIBXIB_BPE <- output7IIIBXIBb$bifpoints[1,];
Eq_K2_IIIBXIB_LP <- output7IIIBXIBc$bifpoints[1,];
Eq_K2_IIIBXIB_LPb <- output7IIIBXIBh$bifpoints[2,];

Eq_K2_IIIBXIB_LPc <- output7IIIBXIBj$bifpoints[2,];
### Codim ========================
### BP 1
cmd = output18IIIBXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIIBXIBe$bifpoints[,c(1:4,6,8)], 1), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output18IIIBXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIIBXIBe$bifpoints[,c(1:4,6,8)], 1), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
output18IIIBXIB<-rbind(output18IIIBXIBa$curvepoints[dim(output18IIIBXIBa$curvepoints)[1]:1,],
                       output18IIIBXIBb$curvepoints[2:dim(output18IIIBXIBb$curvepoints)[1],]);
### BP 2
cmd = output19IIIBXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIIBXIBd$bifpoints[,c(1:4,6,8)], 1), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
cmd = output19IIIBXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIIBXIBd$bifpoints[,c(1:4,6,8)], 1), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
output19IIIBXIB<-rbind(output19IIIBXIBa$curvepoints[dim(output19IIIBXIBa$curvepoints)[1]:1,],
                       output19IIIBXIBb$curvepoints[2:dim(output19IIIBXIBb$curvepoints)[1],]);

cmd = output19IIIBXIBc <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIIBXIBc$bifpoints[2,c(1:5,8)], 2.5), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE);
cmd = output19IIIBXIBd <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIIBXIBc$bifpoints[2,c(1:5,8)], 2.5), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE);
output19IIIBXIB2<-rbind(output19IIIBXIBc$curvepoints[dim(output19IIIBXIBc$curvepoints)[1]:1,],
                        output19IIIBXIBd$curvepoints[2:dim(output19IIIBXIBd$curvepoints)[1],]);
### BPE
cmd = output20IIIBXIBa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IIIBXIBb$bifpoints[1,c(1:3,7,8)], 2.5), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)
cmd = output20IIIBXIBb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IIIBXIBb$bifpoints[1,c(1:3,7,8)], 2.5), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)
output20IIIBXIB<-rbind(output20IIIBXIBa$curvepoints[dim(output20IIIBXIBa$curvepoints)[1]:1,],
                       output20IIIBXIBb$curvepoints[2:dim(output20IIIBXIBb$curvepoints)[1],]);
### LP
cmd = output21IIIBXIBa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIIBXIBc$bifpoints[1,c(1:8)], 2.5), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, NULL, clean=FALSE);
cmd = output21IIIBXIBb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIIBXIBc$bifpoints[1,c(1:8)], 2.5), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, NULL, clean=FALSE);
output21IIIBXIB<-rbind(output21IIIBXIBa$curvepoints[dim(output21IIIBXIBa$curvepoints)[1]:1,],
                       output21IIIBXIBb$curvepoints[2:dim(output21IIIBXIBb$curvepoints)[1],]);

cmd = output21IIIBXIBc <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIIBXIBh$bifpoints[2,c(1:5,7)], 0.5), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE);
cmd = output21IIIBXIBd <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIIBXIBh$bifpoints[2,c(1:5,7)], 0.5), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE);
output21IIIBXIB2<-rbind(output21IIIBXIBc$curvepoints[dim(output21IIIBXIBc$curvepoints)[1]:1,],
                        output21IIIBXIBd$curvepoints[2:dim(output21IIIBXIBd$curvepoints)[1],]);
################## IIIC. K1=4e-4 ========================
### EQ ========================
cmd = output7IIICXIB <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output2.2$curvepoints[87,c(3,6,8)]),
                                 0.5, c(10, 0, 4E-4), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7IIICXIBb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIICXIB$bifpoints[1,c(1:3,5:8)]),
                                  0.1, c(10, 0, 9E-4), NULL, c("envZE", "2"), clean=FALSE)
cmd = output7IIICXIBc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIICXIBb$bifpoints[,c(1:8)]),
                                  -0.1, c(10, 0, 6E-4), NULL, NULL, clean=FALSE)

# with Beta = 1.0 Detection of BP C1 and C2 in presence of P
cmd = output7IIICXIBd <- PSPMequi("Tritrophicmod5_fin", "EQ",c(4.5e-04,output3$bifpoints[2,c(2:6,7,8)]),
                                  -0.1, c(10, 0, 9E-1), NULL, NULL, clean=FALSE)
cmd = output7IIICXIBe <- PSPMequi("Tritrophicmod5_fin", "EQ",c(4.5e-04,output3$bifpoints[2,c(2:6,7,8)]),
                                  0.1, c(10, 0, 9E-1), NULL, NULL, clean=FALSE)

### Starting from P-C1-R1+R2 with Beta=0.5 to detect LP#2 
cmd = output7IIICXIBf <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output3.2$curvepoints[289,c(3,4,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output7IIICXIBg <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIICXIBf$bifpoints[1,c(1:8)]), 0.1, c(10, 0, 2E-3), NULL, NULL, clean=FALSE)
cmd = output7IIICXIBh <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIICXIBg$bifpoints[1,c(1:5,7)]), -0.1, c(10, 0, 2E-3), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

# with Beta = 2.5 Continuing output7IIBXIBc > fixing K2=1E-3 > Beta gradient to detect LP#3
cmd = output7IIICXIBi <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIICXIBc$bifpoints[4,c(1:4,6,8)]), 0.1, c(10, 0, 1.1E-3), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output7IIICXIBj <- PSPMequi("Tritrophicmod5_fin", "EQ",c(2.5, output7IIICXIBi$curvepoints[121,c(2:4,6,8)]), 0.1, c(20, 0, 4.1), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K2_IIICXIB_BP1 <- output7IIICXIBe$bifpoints[1,]; 
Eq_K2_IIICXIB_BP2 <- output7IIICXIB$bifpoints[1,]; #inv of C2 in C1-R1+R2
Eq_K2_IIICXIB_BP2b <- output7IIICXIBc$bifpoints[4,]; #inv of C2 on unstable branch
Eq_K2_IIICXIB_BP2c <- output7IIICXIBd$bifpoints[1,]; #inv of C2 in P-C1-R1+R2
Eq_K2_IIICXIB_BPE <- output7IIICXIBb$bifpoints[1,];
Eq_K2_IIICXIB_LP <- output7IIICXIBc$bifpoints[1,];
Eq_K2_IIICXIB_LPb <- output7IIICXIBh$bifpoints[2,];

Eq_K2_IIICXIB_LPc <- output7IIICXIBj$bifpoints[2,];
### Codim ========================
### BP 1 
cmd = output18IIICXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIICXIBe$bifpoints[,c(1:5,7)], 1), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output18IIICXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIICXIBe$bifpoints[,c(1:5,7)], 1), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
output18IIICXIB<-rbind(output18IIICXIBa$curvepoints[dim(output18IIICXIBa$curvepoints)[1]:1,],
                       output18IIICXIBb$curvepoints[2:dim(output18IIICXIBb$curvepoints)[1],]);
### BP 2
cmd = output19IIICXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIICXIBd$bifpoints[,c(1:5,8)], 1), -0.1, c(10, 0, 9.5E-1, 20, 0.1, 4.5), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
cmd = output19IIICXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIICXIBd$bifpoints[,c(1:5,8)], 1), 0.1, c(10, 0, 9.5E-1, 20, 0.1, 4.5), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
output19IIICXIB<-rbind(output19IIICXIBa$curvepoints[dim(output19IIICXIBa$curvepoints)[1]:1,],
                       output19IIICXIBb$curvepoints[2:dim(output19IIICXIBb$curvepoints)[1],]);

cmd = output19IIICXIBc <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIICXIBc$bifpoints[4,c(1:5,8)], 2.5), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE);
cmd = output19IIICXIBd <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IIICXIBc$bifpoints[4,c(1:5,8)], 2.5), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE);
output19IIICXIB2<-rbind(output19IIICXIBc$curvepoints[dim(output19IIICXIBc$curvepoints)[1]:1,],
                        output19IIICXIBd$curvepoints[2:dim(output19IIICXIBd$curvepoints)[1],]);

### BPE
cmd = output20IIICXIBa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IIICXIBb$bifpoints[1,c(1:3,7,8)], 2.5), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)
cmd = output20IIICXIBb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IIICXIBb$bifpoints[1,c(1:3,7,8)], 2.5), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)
output20IIICXIB<-rbind(output20IIICXIBa$curvepoints[dim(output20IIICXIBa$curvepoints)[1]:1,],
                       output20IIICXIBb$curvepoints[2:dim(output20IIICXIBb$curvepoints)[1],]);
### LP
cmd = output21IIICXIBa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIICXIBc$bifpoints[1,c(1:8)], 2.5), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, NULL, clean=FALSE);
cmd = output21IIICXIBb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIICXIBc$bifpoints[1,c(1:8)], 2.5), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, NULL, clean=FALSE);
output21IIICXIB<-rbind(output21IIICXIBa$curvepoints[dim(output21IIICXIBa$curvepoints)[1]:1,],
                       output21IIICXIBb$curvepoints[2:dim(output21IIICXIBb$curvepoints)[1],]);

cmd = output21IIICXIBc <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIICXIBh$bifpoints[2,c(1:5,7)], 0.5), -0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE);
cmd = output21IIICXIBd <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIICXIBh$bifpoints[2,c(1:5,7)], 0.5), 0.1, c(10, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE);
output21IIICXIB2<-rbind(output21IIICXIBc$curvepoints[dim(output21IIICXIBc$curvepoints)[1]:1,],
                        output21IIICXIBd$curvepoints[2:dim(output21IIICXIBd$curvepoints)[1],]);

cmd = output21IIICXIBe <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIICXIBo$bifpoints[1,c(1:8)], 4E-4), -0.1, c(20, 0, 4.5, 10, 0, 9.5E-01), NULL, NULL, clean=FALSE);
cmd = output21IIICXIBf <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIICXIBn$bifpoints[1,c(1:8,7)], 0.5), 0.1, c(20, 0, 9.5E-01, 20, 0.1, 4.5), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE);
output21IIICXIB3<-rbind(output21IIICXIBe$curvepoints[dim(output21IIICXIBe$curvepoints)[1]:1,],
                        output21IIICXIBf$curvepoints[2:dim(output21IIICXIBf$curvepoints)[1],]);
output7IIICXIBj$bifpoints