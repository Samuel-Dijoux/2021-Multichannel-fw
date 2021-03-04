###################################################################################### 
#                -> Bifurcation analyses over K1 and K2 gradients for varying Beta            
#
# Codimension bifurcation tracking selected BP
# as a function of two parameters : Beta and K1 + K2
# Representation of the coexistence phase of both consumers:  P-C1C2-R1R2 ~ K1+K2
######################################################################################

## source("Tritrophicmod_PCiRi.r")
## source("Tritrophicmod_PC1R1.r")
## source("Tritrophicmod_PC2R2.r")

################## OA. BETA = 0.3 ========================
### EQ ========================
### Along K2 gradient
# Starting from C1-R1+R2 with K1=1E-5
# inv C2 in C1-R1+R2 > inv P > detect LP#1
cmd = output7OAXIAa <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output2.2$curvepoints[4,c(3,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7OAXIAb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7OAXIAa$bifpoints[1,c(1:3,5:8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("envZE", "2"), clean=FALSE)
cmd = output7OAXIAc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7OAXIAb$bifpoints[1,c(1:8)]), -0.1, c(10, 0, 9.5E-4), NULL, NULL, clean=FALSE)

# Starting from P-C2-R1+R2 with K1=2E-5
# inv C1 in P-C1-R1+R2 > ext C1 + detect LP#2 + Bistability 5/6
cmd = output7OAXIBa <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output3.2$curvepoints[28,c(3,4,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output7OAXIBb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7OAXIBa$bifpoints[1,c(1:8)]), 0.1, c(10, 0, 9.5E-4), NULL, NULL, clean=FALSE)
cmd = output7OAXIBc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7OAXIBb$bifpoints[1,c(1:5,7)]), -0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

Eq_K2_OA_BP1 <- output7OAXIBb$bifpoints[1,];
Eq_K2_OA_BP2 <- output7OAXIAa$bifpoints[1,]; #inv of C2 in C1-R1+R2
Eq_K2_OA_BPE <- output7OAXIAb$bifpoints[1,];
Eq_K2_OA_LP <- output7OAXIAc$bifpoints[1,]; # LP1
Eq_K2_OA_LPb <- output7OAXIBc$bifpoints[2,]; # LP1

### Codim ========================
### BP 1
cmd = output10OAXIAa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7OAXIBb$bifpoints[1,c(1:5,7)], 1E-5), -0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output10OAXIAb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7OAXIBb$bifpoints[1,c(1:5,7)], 1E-5), 0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
output10OAXIA<-rbind(output10OAXIAa$curvepoints[dim(output10OAXIAa$curvepoints)[1]:1,],
                     output10OAXIAb$curvepoints[2:dim(output10OAXIAb$curvepoints)[1],]);

### BP 2
cmd = output11OAXIAa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7OAXIAa$bifpoints[1,c(1:4,6,8)], 1E-5), -0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
cmd = output11OAXIAb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7OAXIAa$bifpoints[1,c(1:4,6,8)], 1E-5), 0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
output11OAXIA <-rbind(output11OAXIAa$curvepoints[dim(output11OAXIAa$curvepoints)[1]:1,],
                      output11OAXIAb$curvepoints[2:dim(output11OAXIAb$curvepoints)[1],]);

### BPE
cmd = output12OAXIAa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7OAXIAb$bifpoints[1,c(1:3,5:8)], 1E-5), -0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, c("envBP", "2"), clean=FALSE)
cmd = output12OAXIAb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7OAXIAb$bifpoints[1,c(1:3,5:8)], 1E-5), 0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, c("envBP", "2"), clean=FALSE)
output12OAXIA<-rbind(output12OAXIAa$curvepoints[dim(output12OAXIAa$curvepoints)[1]:1,],
                     output12OAXIAb$curvepoints[2:dim(output12OAXIAb$curvepoints)[1],]);

### LP
cmd = output13OAXIAa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7OAXIAc$bifpoints[1,c(1:8)], 1E-5), -1.1, c(10, 0, 9.5-01, 11, -0.6, 9.5E-01), NULL, NULL, clean=FALSE)
cmd = output13OAXIAb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7OAXIAc$bifpoints[1,c(1:8)], 1E-5), 0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, NULL, clean=FALSE)
cmd = output13OAXIAc <- PSPMequi("Tritrophicmod5_fin", "LP", c(output13OAXIAa$curvepoints[dim(output13OAXIAa$curvepoints)[1]-3,c(1:8)], 1E-5), -0.8, c(10, 0, 9.5-01, 11, 0.9, 9.5E-01), NULL, NULL, clean=FALSE)
cmd = output13OAXIAd <- PSPMequi("Tritrophicmod5_fin", "LP", c(output13OAXIAc$curvepoints[dim(output13OAXIAc$curvepoints)[1]-3,c(1:8)], 1E-5), -0.8, c(10, 0, 9.5-01, 11, 0.9, 9.5E-01), NULL, NULL, clean=FALSE)

output13OAXIA<-rbind(output13OAXIAd$curvepoints[dim(output13OAXIAd$curvepoints)[1]:1,],
                     output13OAXIAc$curvepoints[dim(output13OAXIAc$curvepoints)[1]:1,],
                     output13OAXIAa$curvepoints[dim(output13OAXIAa$curvepoints)[1]:1,],
                     output13OAXIAb$curvepoints[2:dim(output13OAXIAb$curvepoints)[1],]);

cmd = output13OAXIBa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7OAXIBc$bifpoints[2,c(1:5,7)], 2E-5), -0.1, c(10, 0, 9.5-01, 11, -0.6, 9.5E-01), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)
cmd = output13OAXIBb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7OAXIBc$bifpoints[2,c(1:5,7)], 2E-5), 0.1, c(10, 0, 9.5-01, 11, -0.6, 9.5E-01), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)
output13OAXIB<-rbind(output13OAXIBa$curvepoints[dim(output13OAXIBa$curvepoints)[1]:1,],
                     output13OAXIBb$curvepoints[2:dim(output13OAXIBb$curvepoints)[1],]);

################## OB. BETA = 0.5 ========================
### EQ ========================
### Along K2 gradient
# Starting from C1-R1+R2 with K1=1E-5
# inv C2 in C1-R1+R2 > inv P > detect LP#1 + Tristability 2/5/6 + ext C1 > detect LP#2 
cmd = output7OBXIBa <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output2.2$curvepoints[21,c(3,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7OBXIBb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7OBXIBa$bifpoints[1,c(1:3,5:8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("envZE", "2"), clean=FALSE)
cmd = output7OBXIBc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7OBXIBb$bifpoints[1,c(1:8)]), -0.1, c(10, 0, 9.5E-4), NULL, NULL, clean=FALSE)
cmd = output7OBXIBd <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7OBXIBc$bifpoints[3,c(1:5,7)]), -0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

### Along K1 gradient
# Starting from C2-R2+R1 with K2=1E-5
cmd = output8OBXIAa <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[7,c(2:3,5,7)]), 0.1, c(11, 0, 9.5E-4), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = output8OBXIAb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8OBXIAa$bifpoints[1,c(1:3,5:8)]), 0.1, c(11, 0, 9.5E-4), NULL, c("envZE", "2"), clean=FALSE)
cmd = output8OBXIAc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8OBXIAb$bifpoints[1,c(1:8)]), -0.1, c(11, 0, 9.5E-4), NULL, NULL, clean=FALSE)

Eq_K2_OB_BP1 <- output7OBXIBc$bifpoints[3,];
Eq_K2_OB_BP2 <- output7OBXIBa$bifpoints[1,]; #inv of C2 in C1-R1+R2
Eq_K1_OB_BP2 <- output8OBXIAc$bifpoints[2,];
Eq_K2_OB_BPE <- output7OBXIBb$bifpoints[1,];
Eq_K1_OB_BPE <- output8OBXIAb$bifpoints[1,];
Eq_K2_OB_LP <- output7OBXIBc$bifpoints[2,]; # LP1
Eq_K1_OB_LP <- output8OBXIAc$bifpoints[1,]; # LP1
Eq_K2_OB_LPb <- output7OBXIBd$bifpoints[1,]; # LP2

### Codim ========================
### BP 1
cmd = output10OBXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7OBXIBc$bifpoints[3,c(1:4,6,7)], 2E-5), -0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output10OBXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7OBXIBc$bifpoints[3,c(1:4,6,7)], 2E-5), 0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
output10OBXIB<-rbind(output10OBXIBa$curvepoints[dim(output10OBXIBa$curvepoints)[1]:1,],
                     output10OBXIBb$curvepoints[2:dim(output10OBXIBb$curvepoints)[1],]);
### BP 2
cmd = output11OBXIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8OBXIAc$bifpoints[2,c(1:4,6,8)], 1E-5), -0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
cmd = output11OBXIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8OBXIAc$bifpoints[2,c(1:4,6,8)], 1E-5), 0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
output11OBXIB <-rbind(output11OBXIBa$curvepoints[dim(output11OBXIBa$curvepoints)[1]:1,],
                      output11OBXIBb$curvepoints[2:dim(output11OBXIBb$curvepoints)[1],]);
### BPE
cmd = output12OBXIBa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7OBXIBb$bifpoints[1,c(1:3,5:8)], 2E-5), -0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, c("envBP", "2"), clean=FALSE)
cmd = output12OBXIBb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7OBXIBb$bifpoints[1,c(1:3,5:8)], 2E-5), 0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, c("envBP", "2"), clean=FALSE)
output12OBXIB <-rbind(output12OBXIBa$curvepoints[dim(output12OBXIBa$curvepoints)[1]:1,],
                      output12OBXIBb$curvepoints[2:dim(output12OBXIBb$curvepoints)[1],]);
### LP
cmd = output13OBXIBa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7OBXIBc$bifpoints[2,c(1:8)], 2E-5), -0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, NULL, clean=FALSE)
cmd = output13OBXIBb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7OBXIBc$bifpoints[2,c(1:8)], 2E-5), 0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, NULL, clean=FALSE)
output13OBXIB<-rbind(output13OBXIBa$curvepoints[dim(output13OBXIBa$curvepoints)[1]:1,],
                     output13OBXIBb$curvepoints[2:dim(output13OBXIBb$curvepoints)[1],]);

cmd = output13OBXIBc <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7OBXIBd$bifpoints[1,c(1:5,7)], 2E-5), -0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)
cmd = output13OBXIBd <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7OBXIBd$bifpoints[1,c(1:5,7)], 2E-5), 0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)
output13OBXIB2<-rbind(output13OBXIBc$curvepoints[dim(output13OBXIBc$curvepoints)[1]:1,],
                      output13OBXIBd$curvepoints[2:dim(output13OBXIBd$curvepoints)[1],]);
################## OC. BETA = 0.8 ========================
### EQ ========================
### Along K2 gradient
# Starting from C1-R1+R2 with K1=3E-5
cmd = output7OCXICa <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output2.2$curvepoints[27,c(3,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7OCXICb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7OCXICa$bifpoints[1,c(1:3,5:8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("envZE", "2"), clean=FALSE)
cmd = output7OCXICc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7OCXICb$bifpoints[1,c(1:8)]), -0.1, c(10, 0, 9.5E-4), NULL, NULL, clean=FALSE)
cmd = output7OCXICd <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7OCXICc$bifpoints[2,c(1:5,7)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

# Starting from P-C1-R1+R2 with K1=3E-4
cmd = output7OCXIIBa <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output3.2$curvepoints[222,c(3:4,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output7OCXIIBb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7OCXIIBa$bifpoints[1,c(1:8)]), 0.1, c(10, 0, 9.5E-4), NULL, NULL, clean=FALSE)
cmd = output7OCXIIBc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7OCXIIBb$bifpoints[1,c(1:5,7)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

### Along K1 gradient
# Starting from C2-R2+R1 with K2=3E-5
cmd = output8OCXICa <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[27,c(2:3,5,7)]), 0.1, c(11, 0, 9.5E-4), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = output8OCXICb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8OCXICa$bifpoints[1,c(1:3,5:8)]), 0.1, c(11, 0, 9.5E-4), NULL, c("envZE", "2"), clean=FALSE)
cmd = output8OCXICc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8OCXICb$bifpoints[1,c(1:8)]), -0.1, c(11, 0, 9.5E-4), NULL, NULL, clean=FALSE)
cmd = output8OCXICd <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8OCXICc$bifpoints[2,c(1:5,7)]), 0.1, c(11, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K2_OC_BP1 <- output7OCXICc$bifpoints[2,];
Eq_K1_OC_BP2 <- output8OCXICc$bifpoints[2,];
Eq_K1_OC_BPE <- output8OCXICb$bifpoints[1,];
Eq_K1_OC_LP <- output8OCXICc$bifpoints[1,];

### Codim ========================
### BP 1
cmd = output10OCXICa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7OCXICc$bifpoints[2,c(1:4,6,7)], 3E-5), -0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output10OCXICb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7OCXICc$bifpoints[2,c(1:4,6,7)], 3E-5), 0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
output10OCXIC<-rbind(output10OCXICa$curvepoints[dim(output10OCXICa$curvepoints)[1]:1,],
                     output10OCXICb$curvepoints[2:dim(output10OCXICb$curvepoints)[1],]);
### BP 2
cmd = output11OCXICa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8OCXICc$bifpoints[2,c(1:4,6,8)], 3E-5), -0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
cmd = output11OCXICb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8OCXICc$bifpoints[2,c(1:4,6,8)], 3E-5), 0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
output11OCXIC<-rbind(output11OCXICa$curvepoints[dim(output11OCXICa$curvepoints)[1]:1,],
                     output11OCXICb$curvepoints[2:dim(output11OCXICb$curvepoints)[1],]);
### BPE
cmd = output12OCXICa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output8OCXICb$bifpoints[1,c(1:3,5:8)], 3E-5), -0.1, c(11, 0, 9.5-01, 10, -0.1, 9.5E-01), NULL, c("envBP", "2"), clean=FALSE)
cmd = output12OCXICb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output8OCXICb$bifpoints[1,c(1:3,5:8)], 3E-5), 0.1, c(11, 0, 9.5-01, 10, -0.1, 9.5E-01), NULL, c("envBP", "2"), clean=FALSE)
output12OCXIC<-rbind(output12OCXICa$curvepoints[dim(output12OCXICa$curvepoints)[1]:1,],
                     output12OCXICb$curvepoints[2:dim(output12OCXICb$curvepoints)[1],]);
### LP
cmd = output13OCXICa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8OCXICc$bifpoints[1,c(1:8)], 3E-5), -0.1, c(11, 0, 9.5-01, 10, -0.1, 9.5E-01), NULL, NULL, clean=FALSE)
cmd = output13OCXICb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8OCXICc$bifpoints[1,c(1:8)], 3E-5), 0.1, c(11, 0, 9.5-01, 10, -0.1, 9.5E-01), NULL, NULL, clean=FALSE)
output13OCXIC<-rbind(output13OCXICa$curvepoints[dim(output13OCXICa$curvepoints)[1]:1,],
                     output13OCXICb$curvepoints[2:dim(output13OCXICb$curvepoints)[1],]);

################## IA. BETA = 1.0 ========================
### EQ ========================
### Along K2 gradient
# Starting from C1-R1+R2 with K1=3E-5
cmd = output7IAXICa <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output2.2$curvepoints[27,c(3,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7IAXICb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IAXICa$bifpoints[1,c(1:3,5:8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("envZE", "2"), clean=FALSE)
cmd = output7IAXICc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IAXICb$bifpoints[1,c(1:8)]), -0.1, c(10, 0, 9.5E-4), NULL, NULL, clean=FALSE)
cmd = output7IAXICd <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IAXICc$bifpoints[2,c(1:5,7)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

### Along K1 gradient
# Starting from C2-R2+R1 with K2=3E-5
cmd = output8IAXICa <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[27,c(2:3,5,7)]), 0.1, c(11, 0, 9.5E-4), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = output8IAXICb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IAXICa$bifpoints[1,c(1:3,5:8)]), 0.1, c(11, 0, 9.5E-4), NULL, c("envZE", "2"), clean=FALSE)
cmd = output8IAXICc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IAXICb$bifpoints[1,c(1:8)]), -0.1, c(11, 0, 9.5E-4), NULL, NULL, clean=FALSE)

Eq_K2_IA_BP1 <- output7IAXICc$bifpoints[2,];
Eq_K1_IA_BP2 <- output8IAXICc$bifpoints[2,];
Eq_K2_IA_BPE <- output7IAXICb$bifpoints[1,];
Eq_K2_IA_LP <- output7IAXICc$bifpoints[1,];

### Codim ========================
### BP 1
cmd = output10IAXICa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IAXICc$bifpoints[2,c(1:4,6,7)], 3E-5), -0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output10IAXICb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IAXICc$bifpoints[2,c(1:4,6,7)], 3E-5), 0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
output10IAXIC<-rbind(output10IAXICa$curvepoints[dim(output10IAXICa$curvepoints)[1]:1,],
                     output10IAXICb$curvepoints[2:dim(output10IAXICb$curvepoints)[1],]);

### BP 2
cmd = output11IAXICa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IAXICc$bifpoints[2,c(1:4,6,8)], 3E-5), -0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
cmd = output11IAXICb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IAXICc$bifpoints[2,c(1:4,6,8)], 3E-5), 0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
output11IAXIC<-rbind(output11IAXICa$curvepoints[dim(output11IAXICa$curvepoints)[1]:1,],
                     output11IAXICb$curvepoints[2:dim(output11IAXICb$curvepoints)[1],]);

### BPE curve
cmd = output12IAXICa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IAXICb$bifpoints[1,c(1:3,5:8)], 3E-5), -0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, c("envBP", "2"), clean=FALSE)
cmd = output12IAXICb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IAXICb$bifpoints[1,c(1:3,5:8)], 3E-5), 0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, c("envBP", "2"), clean=FALSE)
output12IAXIC<-rbind(output12IAXICa$curvepoints[dim(output12IAXICa$curvepoints)[1]:1,],
                     output12IAXICb$curvepoints[2:dim(output12IAXICb$curvepoints)[1],]);

### LP curve
cmd = output13IAXICa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IAXICc$bifpoints[1,c(1:8)], 3E-5), -0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, NULL, clean=FALSE)
cmd = output13IAXICb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IAXICc$bifpoints[1,c(1:8)], 3E-5), 0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, NULL, clean=FALSE)
output13IAXIC<-rbind(output13IAXICa$curvepoints[dim(output13IAXICa$curvepoints)[1]:1,],
                     output13IAXICb$curvepoints[2:dim(output13IAXICb$curvepoints)[1],]);

################## IAB. BETA = 1.2 ========================
### EQ ========================
### Along K2 gradient
# Starting from C1-R1+R2 with K1=3E-5
cmd = output7IABXICa <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output2.2$curvepoints[27,c(3,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7IABXICb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IABXICa$bifpoints[1,c(1:3,5:8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("envZE", "2"), clean=FALSE)
cmd = output7IABXICc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IABXICb$bifpoints[1,c(1:8)]), -0.1, c(10, 0, 9.5E-4), NULL, NULL, clean=FALSE)
cmd = output7IABXICd <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IABXICc$bifpoints[2,c(1:5,7)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

### Along K1 gradient
# Starting from C2-R2+R1 with K2=3E-5
cmd = output8IABXICa <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[27,c(2:3,5,7)]), 0.1, c(11, 0, 9.5E-4), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = output8IABXICb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IABXICa$bifpoints[1,c(1:3,5:8)]), 0.1, c(11, 0, 9.5E-4), NULL, c("envZE", "2"), clean=FALSE)
cmd = output8IABXICc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IABXICb$bifpoints[1,c(1:8)]), -0.1, c(11, 0, 9.5E-4), NULL, NULL, clean=FALSE)

Eq_K2_IAB_BP1 <- output7IABXICc$bifpoints[2,];
Eq_K1_IAB_BP2 <- output8IABXICc$bifpoints[2,];
Eq_K2_IAB_BPE <- output7IABXICb$bifpoints[1,];
Eq_K2_IAB_LP <- output7IABXICc$bifpoints[1,];

### Codim ========================
### BP 1
cmd = output10IABXICa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IABXICc$bifpoints[2,c(1:4,6,7)], 3E-5), -0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output10IABXICb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output7IABXICc$bifpoints[2,c(1:4,6,7)], 3E-5), 0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
output10IABXIC<-rbind(output10IABXICa$curvepoints[dim(output10IABXICa$curvepoints)[1]:1,],
                      output10IABXICb$curvepoints[2:dim(output10IABXICb$curvepoints)[1],]);

### BP 2
cmd = output11IABXICa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IABXICc$bifpoints[2,c(1:4,6,8)], 3E-5), -0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
cmd = output11IABXICb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IABXICc$bifpoints[2,c(1:4,6,8)], 3E-5), 0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
output11IABXIC<-rbind(output11IABXICa$curvepoints[dim(output11IABXICa$curvepoints)[1]:1,],
                      output11IABXICb$curvepoints[2:dim(output11IABXICb$curvepoints)[1],]);

### BPE curve
cmd = output12IABXICa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IABXICb$bifpoints[1,c(1:3,5:8)], 3E-5), -0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, c("envBP", "2"), clean=FALSE)
cmd = output12IABXICb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IABXICb$bifpoints[1,c(1:3,5:8)], 3E-5), 0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, c("envBP", "2"), clean=FALSE)
output12IABXIC<-rbind(output12IABXICa$curvepoints[dim(output12IABXICa$curvepoints)[1]:1,],
                      output12IABXICb$curvepoints[2:dim(output12IABXICb$curvepoints)[1],]);

### LP curve
cmd = output13IABXICa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IABXICc$bifpoints[1,c(1:8)], 3E-5), -0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, NULL, clean=FALSE)
cmd = output13IABXICb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IABXICc$bifpoints[1,c(1:8)], 3E-5), 0.1, c(10, 0, 9.5-01, 11, -0.1, 9.5E-01), NULL, NULL, clean=FALSE)
output13IABXIC<-rbind(output13IABXICa$curvepoints[dim(output13IABXICa$curvepoints)[1]:1,],
                      output13IABXICb$curvepoints[2:dim(output13IABXICb$curvepoints)[1],]);

################## IB. BETA = 1.5 ========================
### EQ ========================
### Along K2 gradient
# Starting from C1-R1+R2 with K1=3E-5
cmd = output7IBXICa <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output2.2$curvepoints[27,c(3,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7IBXICb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IBXICa$bifpoints[1,c(1:3,5:8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("envZE", "2"), clean=FALSE)
cmd = output7IBXICc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IBXICb$bifpoints[1,c(1:8)]), -0.1, c(10, 0, 9.5E-4), NULL, NULL, clean=FALSE)

### Along K1 gradient
# Starting from C2-R2+R1 with K2=3E-5
cmd = output8IBXICa <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[27,c(2:3,5,7)]), 0.1, c(11, 0, 9.5E-4), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = output8IBXICb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IBXICa$bifpoints[1,c(1:3,5:8)]), 0.1, c(11, 0, 9.5E-4), NULL, c("envZE", "2"), clean=FALSE)
cmd = output8IBXICc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IBXICb$bifpoints[1,c(1:8)]), -0.1, c(11, 0, 9.5E-4), NULL, NULL, clean=FALSE)
cmd = output8IBXICd <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IBXICc$bifpoints[2,c(1:5,7)]), 0.1, c(11, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

### Along K1 gradient
# Starting from P-C2-R2+R1 with K2=3E-4
cmd = output8IBXIIIBa <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output3$curvepoints[426,c(2:5,7)]), 0.1, c(11, 0, 9.5E-4), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)
cmd = output8IBXIIIBb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IBXIIIBa$bifpoints[1,c(1:8)]), 0.1, c(11, 0, 9.5E-4), NULL, NULL, clean=FALSE)
cmd = output8IBXIIIBc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IBXIIIBb$bifpoints[1,c(1:4,6,8)]), -0.1, c(11, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K1_IB_BP1 <- output8IBXICa$bifpoints[1,];
Eq_K1_IB_BP1b <- output8IBXIIIBa$bifpoints[1,];
Eq_K1_IB_BP2 <- output8IBXICc$bifpoints[2,];
Eq_K1_IB_BP2b <- output8IBXIIIBb$bifpoints[1,];
Eq_K2_IB_BPE <- output7IBXICb$bifpoints[1,];
Eq_K1_IB_LP <- output8IBXICc$bifpoints[1,];
Eq_K1_IB_LPb <- output8IBXIIIBc$bifpoints[2,];

### Codim ========================
### BP 1
cmd = output10IBXICa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IBXIIIBa$bifpoints[1,c(1:5,7)], 3E-4), -0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output10IBXICb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IBXIIIBa$bifpoints[1,c(1:5,7)], 3E-4), 0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
output10IBXIC<-rbind(output10IBXICa$curvepoints[dim(output10IBXICa$curvepoints)[1]:1,],
                     output10IBXICb$curvepoints[2:dim(output10IBXICb$curvepoints)[1],]);

### BP 2
cmd = output11IBXICa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IBXICc$bifpoints[2,c(1:4,6,8)], 3E-5), -0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
cmd = output11IBXICb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IBXICc$bifpoints[2,c(1:4,6,8)], 3E-5), 0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
output11IBXIC<-rbind(output11IBXICa$curvepoints[dim(output11IBXICa$curvepoints)[1]:1,],
                     output11IBXICb$curvepoints[2:dim(output11IBXICb$curvepoints)[1],]);

### BPE
cmd = output12IBXICa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IBXICb$bifpoints[,c(1:3,5:8)], 3E-5), -0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envBP", "2"), clean=FALSE)
cmd = output12IBXICb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IBXICb$bifpoints[,c(1:3,5:8)], 3E-5), 0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envBP", "2"), clean=FALSE)
output12IBXIC<-rbind(output12IBXICa$curvepoints[dim(output12IBXICa$curvepoints)[1]:1,],
                     output12IBXICb$curvepoints[2:dim(output12IBXICb$curvepoints)[1]-1,])

### LP
cmd = output13IBXICa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IBXICc$bifpoints[1,c(1:8)], 3E-5), -1.2, c(11, 0, 9.5E-01, 10, 0.1, 9.5E-01), NULL, NULL, clean=FALSE)
cmd = output13IBXICb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IBXICc$bifpoints[1,c(1:8)], 3E-5), 0.1, c(11, 0, 9.5E-01, 10, 0.1, 9.5E-01), NULL, NULL, clean=FALSE)
cmd = output13IBXICc <- PSPMequi("Tritrophicmod5_fin", "LP", c(output13IBXICa$curvepoints[dim(output13IBXICa$curvepoints)[1],c(1:8)], 3E-5), -0.3, c(11, 0, 9.5E-01, 10, 0.1, 9.5E-01), NULL, NULL, clean=FALSE)

output13IBXIC<-rbind(output13IBXICc$curvepoints[dim(output13IBXICc$curvepoints)[1]:1,],
                     output13IBXICa$curvepoints[dim(output13IBXICa$curvepoints)[1]:1,],
                     output13IBXICb$curvepoints[2:dim(output13IBXICb$curvepoints)[1]-1,])

cmd = output13IBXICc <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IBXIIIBc$bifpoints[2,c(1:4,6,8)], 3E-4), -0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output13IBXICd <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IBXIIIBc$bifpoints[2,c(1:4,6,8)], 3E-4), 0.1, c(11, 0, 9.5E-01, 10, -0.7, 9.5E-01), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
output13IBXIC2<-rbind(output13IBXICc$curvepoints[dim(output13IBXICc$curvepoints)[1]:1,],
                      output13IBXICd$curvepoints[2:dim(output13IBXICd$curvepoints)[1],]);

################## IC. BETA = 2.0 ========================
### EQ ========================
### Along K2 gradient
# Starting from C1-R1+R2 with K1=3E-5
cmd = output7ICXICa <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output2.2$curvepoints[27,c(3,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7ICXICb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IBXICa$bifpoints[1,c(1:3,5:8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("envZE", "2"), clean=FALSE)
cmd = output7ICXICc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IBXICb$bifpoints[1,c(1:8)]), -0.1, c(10, 0, 9.5E-4), NULL, NULL, clean=FALSE)

### Along K1 gradient
# Starting from P-C2-R2+R1 with K2=3E-4
cmd = output8ICXIIIBa <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output3$curvepoints[426,c(2:5,7)]), 0.1, c(11, 0, 9.5E-4), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)
cmd = output8ICXIIIBb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8ICXIIIBa$bifpoints[1,c(1:8)]), 0.1, c(11, 0, 9.5E-4), NULL, NULL, clean=FALSE)
cmd = output8ICXIIIBc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8ICXIIIBb$bifpoints[1,c(1:4,6,8)]), -0.1, c(11, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K1_IC_BP1 <- output8ICXIIIBa$bifpoints[1,];
Eq_K1_IC_BP2 <- output8ICXIIIBb$bifpoints[1,];
Eq_K2_IC_BPE <- output7ICXICb$bifpoints[1,];
Eq_K2_IC_LP <- output7ICXICc$bifpoints[1,];
Eq_K1_IC_LPb <- output8ICXIIIBc$bifpoints[1,];

### Codim ========================
### BP 1
cmd = output10ICXIIIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8ICXIIIBa$bifpoints[1,c(1:5,7)], 3E-4), -0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output10ICXIIIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8ICXIIIBa$bifpoints[1,c(1:5,7)], 3E-4), 0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
output10ICXIIIB<-rbind(output10ICXIIIBa$curvepoints[dim(output10ICXIIIBa$curvepoints)[1]:1,],
                       output10ICXIIIBb$curvepoints[2:dim(output10ICXIIIBb$curvepoints)[1],]);

### BP 2
cmd = output11ICXIIIBa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8ICXIIIBb$bifpoints[1,c(1:4,6,8)], 3E-4), -0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
cmd = output11ICXIIIBb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8ICXIIIBb$bifpoints[1,c(1:4,6,8)], 3E-4), 0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
output11ICXIIIB<-rbind(output11ICXIIIBa$curvepoints[dim(output11ICXIIIBa$curvepoints)[1]:1,],
                       output11ICXIIIBb$curvepoints[2:dim(output11ICXIIIBb$curvepoints)[1],]);

### BPE
cmd = output12ICXIIIBa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7ICXICb$bifpoints[1,c(1:3,5:8)], 3E-5), -0.1, c(10, 0, 9.5E-01, 11, -0.1, 9.5E-01), NULL, c("envBP", "2"), clean=FALSE)
cmd = output12ICXIIIBb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7ICXICb$bifpoints[1,c(1:3,5:8)], 3E-5), 0.1, c(10, 0, 9.5E-01, 11, -0.1, 9.5E-01), NULL, c("envBP", "2"), clean=FALSE)
output12ICXIIIB<-rbind(output12ICXIIIBa$curvepoints[dim(output12ICXIIIBa$curvepoints)[1]:1,],
                       output12ICXIIIBb$curvepoints[2:dim(output12ICXIIIBb$curvepoints)[1],]);

### LP
cmd = output13ICXIIIBa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7ICXICc$bifpoints[1,c(1:8)], 3E-5), -0.1, c(10, 0, 9.5E-01, 11, -0.1, 9.5E-01), NULL, NULL, clean=FALSE)
cmd = output13ICXIIIBb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7ICXICc$bifpoints[1,c(1:8)], 3E-5), 0.1, c(10, 0, 9.5E-01, 11, -0.7, 9.5E-01), NULL, NULL, clean=FALSE)
output13ICXIIIB<-rbind(output13ICXIIIBa$curvepoints[dim(output13ICXIIIBa$curvepoints)[1]:1,],
                       output13ICXIIIBb$curvepoints[2:dim(output13ICXIIIBb$curvepoints)[1],]);

cmd = output13ICXIIIBc <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8ICXIIIBc$bifpoints[1,c(1:4,6,8)], 3E-4), -0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output13ICXIIIBd <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8ICXIIIBc$bifpoints[1,c(1:4,6,8)], 3E-4), 0.1, c(11, 0, 9.5E-01, 10, -0.7, 9.5E-01), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
output13ICXIIIB2<-rbind(output13ICXIIIBc$curvepoints[dim(output13ICXIIIBc$curvepoints)[1]:1,],
                        output13ICXIIIBd$curvepoints[2:dim(output13ICXIIIBd$curvepoints)[1],]);

################## IIA. BETA = 2.5 ========================
### EQ ========================
### Along K2 gradient
# Starting from C1-R1+R2 with K1=3E-5
cmd = output7IIAXICa <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output2.2$curvepoints[27,c(3,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7IIAXICb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIAXICa$bifpoints[1,c(1:3,5:8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("envZE", "2"), clean=FALSE)
cmd = output7IIAXICc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIAXICb$bifpoints[1,c(1:8)]), -0.1, c(10, 0, 9.5E-4), NULL, NULL, clean=FALSE)

### Along K1 gradient
# Starting from P-C2-R2+R1 with K2=3E-4
cmd = output8IIAXIIIBa <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output3$curvepoints[426,c(2:5,7)]), 0.1, c(11, 0, 9.5E-4), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)
cmd = output8IIAXIIIBb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IIAXIIIBa$bifpoints[1,c(1:8)]), 0.1, c(11, 0, 9.5E-4), NULL, NULL, clean=FALSE)
cmd = output8IIAXIIIBc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IIAXIIIBb$bifpoints[1,c(1:4,6,8)]), -0.1, c(11, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K1_IIA_BP1 <- output8IIAXIIIBa$bifpoints[1,];
Eq_K1_IIA_BP2 <- output8IIAXIIIBb$bifpoints[1,];
Eq_K2_IIA_BPE <- output7IIAXICb$bifpoints[1,];
Eq_K2_IIA_LP <- output7IIAXICc$bifpoints[2,];
Eq_K1_IIA_LPb <- output8IIAXIIIBc$bifpoints[1,];

### Codim ========================
### BP 1
cmd = output10IIAXICa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IIAXIIIBa$bifpoints[1,c(1:5,7)], 3E-4), -0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output10IIAXICb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IIAXIIIBa$bifpoints[1,c(1:5,7)], 3E-4), 0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
output10IIAXIC<-rbind(output10IIAXICa$curvepoints[dim(output10IIAXICa$curvepoints)[1]:1,],
                      output10IIAXICb$curvepoints[2:dim(output10IIAXICb$curvepoints)[1],]);

### BP 2
cmd = output11IIAXICa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IIAXIIIBb$bifpoints[1,c(1:4,6,8)], 3E-4), -0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
cmd = output11IIAXICb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IIAXIIIBb$bifpoints[1,c(1:4,6,8)], 3E-4), 0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
output11IIAXIC<-rbind(output11IIAXICa$curvepoints[dim(output11IIAXICa$curvepoints)[1]:1,],
                      output11IIAXICb$curvepoints[2:dim(output11IIAXICb$curvepoints)[1],]);

### BPE
cmd = output12IIAXICa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IIAXICb$bifpoints[1,c(1:3,5:8)], 3E-5), -0.1, c(10, 0, 9.5E-01, 11, -0.1, 9.5E-01), NULL, c("envBP", "2"), clean=FALSE)
cmd = output12IIAXICb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IIAXICb$bifpoints[1,c(1:3,5:8)], 3E-5), 0.1, c(10, 0, 9.5E-01, 11, -0.1, 9.5E-01), NULL, c("envBP", "2"), clean=FALSE)
output12IIAXIC<-rbind(output12IIAXICa$curvepoints[dim(output12IIAXICa$curvepoints)[1]:1,],
                      output12IIAXICb$curvepoints[2:dim(output12IIAXICb$curvepoints)[1],]);

### LP
cmd = output13IIAXICa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIAXICc$bifpoints[2,c(1:8)], 3E-5), -0.1, c(10, 0, 9.5E-01, 11, -0.1, 9.5E-01), NULL, NULL, clean=FALSE)
cmd = output13IIAXICb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIAXICc$bifpoints[2,c(1:8)], 3E-5), 0.1, c(10, 0, 9.5E-01, 11, -0.1, 9.5E-01), NULL, NULL, clean=FALSE)
output13IIAXIC<-rbind(output13IIAXICa$curvepoints[dim(output13IIAXICa$curvepoints)[1]:1,],
                      output13IIAXICb$curvepoints[2:dim(output13IIAXICb$curvepoints)[1],]);

cmd = output13IIAXICc <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IIAXIIIBc$bifpoints[1,c(1:4,6,8)], 3E-4), -0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output13IIAXICd <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IIAXIIIBc$bifpoints[1,c(1:4,6,8)], 3E-4), 0.1, c(11, 0, 9.5E-01, 10, -0.7, 9.5E-01), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
output13IIAXIC2<-rbind(output13IIAXICc$curvepoints[dim(output13IIAXICc$curvepoints)[1]:1,],
                       output13IIAXICd$curvepoints[2:dim(output13IIAXICd$curvepoints)[1],]);

################## IIB. BETA = 3.0 ========================
### EQ ========================
### Along K2 gradient
# Starting from C1-R1+R2 with K1=3E-5
cmd = output7IIBXICa <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output2.2$curvepoints[27,c(3,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7IIBXICb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIBXICa$bifpoints[1,c(1:3,5:8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("envZE", "2"), clean=FALSE)
cmd = output7IIBXICc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IIBXICb$bifpoints[1,c(1:8)]), -0.1, c(10, 0, 9.5E-4), NULL, NULL, clean=FALSE)

### Along K1 gradient
# Starting from P-C2-R2+R1 with K2=3E-4
cmd = output8IIBXIIIBa <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output3$curvepoints[426,c(2:5,7)]), 0.1, c(11, 0, 9.5E-4), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)
cmd = output8IIBXIIIBb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IIBXIIIBa$bifpoints[1,c(1:8)]), 0.1, c(11, 0, 9.5E-2), NULL, NULL, clean=FALSE)
cmd = output8IIBXIIIBc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IIBXIIIBb$bifpoints[1,c(1:4,6,8)]), -0.1, c(11, 0, 9.5E-2), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K1_IIB_BP1 <- output8IIBXIIIBa$bifpoints[1,];
Eq_K1_IIB_BP2 <- output8IIBXIIIBb$bifpoints[1,];
Eq_K2_IIB_BPE <- output7IIBXICb$bifpoints[1,];
Eq_K2_IIB_LP <- output7IIBXICc$bifpoints[2,];
Eq_K1_IIB_LPb <- output8IIBXIIIBc$bifpoints[2,];
### Codim ========================
### BP 1
cmd = output10IIBXICa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IIBXIIIBa$bifpoints[1,c(1:5,7)], 3E-4), -0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output10IIBXICb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IIBXIIIBa$bifpoints[1,c(1:5,7)], 3E-4), 0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
output10IIBXIC<-rbind(output10IIBXICa$curvepoints[dim(output10IIBXICa$curvepoints)[1]:1,],
                      output10IIBXICb$curvepoints[2:dim(output10IIBXICb$curvepoints)[1],]);

### BP 2
cmd = output11IIBXICa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IIBXIIIBb$bifpoints[1,c(1:4,6,8)], 3E-4), -0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
cmd = output11IIBXICb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IIBXIIIBb$bifpoints[1,c(1:4,6,8)], 3E-4), 0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
output11IIBXIC<-rbind(output11IIBXICa$curvepoints[dim(output11IIBXICa$curvepoints)[1]:1,],
                      output11IIBXICb$curvepoints[2:dim(output11IIBXICb$curvepoints)[1],]);

### BPE
cmd = output12IIBXICa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IIBXICb$bifpoints[1,c(1:3,5:8)], 3E-5), -0.1, c(10, 0, 9.5E-01, 11, -0.1, 9.5E-01), NULL, c("envBP", "2"), clean=FALSE)
cmd = output12IIBXICb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IIBXICb$bifpoints[1,c(1:3,5:8)], 3E-5), 0.1, c(10, 0, 9.5E-01, 11, -0.1, 9.5E-01), NULL, c("envBP", "2"), clean=FALSE)
output12IIBXIC<-rbind(output12IIBXICa$curvepoints[dim(output12IIBXICa$curvepoints)[1]:1,],
                      output12IIBXICb$curvepoints[2:dim(output12IIBXICb$curvepoints)[1],]);

### LP
cmd = output13IIBXICa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIBXICc$bifpoints[2,c(1:8)], 3E-5), -0.1, c(10, 0, 9.5E-01, 11, -0.1, 9.5E-01), NULL, NULL, clean=FALSE)
cmd = output13IIBXICb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IIBXICc$bifpoints[2,c(1:8)], 3E-5), 0.1, c(10, 0, 9.5E-01, 11, -0.1, 9.5E-01), NULL, NULL, clean=FALSE)
output13IIBXIC<-rbind(output13IIBXICa$curvepoints[dim(output13IIBXICa$curvepoints)[1]:1,],
                      output13IIBXICb$curvepoints[2:dim(output13IIBXICb$curvepoints)[1],]);

cmd = output13IIBXICc <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IIBXIIIBc$bifpoints[2,c(1:4,6,8)], 3E-4), -0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output13IIBXICd <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IIBXIIIBc$bifpoints[2,c(1:4,6,8)], 3E-4), 0.1, c(11, 0, 9.5E-01, 10, -0.7, 9.5E-01), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
output13IIBXIC2<-rbind(output13IIBXICc$curvepoints[dim(output13IIBXICc$curvepoints)[1]:1,],
                       output13IIBXICd$curvepoints[2:dim(output13IIBXICd$curvepoints)[1],]);
################## IIC. BETA = 3.5 ========================
### EQ ========================
### Along K2 gradient
# Starting from C1-R1+R2 with K1=3E-5
cmd = output7IICXICa <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output2.2$curvepoints[27,c(3,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = output7IICXICb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IICXICa$bifpoints[1,c(1:3,5:8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("envZE", "2"), clean=FALSE)
cmd = output7IICXICc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output7IICXICb$bifpoints[1,c(1:8)]), -0.1, c(10, 0, 9.5E-4), NULL, NULL, clean=FALSE)

### Along K1 gradient
# Starting from P-C2-R2+R1 with K2=3E-4
cmd = output8IICXIIIBa <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output3$curvepoints[426,c(2:5,7)]), 0.1, c(11, 0, 9.5E-4), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)
cmd = output8IICXIIIBb <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IICXIIIBa$bifpoints[1,c(1:8)]), 0.1, c(11, 0, 9.5E-2), NULL, NULL, clean=FALSE)
cmd = output8IICXIIIBc <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output8IICXIIIBb$bifpoints[1,c(1:4,6,8)]), -0.1, c(11, 0, 9.5E-2), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K1_IIC_BP1 <- output8IICXIIIBa$bifpoints[1,];
Eq_K1_IIC_BP2 <- output8IICXIIIBb$bifpoints[1,];
Eq_K2_IIC_BPE <- output7IICXICb$bifpoints[1,];
Eq_K2_IIC_LP <- output7IICXICc$bifpoints[2,];
Eq_K1_IIC_LPb <- output8IICXIIIBc$bifpoints[2,];
### Codim ========================
### BP 1
cmd = output10IICXICa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IICXIIIBa$bifpoints[1,c(1:5,7)], 3E-4), -0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output10IICXICb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IICXIIIBa$bifpoints[1,c(1:5,7)], 3E-4), 0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "4", "popBP", "1"), clean=FALSE)
output10IICXIC<-rbind(output10IICXICa$curvepoints[dim(output10IICXICa$curvepoints)[1]:1,],
                      output10IICXICb$curvepoints[2:dim(output10IICXICb$curvepoints)[1],]);

### BP 2
cmd = output11IICXICa <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IICXIIIBb$bifpoints[1,c(1:4,6,8)], 3E-4), -0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
cmd = output11IICXICb <- PSPMequi("Tritrophicmod5_fin", "BP", c(output8IICXIIIBb$bifpoints[1,c(1:4,6,8)], 3E-4), 0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("envZE", "3", "popBP", "0"), clean=FALSE)
output11IICXIC<-rbind(output11IICXICa$curvepoints[dim(output11IICXICa$curvepoints)[1]:1,],
                      output11IICXICb$curvepoints[2:dim(output11IICXICb$curvepoints)[1],]);

### BPE
cmd = output12IICXICa <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IICXICb$bifpoints[1,c(1:3,5:8)], 3E-5), -0.1, c(10, 0, 9.5E-01, 11, -0.1, 9.5E-01), NULL, c("envBP", "2"), clean=FALSE)
cmd = output12IICXICb <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output7IICXICb$bifpoints[1,c(1:3,5:8)], 3E-5), 0.1, c(10, 0, 9.5E-01, 11, -0.1, 9.5E-01), NULL, c("envBP", "2"), clean=FALSE)
output12IICXIC<-rbind(output12IICXICa$curvepoints[dim(output12IICXICa$curvepoints)[1]:1,],
                      output12IICXICb$curvepoints[2:dim(output12IICXICb$curvepoints)[1],]);

### LP
cmd = output13IICXICa <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IICXICc$bifpoints[2,c(1:8)], 3E-5), -0.1, c(10, 0, 9.5E-01, 11, -0.1, 9.5E-01), NULL, NULL, clean=FALSE)
cmd = output13IICXICb <- PSPMequi("Tritrophicmod5_fin", "LP", c(output7IICXICc$bifpoints[2,c(1:8)], 3E-5), 0.1, c(10, 0, 9.5E-01, 11, -0.1, 9.5E-01), NULL, NULL, clean=FALSE)
output13IICXIC<-rbind(output13IICXICa$curvepoints[dim(output13IICXICa$curvepoints)[1]:1,],
                      output13IICXICb$curvepoints[2:dim(output13IICXICb$curvepoints)[1],]);

cmd = output13IICXICc <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IICXIIIBc$bifpoints[2,c(1:4,6,8)], 3E-4), -0.1, c(11, 0, 9.5E-01, 10, -0.1, 9.5E-01), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = output13IICXICd <- PSPMequi("Tritrophicmod5_fin", "LP", c(output8IICXIIIBc$bifpoints[2,c(1:4,6,8)], 3E-4), 0.1, c(11, 0, 9.5E-01, 10, -0.7, 9.5E-01), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
output13IICXIC2<-rbind(output13IICXICc$curvepoints[dim(output13IICXICc$curvepoints)[1]:1,],
                       output13IICXICd$curvepoints[2:dim(output13IICXICd$curvepoints)[1],]);
