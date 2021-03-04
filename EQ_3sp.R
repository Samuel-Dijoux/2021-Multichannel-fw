###########################################################################################
#               -> Equilibrium analyses of tritrophic chain "P-C1-R1" & "P-C2-R2"
#                     along productivity gradients
#
#               -> Bifurcation analyses over K1 and Beta gradients for P-C1-R1
###########################################################################################
library(PSPManalysis)
if (!exists("PSPMsrc.fullpath")) source("../PSPMworkshop/PSPManalysis/PSPManalysis.r")

################## Computation of P-C1-R1  ================================================

cat("\n\n\nStarting from the (known) trivial equilibrium with only resources, R1\n\n")
cmd = output1.2 <- PSPMequi("Tritrophicmod5_fin", "EQ", c(1.0E-08, 1.0E-08, 1.0E-08), 0.5, c(11, 0, 9e-2), NULL,c("popZE", "0", "popZE", "1", "envZE", "2", "envZE", "3", "envZE", "4"), clean=FALSE)

str = readline(paste0("\n> ", cmd))
if ((str == '') || (str == 'y') || (str == 'Y')) eval(parse(text=cmd)) else stop('Test script interrupted by user')

cat("\n\n\nStart from the detected branching point to compute the C1-R1 equilibrium\n\n")
cmd = output2.2 <- PSPMequi("Tritrophicmod5_fin", "EQ", output1.2$bifpoint[,c(1, 2, 3, 8)], 0.5, c(11, 0, 9e-2), NULL,c("popZE", "0", "envZE", "2", "envZE", "3", "envZE", "4"), clean=FALSE)

str = readline(paste0("\n> ", cmd))
if ((str == '') || (str == 'y') || (str == 'Y')) eval(parse(text=cmd)) else stop('Test script interrupted by user')

cat("\n\n\nStart from the detected branching point of the predator to compute the (P)-C1-R1 equilibrium\n\n")
cmd = output3.2 <- PSPMequi("Tritrophicmod5_fin", "EQ", output2.2$bifpoint[1,c(1,2,3,4,6,8)], -0.1, c(11, 0, 9e-3), NULL,c("popZE", "0", "envZE", "3"), clean=FALSE)

str = readline(paste0("\n> ", cmd))
if ((str == '') || (str == 'y') || (str == 'Y')) eval(parse(text=cmd)) else stop('Test script interrupted by user')

cmd = output3.2b <- PSPMequi("Tritrophicmod5_fin", "EQ", output2.2$bifpoint[1,c(1,2,3,4,6,8)], 0.1, c(11, 0, 9e-2), NULL,c("popZE", "0", "envZE", "3"), clean=FALSE)
str = readline(paste0("\n> ", cmd))
if ((str == '') || (str == 'y') || (str == 'Y')) eval(parse(text=cmd)) else stop('Test script interrupted by user')

################## Computation of P-C2-R2  ================================================

cat("\n\n\nStarting from the (known) trivial equilibrium with only resources, R2\n\n")
cmd = output1 <- PSPMequi("Tritrophicmod5_fin", "EQ", c(1.0E-06, 1.0E-06, 1.0E-06), 0.5, c(10, 0, 9e-2), NULL,c("popZE", "0", "popZE", "1", "envZE", "2", "envZE", "3", "envZE", "4"), clean=FALSE)

str = readline(paste0("\n> ", cmd))
if ((str == '') || (str == 'y') || (str == 'Y')) eval(parse(text=cmd)) else stop('Test script interrupted by user')

cat("\n\n\nStart from the detected branching point to compute the C2-R2 equilibrium\n\n")
cmd = output2 <- PSPMequi("Tritrophicmod5_fin", "EQ", output1$bifpoint[,c(1, 2, 3, 7)], 0.5, c(10, 0, 9e-2), NULL,c("popZE", "1", "envZE", "2", "envZE", "3", "envZE", "4"), clean=FALSE)

str = readline(paste0("\n> ", cmd))
if ((str == '') || (str == 'y') || (str == 'Y')) eval(parse(text=cmd)) else stop('Test script interrupted by user')

cat("\n\n\nStart from the detected branching point of the predator to compute the (P)-C2-R2 equilibrium\n\n")
cmd = output3 <- PSPMequi("Tritrophicmod5_fin", "EQ", output2$bifpoint[,c(1,2,3,5,12,7)], -0.1, c(10, 0, 9e-2), NULL,c("popZE", "1", "envZE", "4"), clean=FALSE)

str = readline(paste0("\n> ", cmd))
if ((str == '') || (str == 'y') || (str == 'Y')) eval(parse(text=cmd)) else stop('Test script interrupted by user')

cmd = 'output3b <- PSPMequi("Tritrophicmod5_fin", "EQ", output2$bifpoint[,c(1,2,3,5,12,7)], 0.1, c(10, 0, 9e-2), NULL,c("popZE", "1", "envZE", "4"), clean=FALSE)'
str = readline(paste0("\n> ", cmd))
if ((str == '') || (str == 'y') || (str == 'Y')) eval(parse(text=cmd)) else stop('Test script interrupted by user')

################## Bifurcation analyses P-C1-R1 ===========================================
###  BP-curve Computation  ========================
cmd = output4 <- PSPMequi("Tritrophicmod5_fin", "BP", c(output1.2$bifpoints[1:3], 1.2), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4), NULL,
                          c("popZE", "0", "envZE", "2", "envZE", "3", "envZE", "4", "popBP", "1"), clean=FALSE)
cmd = output4b <- PSPMequi("Tritrophicmod5_fin", "BP", c(output1.2$bifpoints[1:3], 1.2), 0.1, c(11, 0, 9.5E-01, 20, 0.1, 4), NULL,
                           c("popZE", "0", "envZE", "2", "envZE", "3", "envZE", "4", "popBP", "1"), clean=FALSE)

o_4<-rbind(output4$curvepoints[dim(output4$curvepoints)[1]:1,],
           output4b$curvepoints[2:dim(output4b$curvepoints)[1],]);

###  BPE-curve Computation  ========================
cmd = output5 <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output2.2$bifpoints[1,c(1:3, 8)], 1.2), 0.1, c(11, 0, 9.5E-01, 20, 0.1, 4), NULL,
                          c("popZE", "0", "envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)
cmd = output5b <- PSPMequi("Tritrophicmod5_fin", "BPE", c(output2.2$bifpoints[1,c(1:3, 8)], 1.2), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4), NULL,
                           c("popZE", "0", "envZE", "3", "envZE", "4", "envBP", "2"), clean=FALSE)

o_5<-rbind(output5$curvepoints[dim(output5$curvepoints)[1]:1,],
           output5b$curvepoints[2:dim(output5b$curvepoints)[1],]);

###  LP-curve Computation  ========================
cmd = output6 <- PSPMequi("Tritrophicmod5_fin", "LP", c(output3.2$bifpoints[2,c(1:4,6,8)], 1.2), 0.5, c(11, 0, 9.5E-01, 20, 0.1, 4.0), NULL,
                          c("popZE","0", "envZE","3"),clean=FALSE)
cmd = output6b <- PSPMequi("Tritrophicmod5_fin", "LP", c(output3.2$bifpoints[2,c(1:4,6,8)], 1.2), -0.9, c(11, 0, 9.5E-01, 20, 0.1, 4.0), NULL,
                           c("popZE","0", "envZE","3"),clean=FALSE)

o_6<-rbind(output6b$curvepoints[dim(output6b$curvepoints)[1]:1,],
           output6$curvepoints[2:dim(output6$curvepoints)[1],]);
