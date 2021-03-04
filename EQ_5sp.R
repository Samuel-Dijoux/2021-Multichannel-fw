###########################################################################################
#                -> Equilibrium analyses of the Multichannel food web "P-C1C2-R1R2"
#
# Selected analyses Equilibrium analyses   
# Invasion of C1 in R2-C2-P+R1 over K1 gradient for varied K2 and Beta
# Invasion of C2 in R1-C1-P+R2 over K2 gradient for varied K1 and Beta
###########################################################################################
library(PSPManalysis)
if (!exists("PSPMsrc.fullpath")) source("../PSPMworkshop/PSPManalysis/PSPManalysis.r")

################## Equilibrium analyses over K1 gradient for varied K2 and Beta ===========
############## K1-0 , K2=5E-6 + Beta = 1.2 ========================
cmd = o_K1_0a <- PSPMequi("Tritrophicmod5_fin", "EQ", c(1.0E-08, 1.0E-08, 1.0E-08), 0.5, c(11, 0, 4E-03), NULL,
                          c("popZE", "0", "popZE", "1", "envZE", "2", "envZE", "3", "envZE", "4"), clean=FALSE)
cmd = o_K1_0b <- PSPMequi("Tritrophicmod5_fin", "EQ", o_K1_0a$bifpoint[,c(1, 2, 3, 8)], 0.2, c(11, 0, 4E-03), NULL,
                          c("popZE", "0", "envZE", "2", "envZE", "3", "envZE", "4"), clean=FALSE)
cmd = o_K1_0c <- PSPMequi("Tritrophicmod5_fin", "EQ", o_K1_0b$bifpoint[,c(1,2,3,5,17,8)], -0.1, c(11, 0, 4E-03), NULL,
                          c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K1_0_BP1 <- o_K1_0a$bifpoints[1,];
Eq_K1_0_BPE <- o_K1_0b$bifpoints[1,];
Eq_K1_0_LP  <- o_K1_0c$bifpoints[2,];
############## K1-A, K2=1E-5 + Beta = 1.2 ========================
cmd = o_K1_1a <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[7,2], 1E-8, output2$curvepoints[7,c(5,7)]),
                          0.1, c(11, 0, 4E-4), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = o_K1_1b <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_1a$bifpoints[1,c(1:3,5:8)]),
                          0.1, c(11, 0, 5E-3), NULL, c("envZE", "2"), clean=FALSE)
cmd = o_K1_1c <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_1b$bifpoints[,c(1:8)]),
                          -0.5, c(11, 0, 5E-3), NULL, NULL, clean=FALSE)
# After extinction of C2
cmd = o_K1_1d <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_1c$bifpoints[1,c(1:4,6,8)]),
                          -0.7, c(11, 0, 5E-3), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K1_A_BP1 <- o_K1_1a$bifpoints[1,]; #inv of C1 in C2-R2+R1
Eq_K1_A_BP2 <- o_K1_1c$bifpoints[1,];
Eq_K1_A_BPE <- o_K1_1b$bifpoints[1,];
Eq_K1_A_LP <- o_K1_1d$bifpoints[1,];
############## K1-B, K2=3E-5 + Beta = 1.2 ========================
cmd = o_K1_2a <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[21,2], 1E-8, output2$curvepoints[21,c(5,7)]),
                          0.1, c(11, 0, 4E-4), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = o_K1_2b <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_2a$bifpoints[1,c(1:3,5:8)]),
                          0.1, c(11, 0, 5E-3), NULL, c("envZE", "2"), clean=FALSE)
cmd = o_K1_2c <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_2b$bifpoints[,c(1:8)]),
                          -0.5, c(11, 0, 5E-3), NULL, NULL, clean=FALSE)
# After extinction of C2
cmd = o_K1_2d <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_2c$bifpoints[2,c(1:4,6,8)]),
                          0.1, c(11, 0, 5E-3), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K1_B_BP1 <- o_K1_2a$bifpoints[1,]; #inv of C1 in C2-R2+R1
Eq_K1_B_BP2 <- o_K1_2c$bifpoints[2,];
Eq_K1_B_BPE <- o_K1_2b$bifpoints[1,];
Eq_K1_B_LP <- o_K1_2c$bifpoints[1,];
############## K1-C, K2=3E-5 + Beta = 2 ========================
cmd = o_K1_3a <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[27,2], 1E-8, output2$curvepoints[27,c(5,7)]),
                          0.1, c(11, 0, 4E-4), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = o_K1_3b <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_3a$bifpoints[1,c(1:3,5:8)]),
                          0.1, c(11, 0, 1E-3), NULL, c("envZE", "2"), clean=FALSE)
cmd = o_K1_3c <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_3b$bifpoints[,c(1:8)]),
                          -0.6, c(11, 0, 1E-3), NULL, NULL, clean=FALSE)
# After extinction of C2
cmd = o_K1_3d <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_3c$bifpoints[2,c(1:4,6,8)]),
                          -0.1, c(11, 0, 5E-3), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K1_C_BP1 <- o_K1_3a$bifpoints[1,]; #inv of C1 in C2-R2+R1
Eq_K1_C_BP2 <- o_K1_3c$bifpoints[2,];
Eq_K1_C_BPE <- o_K1_3b$bifpoints[1,];
Eq_K1_C_LP <- o_K1_3d$bifpoints[2,];
############## K1-D, K2=3E-5 + Beta = 3.5  ========================
cmd = o_K1_4a <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[27,2], 1E-8, output2$curvepoints[27,c(5,7)]),
                          0.1, c(11, 0, 4E-4), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = o_K1_4b <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_4a$bifpoints[1,c(1:3,5:8)]),
                          0.1, c(11, 0, 5E-2), NULL, c("envZE", "2"), clean=FALSE)
cmd = o_K1_4c <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_4b$bifpoints[,c(1:8)]),
                          -0.6, c(11, 0, 5E-2), NULL, NULL, clean=FALSE)
# After extinction of C2
cmd = o_K1_4d <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_4c$bifpoints[3,c(1:4,6,8)]),
                          -0.1, c(11, 0, 5E-2), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K1_D_BP1 <- o_K1_4a$bifpoints[1,]; #inv of C1 in C2-R2+R1
Eq_K1_D_BP2 <- o_K1_4c$bifpoints[3,];
Eq_K1_D_BPE <- o_K1_4b$bifpoints[1,];
Eq_K1_D_LP <- o_K1_4c$bifpoints[2,];
Eq_K1_D_LP2 <- o_K1_4d$bifpoints[1,];
############## K1-E, K2=8E-5 + Beta = 2 ========================
cmd = o_K1_5a <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[52,2], 1E-8, output2$curvepoints[52,c(5,7)]),
                          0.1, c(11, 0, 4E-4), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = o_K1_5b <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_5a$bifpoints[1,c(1:3,5:8)]),
                          0.1, c(11, 0, 5E-2), NULL, c("envZE", "2"), clean=FALSE)
cmd = o_K1_5c <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_5b$bifpoints[,c(1:8)]),
                          -0.6, c(11, 0, 5E-2), NULL, NULL, clean=FALSE)
# After extinction of C2
cmd = o_K1_5d <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_5c$bifpoints[3,c(1:4,6,8)]),
                          -0.1, c(11, 0, 5E-2), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K1_E_BP1 <- o_K1_5a$bifpoints[1,]; #inv of C1 in C2-R2+R1
Eq_K1_E_BP2 <- o_K1_5c$bifpoints[3,];
Eq_K1_E_BPE <- o_K1_5b$bifpoints[1,];
Eq_K1_E_LP <- o_K1_5c$bifpoints[2,];
Eq_K1_E_LP2 <- o_K1_5d$bifpoints[2,];
############## K1-F, K2=1E-4 + Beta = 0.5 ========================
# Starting from C2-R2+R1
cmd = o_K1_6a <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[52,2], 1E-8, output2$curvepoints[52,c(5,7)]),
                          0.1, c(11, 0, 4E-4), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = o_K1_6b <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_6a$bifpoints[1,c(1:3,5:8)]),
                          0.1, c(11, 0, 5E-2), NULL, c("envZE", "2"), clean=FALSE)
cmd = o_K1_6c <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_6b$bifpoints[,c(1:8)]),
                          -0.6, c(11, 0, 5E-2), NULL, NULL, clean=FALSE)
# After extinction of C1
cmd = o_K1_6d <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_6c$bifpoints[2,c(1:4,5,7)]),
                          -0.1, c(11, 0, 5E-2), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)
cmd = o_K1_6e <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_6c$bifpoints[2,c(1:4,5,7)]),
                          0.1, c(11, 0, 5E-2), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

#Starting with from P-C2-R2+R1
cmd = o_K1_6f <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output3$curvepoints[276,2], 1E-8, output3$curvepoints[276,c(4,5,7)]),
                          0.1, c(11, 0, 1.1E-3), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

Eq_K1_F_BP1 <- o_K1_6a$bifpoints[1,]; #inv of C1 in C2-R2+R1
Eq_K1_F_BP1b <- o_K1_6c$bifpoints[2,];
Eq_K1_F_BPE <- o_K1_6b$bifpoints[1,];
Eq_K1_F_LP <- o_K1_6c$bifpoints[1,];
############## K1-G, K2=1E-4 + Beta = 1 ========================
# Starting from C2-R2+R1
cmd = o_K1_7a <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[62,2], 1E-8, output2$curvepoints[62,c(5,7)]),
                          0.1, c(11, 0, 4E-4), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = o_K1_7b <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_7a$bifpoints[1,c(1:3,5:8)]),
                          0.1, c(11, 0, 1E-3), NULL, c("envZE", "2"), clean=FALSE)
cmd = o_K1_7c <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_7b$bifpoints[,c(1:8)]),
                          -2.99, c(11, 0, 1E-3), NULL, NULL, clean=FALSE)
# After extinction of C2
cmd = o_K1_7d <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_7c$bifpoints[4,c(1:5,7)]),
                          0.1, c(11, 0, 1E-3), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

#Starting with from P-C2-R2+R1
cmd = o_K1_7e <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output3$curvepoints[275,2], 1E-8, output3$curvepoints[275,c(4,5,7)]),
                          0.1, c(11, 0, 1.1E-3), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

Eq_K1_G_BP1 <- o_K1_7a$bifpoints[1,]; #inv of C1 in C2-R2+R1
Eq_K1_G_BP1b <- o_K1_7c$bifpoints[1,];#unstable branch
Eq_K1_G_BP1c <- o_K1_7c$bifpoints[3,];
Eq_K1_G_BP2 <- o_K1_7c$bifpoints[4,];
Eq_K1_G_BPE <- o_K1_7b$bifpoints[1,];
Eq_K1_G_LP <- o_K1_7c$bifpoints[2,];
############## K1-H, K2=1E-4 + Beta = 2 ====
# Starting from C2-R2+R1
cmd = o_K1_8a <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[62,2], 1E-8, output2$curvepoints[62,c(5,7)]),
                          0.1, c(11, 0, 4E-4), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = o_K1_8b <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_8a$bifpoints[1,c(1:3,5:8)]),
                          0.1, c(11, 0, 2E-3), NULL, c("envZE", "2"), clean=FALSE)
cmd = o_K1_8c <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_8b$bifpoints[,c(1:8)]),
                          -0.1, c(11, 0, 2E-3), NULL, NULL, clean=FALSE)
# After extinction of C1
cmd = o_K1_8d <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_8c$bifpoints[2,c(1:8)]),
                          -0.3, c(11, 0, 2E-3), NULL, NULL, clean=FALSE)
# Starting from P-C2-R2+R1
cmd = o_K1_8e <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output3$curvepoints[277,2], 1E-8, output3$curvepoints[277,c(4,5,7)]),
                          0.1, c(11, 0, 2E-3), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)
cmd = o_K1_8f <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_8e$bifpoints[1,c(1:8)]),
                          0.1, c(11, 0, 2E-3), NULL, NULL, clean=FALSE)
# After extinction of C2
cmd = o_K1_8g <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_8f$bifpoints[1,c(1:4,6,8)]),
                          -0.1, c(11, 0, 2E-3), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K1_H_BP1 <- o_K1_8a$bifpoints[1,]; #inv of C1 in C2-R2+R1
Eq_K1_H_BP1b <- o_K1_8c$bifpoints[2,];#unstable branch
Eq_K1_H_BP1c <- o_K1_8e$bifpoints[1,];#inv of C1 in P-C2-R2+R1
Eq_K1_H_BP2 <- o_K1_8f$bifpoints[1,];
Eq_K1_H_BPE <- o_K1_8b$bifpoints[1,];
Eq_K1_H_LP <- o_K1_8g$bifpoints[1,];
############## K1-I, K2=2E-4 + Beta = 1 ========================
# Starting from C2-R2+R1
cmd = o_K1_9a <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[78,2], 1E-8, output2$curvepoints[78,c(5,7)]),
                          0.1, c(11, 0, 4E-4), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = o_K1_9b <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_9a$bifpoints[1,c(1:3,5:8)]),
                          0.1, c(11, 0, 2E-3), NULL, c("envZE", "2"), clean=FALSE)
cmd = o_K1_9c <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_9b$bifpoints[,c(1:8)]),
                          -0.6, c(11, 0, 2E-3), NULL, NULL, clean=FALSE)
# After extinction of C1
cmd = o_K1_9d <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_9c$bifpoints[,c(1:5,7)]),
                          -0.1, c(11, 0, 2E-3), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

# Starting from P-C2-R2+R1
cmd = o_K1_9e <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output3$curvepoints[372,2], 1E-8, output3$curvepoints[372,c(4,5,7)]),
                          0.1, c(11, 0, 4E-4), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)
cmd = o_K1_9f <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_9e$bifpoints[1,c(1:8)]),
                          0.1, c(11, 0, 2E-3), NULL, NULL, clean=FALSE)
# After extinction of C2
cmd = o_K1_9g <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_9f$bifpoints[1,c(1:4,6,8)]),
                          0.1, c(11, 0, 2E-3), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K1_I_BP1 <- o_K1_9a$bifpoints[1,]; #inv of C1 in C2-R2+R1
Eq_K1_I_BP1b <- o_K1_9c$bifpoints[,]; #unstable branch
Eq_K1_I_BP1c <- o_K1_9e$bifpoints[1,]; #inv of C1 in P-C2-R2+R1
Eq_K1_I_BP2 <- o_K1_9f$bifpoints[1,];
Eq_K1_I_BPE <- o_K1_9b$bifpoints[1,];
############## K1-J, K2=2E-4 + Beta = 3 ========================
# Starting from C2-R2+R1
cmd = o_K1_10a <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output2$curvepoints[78,2], 1E-8, output2$curvepoints[78,c(5,7)]),
                           0.1, c(11, 0, 4E-4), NULL, c("popZE", "1", "envZE", "2", "envZE", "4"), clean=FALSE)
cmd = o_K1_10b <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_10a$bifpoints[1,c(1:3,5:8)]),
                           0.1, c(11, 0, 1E-2), NULL, c("envZE", "2"), clean=FALSE)
cmd = o_K1_10c <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_10b$bifpoints[1,c(1:8)]),
                           -0.1, c(11, 0, 1E-2), NULL, NULL, clean=FALSE)

cmd = o_K1_10d <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_10c$bifpoints[2,c(1:5,7)]),
                           0.1, c(11, 0, 1E-2), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)
# Starting from P-C2-R2+R1
cmd = o_K1_10e <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output3$curvepoints[371,2], 1E-8, output3$curvepoints[371,c(4,5,7)]),
                           0.1, c(11, 0, 5E-2), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)
cmd = o_K1_10f <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_10e$bifpoints[1,c(1:8)]),
                           0.1, c(11, 0, 5E-2), NULL, NULL, clean=FALSE)
# After extinction of C2
cmd = o_K1_10g <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_10f$bifpoints[,c(1:4,6,8)]),
                           -0.1, c(11, 0, 1E-2), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K1_J_BP1 <- o_K1_10a$bifpoints[1,]; #inv of C1 in C2-R2+R1
Eq_K1_J_BP1b <- o_K1_10c$bifpoints[2,]; #unstable branch
Eq_K1_J_BP1c <- o_K1_10e$bifpoints[1,]; #inv of C1 in P-C2-R2+R1
Eq_K1_J_BP2 <- o_K1_10f$bifpoints[1,];
Eq_K1_J_BPE <- o_K1_10b$bifpoints[1,];
Eq_K1_J_LP <- o_K1_10g$bifpoints[1,];
############## K1-K, K2=3E-4 + Beta = 1.2 ========================
# Starting from P-C2-R2+R1
cmd = o_K1_11a <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, output3$curvepoints[425,2], 1E-8, output3$curvepoints[425,c(4,5,7)]),
                           0.1, c(11, 0, 1E-2), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)
cmd = o_K1_11b <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_11a$bifpoints[1,c(1:8)]),
                           0.1, c(11, 0, 1E-2), NULL, NULL, clean=FALSE)
# After extinction of C2
cmd = o_K1_11c <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_11b$bifpoints[,c(1:4,6,8)]),
                           0.3, c(11, 0, 1E-2), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = o_K1_11d <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K1_11b$bifpoints[,c(1:4,6,8)]),
                           -0.3, c(11, 0, 1E-2), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K1_K_BP1 <- o_K1_11a$bifpoints[1,];
Eq_K1_K_BP2 <- o_K1_11b$bifpoints[1,];
################## Equilibrium analyses over K2 gradient for varied K1 and Beta ===========
############## K2-A, K1=1E-5 + Beta = 0.5 ========================
cmd = o_K2_1a <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output2.2$curvepoints[6,c(3,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = o_K2_1b <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K2_1a$bifpoints[1,c(1:3,5:8)]), 0.1, c(10, 0, 2E-3), NULL, c("envZE", "2"), clean=FALSE)
cmd = o_K2_1c <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K2_1b$bifpoints[1,c(1:8)]), -0.1, c(10, 0, 2E-3), NULL, NULL, clean=FALSE)
# After extinction of C1
cmd = o_K2_1d <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K2_1c$bifpoints[1,c(1:5,7)]), -0.1, c(10, 0, 2E-3), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

Eq_K2_A_BP1<-o_K2_1d$bifpoints[1,];
Eq_K2_A_BP2<-o_K2_1a$bifpoints[1,];
Eq_K2_A_BPE<-o_K2_1b$bifpoints[1,];
Eq_K2_A_LP<-o_K2_1d$bifpoints[2,];

############## K2-B, K1=1E-5 + Beta = 2  ========================
cmd = o_K2_2a <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output2.2$curvepoints[9,c(3,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = o_K2_2b <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K2_2a$bifpoints[1,c(1:3,5:8)]), 0.1, c(10, 0, 2E-3), NULL, c("envZE", "2"), clean=FALSE)
cmd = o_K2_2c <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K2_2b$bifpoints[1,c(1:8)]), -0.1, c(10, 0, 2E-3), NULL, NULL, clean=FALSE)
# After extinction of C1
cmd = o_K2_2d <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K2_2c$bifpoints[3,c(1:5,7)]), 0.1, c(10, 0, 2E-3), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

Eq_K2_B_BP1<-o_K2_2c$bifpoints[3,];
Eq_K2_B_BP2<-o_K2_2a$bifpoints[1,];
Eq_K2_B_BPE<-o_K2_2b$bifpoints[1,];
Eq_K2_B_LP<-o_K2_2c$bifpoints[2,];

############## K2-C, K1=2E-5 + Beta = 0.27  ========================
cmd = o_K2_3a <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output3.2$curvepoints[26,c(3,4,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = o_K2_3b <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K2_3a$bifpoints[1,c(1:8)]), 0.1, c(10, 0, 9.5E-4), NULL, NULL, clean=FALSE)
cmd = o_K2_3c <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K2_3b$bifpoints[1,c(1:5,7)]), -0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

Eq_K2_C_BP1<-o_K2_3b$bifpoints[1,];
Eq_K2_C_BP2<-o_K2_3a$bifpoints[1,];
Eq_K2_C_LP<-o_K2_3c$bifpoints[1,];

############## K2-D, K1=2E-5 + Beta = 0.5 ========================
cmd = o_K2_4a <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output2.2$curvepoints[21,c(3,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = o_K2_4b <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K2_4a$bifpoints[1,c(1:3,5:8)]), 0.1, c(10, 0, 2E-3), NULL, c("envZE", "2"), clean=FALSE)
cmd = o_K2_4c <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K2_4b$bifpoints[1,c(1:8)]), -0.1, c(10, 0, 2E-3), NULL, NULL, clean=FALSE)
# After extinction of C1
cmd = o_K2_4d <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K2_4c$bifpoints[3,c(1:5,7)]), -0.1, c(10, 0, 2E-3), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

Eq_K2_D_BP1<-o_K2_4c$bifpoints[3,];
Eq_K2_D_BP2<-o_K2_4a$bifpoints[1,];
Eq_K2_D_BPE<-o_K2_4b$bifpoints[1,];
Eq_K2_D_LP<-o_K2_4c$bifpoints[2,];
Eq_K2_D_LPb<-o_K2_4d$bifpoints[1,];
############## K2-E, K1=1E-4 + Beta = 0.5 ========================
cmd = o_K2_5a <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output3.2$curvepoints[157,c(3,4,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = o_K2_5b <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K2_5a$bifpoints[1,c(1:8)]), 0.1, c(10, 0, 2E-3), NULL, NULL, clean=FALSE)
cmd = o_K2_5c <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K2_5b$bifpoints[1,c(1:5,7)]), 0.1, c(10, 0, 2E-3), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

Eq_K2_E_BP1<-o_K2_5b$bifpoints[1,];
Eq_K2_E_BP2<-o_K2_5a$bifpoints[1,];
Eq_K2_E_LP<-o_K2_5c$bifpoints[1,];
############## K2-F, K1=1E-4 + Beta = 1.0 ========================
## Starting from EQ C1-R1+R2
cmd = o_K2_6a <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output2.2$curvepoints[62,c(3,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = o_K2_6b <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K2_6a$bifpoints[1,c(1:3,5:8)]), 0.1, c(10, 0, 2E-3), NULL, c("envZE", "2"), clean=FALSE)
cmd = o_K2_6c <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K2_6b$bifpoints[1,c(1:8)]), -0.1, c(10, 0, 2E-3), NULL, NULL, clean=FALSE)
# After extinction of C2
cmd = o_K2_6d <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K2_6c$bifpoints[1,c(1:4,6,8)]), -0.1, c(10, 0, 2E-3), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

# Starting from EQ P-C1-R1+R2 
cmd = o_K2_6e <- PSPMequi("Tritrophicmod5_fin", "EQ", c(1E-8, 1E-8, output3.2$curvepoints[182,c(3,4,6,8) ]), 0.1, c(10, 0, 4E-03), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = o_K2_6f <- PSPMequi("Tritrophicmod5_fin", "EQ", c(o_K2_6e$bifpoints[1,c(1:8)]), 0.1, c(10, 0, 4E-03), NULL, NULL, clean=FALSE)
cmd = o_K2_6g <- PSPMequi("Tritrophicmod5_fin", "EQ", c(o_K2_6f$bifpoints[1,c(1:5,7)]), 0.1, c(10, 0, 4E-03), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)

Eq_K2_F_BP1<-o_K2_6f$bifpoints[1,];
Eq_K2_F_BP2<-o_K2_6a$bifpoints[1,];#inv in C1-R1+R2
Eq_K2_F_BP2b<-o_K2_6c$bifpoints[1,];#unstable branch
Eq_K2_F_BP2c<-o_K2_6e$bifpoints[1,];#inv in P-C1-R1+R2
Eq_K2_F_BPE<-o_K2_6b$bifpoints[1,];
############## K2-G, K1=1E-4 + Beta = 1.5 ========================
## Starting from EQ C1-R1+R2
cmd = o_K2_7a <- PSPMequi("Tritrophicmod5_fin", "EQ",c(1E-8, 1E-8, output2.2$curvepoints[63,c(3,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = o_K2_7b <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K2_7a$bifpoints[1,c(1:3,5:8)]), 0.1, c(10, 0, 2E-3), NULL, c("envZE", "2"), clean=FALSE)
cmd = o_K2_7c <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K2_7b$bifpoints[1,c(1:8)]), -0.1, c(10, 0, 2E-3), NULL, NULL, clean=FALSE)
# After extinction of C2
cmd = o_K2_7d <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K2_7c$bifpoints[1,c(1:4,6,8)]), -0.1, c(10, 0, 2E-3), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = o_K2_7e <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K2_7c$bifpoints[1,c(1:4,6,8)]), 0.1, c(10, 0, 2E-3), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K2_G_BP2<-o_K2_7a$bifpoints[1,];#inv in C1-R1+R2
Eq_K2_G_BP2b<-o_K2_7c$bifpoints[1,];#unstable branch
Eq_K2_G_BPE<-o_K2_7b$bifpoints[1,];

############## K2-H, K1=1E-4 + Beta = 2.0 ========================
cmd = o_K2_8a <- PSPMequi("Tritrophicmod5_fin", "EQ", c(1E-8, 1E-8, output2.2$curvepoints[63,c(3,6,8) ]), 0.1, c(10, 0, 4E-03), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = o_K2_8b <- PSPMequi("Tritrophicmod5_fin", "EQ", c(o_K2_8a$bifpoints[1,c(1:3,5:8) ]), 0.1, c(10, 0, 4E-03), NULL, c("envZE", "2"), clean=FALSE)
cmd = o_K2_8c <- PSPMequi("Tritrophicmod5_fin", "EQ", c(o_K2_8b$bifpoints[1,c(1:8) ]), -0.1, c(10, 0, 4E-03), NULL, NULL, clean=FALSE)

cmd = o_K2_8d <- PSPMequi("Tritrophicmod5_fin", "EQ", c(o_K2_8c$bifpoints[4,c(1:4,6,8) ]), 0.1, c(10, 0, 4E-03), NULL,  c("popZE", "0", "envZE", "3"), clean=FALSE)
cmd = o_K2_8e <- PSPMequi("Tritrophicmod5_fin", "EQ", c(o_K2_8c$bifpoints[4,c(1:4,6,8) ]), -0.1, c(10, 0, 4E-03), NULL,  c("popZE", "0", "envZE", "3"), clean=FALSE)

Eq_K2_H_BP2<-o_K2_8a$bifpoints[1,];
Eq_K2_H_BP2b<-o_K2_8c$bifpoints[4,];
Eq_K2_H_BPE<-o_K2_8b$bifpoints[1,];
Eq_K2_H_LP<-o_K2_8c$bifpoints[1,];
Eq_K2_H_LPb<-o_K2_8c$bifpoints[2,];
Eq_K2_H_LPc<-o_K2_8c$bifpoints[3,];

################# Additionnal Eq ====
############## K2-G, K1=8E-5 + Beta = 1.0 ========================
cmd = o_K2_7a <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output3$bifpoints[3,c(1:6,7,8)]), -0.7, c(10, 0, 9.5E-4), NULL, NULL, clean=FALSE)
#After extinction of C2
cmd = o_K2_7b <- PSPMequi("Tritrophicmod5_fin", "EQ",c(o_K2_7a$bifpoints[1,c(1:5,7)]), -0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "0", "envZE", "3"), clean=FALSE)
#After extinction of C1
cmd = o_K2_7c <- PSPMequi("Tritrophicmod5_fin", "EQ",c(output3$bifpoints[3,c(1:4,6,8)]), 0.1, c(10, 0, 9.5E-4), NULL, c("popZE", "1", "envZE", "4"), clean=FALSE)
# Bistable state in which C1 is limited
# Starting from EQ C1-R1
cmd = o_K2_7d <- PSPMequi("Tritrophicmod5_fin", "EQ", c(1E-8, 1E-8, output2.2$curvepoints[52,c(3,6,8) ]), 0.1, c(10, 0, 4E-03), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = o_K2_7e <- PSPMequi("Tritrophicmod5_fin", "EQ", c(o_K2_7d$bifpoints[,c(1:3,6:8) ]), 0.1, c(10, 0, 4E-03), NULL, c("envZE", "2", "envZE", "3"), clean=FALSE)
cmd = o_K2_7f <- PSPMequi("Tritrophicmod5_fin", "EQ", c(o_K2_7e$bifpoints[,c(1:8) ]), -0.1, c(10, 0, 4E-03), NULL, NULL, clean=FALSE)

Eq_K2_G_BP1<-output3$bifpoints[3,];
Eq_K2_G_BP2<-o_K2_7a$bifpoints[1,];
Eq_K2_G_BP2b<-o_K2_7d$bifpoints[1,];
Eq_K2_G_BP2c<-o_K2_7f$bifpoints[2,];
Eq_K2_G_BPE<-o_K2_7f$bifpoints[1,];
############## K2-I, K1=8E-5 + Beta = 2.0 ========================
cmd = o_K2_9a <- PSPMequi("Tritrophicmod5_fin", "EQ", c(1E-8, 1E-8, output2.2$curvepoints[53,c(3,6,8) ]), 0.1, c(10, 0, 4E-03), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = o_K2_9b <- PSPMequi("Tritrophicmod5_fin", "EQ", c(o_K2_9a$bifpoints[1,c(1:3,6:8) ]), 0.1, c(10, 0, 4E-03), NULL, c("envZE", "2", "envZE", "3"), clean=FALSE)
cmd = o_K2_9c <- PSPMequi("Tritrophicmod5_fin", "EQ", c(o_K2_9b$bifpoints[1,c(1:8) ]), -0.1, c(10, 0, 4E-03), NULL, NULL, clean=FALSE)

Eq_K2_I_BP2<-o_K2_9a$bifpoints[1,];
Eq_K2_I_BP2b<-o_K2_9c$bifpoints[2,];
Eq_K2_I_BPE<-o_K2_9b$bifpoints[1,];
Eq_K2_I_LP<-o_K2_9c$bifpoints[1,];
############## K2-J, K1=2E-4 + Beta = 2.3 ========================
cmd = o_K2_10a <- PSPMequi("Tritrophicmod5_fin", "EQ", c(1E-8, 1E-8, output2.2$curvepoints[77,c(3,6,8) ]), 0.1, c(10, 0, 4E-03), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = o_K2_10b <- PSPMequi("Tritrophicmod5_fin", "EQ", c(o_K2_8a$bifpoints[1,c(1:3,6:8) ]), 0.1, c(10, 0, 4E-03), NULL, c("envZE", "2", "envZE", "3"), clean=FALSE)
cmd = o_K2_10c <- PSPMequi("Tritrophicmod5_fin", "EQ", c(o_K2_8b$bifpoints[1,c(1:8) ]), -0.1, c(10, 0, 4E-03), NULL, NULL, clean=FALSE)

Eq_K2_J_BP2<-o_K2_10a$bifpoints[1,];
Eq_K2_J_BP2b<-o_K2_10c$bifpoints[5,];
Eq_K2_J_BPE<-o_K2_10b$bifpoints[1,];
Eq_K2_J_LP1<-o_K2_10c$bifpoints[2,];
Eq_K2_J_LP2<-o_K2_10c$bifpoints[3,];
Eq_K2_J_LP3<-o_K2_10c$bifpoints[3,];
############## K2-K, K1=3E-4 + Beta = 2.5 ========================
cmd = o_K2_11a <- PSPMequi("Tritrophicmod5_fin", "EQ", c(1E-8, 1E-8, output2.2$curvepoints[82,c(3,6,8) ]), 0.1, c(10, 0, 4E-03), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = o_K2_11b <- PSPMequi("Tritrophicmod5_fin", "EQ", c(o_K2_11a$bifpoints[1,c(1:3,6:8) ]), 0.1, c(10, 0, 4E-03), NULL, c("envZE", "2", "envZE", "3"), clean=FALSE)
cmd = o_K2_11c <- PSPMequi("Tritrophicmod5_fin", "EQ", c(o_K2_11b$bifpoints[1,c(1:8) ]), -0.1, c(10, 0, 4E-03), NULL, NULL, clean=FALSE)

Eq_K2_K_BP2<-o_K2_11a$bifpoints[1,];
Eq_K2_K_BP2b<-o_K2_11c$bifpoints[2,];
Eq_K2_K_BPE<-o_K2_11b$bifpoints[1,];
Eq_K2_K_LP<-o_K2_11c$bifpoints[1,];
############## K2-L, K1=4E-4 + Beta = 2.5 ========================
cmd = o_K2_12a <- PSPMequi("Tritrophicmod5_fin", "EQ", c(1E-8, 1E-8, output2.2$curvepoints[82,c(3,6,8) ]), 0.1, c(10, 0, 4E-03), NULL, c("popZE", "0", "envZE", "2", "envZE", "3"), clean=FALSE)
cmd = o_K2_12b <- PSPMequi("Tritrophicmod5_fin", "EQ", c(o_K2_12a$bifpoints[1,c(1:3,6:8) ]), 0.1, c(10, 0, 4E-03), NULL, c("envZE", "2", "envZE", "3"), clean=FALSE)
cmd = o_K2_12c <- PSPMequi("Tritrophicmod5_fin", "EQ", c(o_K2_12b$bifpoints[1,c(1:8) ]), -0.1, c(10, 0, 4E-03), NULL, NULL, clean=FALSE)

Eq_K2_L_BP2<-o_K2_12a$bifpoints[1,];
Eq_K2_L_BP2b<-o_K2_12c$bifpoints[4,];
Eq_K2_L_BPE<-o_K2_12c$bifpoints[1,];
Eq_K2_L_LP1<-o_K2_12c$bifpoints[1,]
Eq_K2_L_LP2<-o_K2_12c$bifpoints[2,]
Eq_K2_L_LP3<-o_K2_12c$bifpoints[3,]
