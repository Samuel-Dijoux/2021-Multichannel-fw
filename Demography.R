###########################################################################################
#               -> Demographic analysis of size-structured consumers
#                 Over K1 productivity gradient and Beta
###########################################################################################
library(PSPManalysis)
if (!exists("PSPMsrc.fullpath")) source("../PSPMworkshop/PSPManalysis/PSPManalysis.r")
###################################################### Consumer Demography ~ Beta
################## Birth rate ====
### Computation of R-C-P EQ using Tritrophicmod5_demography.h ====
cmd = o_1 <- PSPMequi("Tritrophicmod5_demography", "EQ", c(1.0E-08, 1.0E-08),
                      0.1, c(9, 0, 4e-4), NULL, c("popZE", "0", "envZE", "1", "envZE", "2"), clean=FALSE)
cmd = o_2 <- PSPMequi("Tritrophicmod5_demography", "EQ", c(o_1$bifpoint[,c(1, 2, 5)]),
                      0.1, c(9, 0, 4e-4), NULL,c("envZE", "1", "envZE", "2"), clean=FALSE)
cmd = o_3 <- PSPMequi("Tritrophicmod5_demography", "EQ", c(o_2$bifpoint[1,c(1:5)]),
                      -0.3, c(9, 0, 4e-3), NULL, NULL, clean=FALSE)
cmd = o_3b <- PSPMequi("Tritrophicmod5_demography", "EQ", c(o_2$bifpoint[1,c(1:5)]),
                       0.3, c(9, 0, 4e-3), NULL, NULL, clean=FALSE)
### Selection for K=3E-5 ====
#without Predator
cmd = o_BR_R1a <- PSPMequi("Tritrophicmod5_demography", "EQ", c(0.5, o_2$curvepoints[55,c(2, 5)]), 0.1, c(16, 0, 4), NULL, c("envZE", "1", "envZE", "2"), clean=FALSE)
cmd = o_BR_R1b <- PSPMequi("Tritrophicmod5_demography", "EQ", c(0.5, o_2$curvepoints[55,c(2, 5)]), -0.1, c(16, 0, 4), NULL, c("envZE", "1", "envZE", "2"), clean=FALSE)
o_BR_R1<-rbind(o_BR_R1b$curvepoints[dim(o_BR_R1b$curvepoints)[1]:1,], o_BR_R1a$curvepoints)

#with Predator
cmd = PBR_R1a <- PSPMequi("Tritrophicmod5_demography", "EQ", c(0.5, o_3$curvepoints[58,c(2:5)]), 0.1, c(16, 0, 4), NULL, NULL, clean=FALSE)
cmd = PBR_R1b <- PSPMequi("Tritrophicmod5_demography", "EQ", c(0.5, o_3$curvepoints[58,c(2:5)]), -0.1, c(16, 0, 4), NULL, NULL, clean=FALSE)
P_BR_R1<-rbind(PBR_R1b$curvepoints[dim(PBR_R1b$curvepoints)[1]:1,], PBR_R1a$curvepoints)

### Selection for K=8E-5 ====
#without Predator
cmd = o_BR_R2a <- PSPMequi("Tritrophicmod5_demography", "EQ", c(0.5, o_2$curvepoints[117,c(2, 5)]), 0.1, c(16, 0, 4), NULL, c("envZE", "1", "envZE", "2"), clean=FALSE)
cmd = o_BR_R2b <- PSPMequi("Tritrophicmod5_demography", "EQ", c(0.5, o_2$curvepoints[117,c(2, 5)]), -0.1, c(16, 0, 4), NULL, c("envZE", "1", "envZE", "2"), clean=FALSE)
o_BR_R2<-rbind(o_BR_R2b$curvepoints[dim(o_BR_R2b$curvepoints)[1]:1,], o_BR_R2a$curvepoints)

#with Predator
cmd = PBR_R2a <- PSPMequi("Tritrophicmod5_demography", "EQ", c(0.5, o_3$curvepoints[95,c(2:5)]), 0.1, c(16, 0, 4), NULL, NULL, clean=FALSE)
cmd = PBR_R2b <- PSPMequi("Tritrophicmod5_demography", "EQ", c(0.5, o_3$curvepoints[95,c(2:5)]), -0.1, c(16, 0, 4), NULL, NULL, clean=FALSE)
P_BR_R2<-rbind(PBR_R2b$curvepoints[dim(PBR_R2b$curvepoints)[1]:1,], PBR_R2a$curvepoints)

### Selection for K=3E-4 ====
#without Predator
cmd = o_BR_R3a <- PSPMequi("Tritrophicmod5_demography", "EQ", c(0.5, o_2$curvepoints[222,c(2,5)]), 0.1, c(16, 0, 5), NULL,c("envZE", "1", "envZE", "2"), clean=FALSE)
cmd = o_BR_R3b <- PSPMequi("Tritrophicmod5_demography", "EQ", c(0.5, o_2$curvepoints[222,c(2,5)]), -0.1, c(16, 0, 5), NULL,c("envZE", "1", "envZE", "2"), clean=FALSE)
o_BR_R3<-rbind(o_BR_R3b$curvepoints[dim(o_BR_R3b$curvepoints)[1]:1,], o_BR_R3a$curvepoints)

#with Predator
cmd = PBR_R3a <- PSPMequi("Tritrophicmod5_demography", "EQ", c(0.5, o_3$curvepoints[136,c(2:5)]), 0.1, c(16, 0, 4), NULL, NULL, clean=FALSE)
cmd = PBR_R3b <- PSPMequi("Tritrophicmod5_demography", "EQ", c(0.5, o_3$curvepoints[136,c(2:5)]), -0.1, c(16, 0, 4), NULL, NULL, clean=FALSE)
P_BR_R3<-rbind(PBR_R3b$curvepoints[dim(PBR_R3b$curvepoints)[1]:1,], PBR_R3a$curvepoints)

################## Population growth rate ====
####No predation pressure, Computated using Tritrophicmod5_demography2.h 
### Selection for R=K=3E-5
cmd = O_PGR_R1a <- PSPMdemo("Tritrophicmod5_demography2", c(16, 0.5, 0.01, 0.0, 4), clean=F);
cmd = O_PGR_R1b <- PSPMdemo("Tritrophicmod5_demography2", c(16, 0.5, -0.01, 0.0, 4), clean=F);
O_PGR_R1<-rbind(O_PGR_R1b$curvepoints[dim(O_PGR_R1b$curvepoints)[1]:1,], O_PGR_R1a$curvepoints[2:dim(O_PGR_R1a$curvepoints)[1],]);
### Selection for R=K=8E-5
cmd = O_PGR_R2a <- PSPMdemo("Tritrophicmod5_demography2",c(16, 0.5, 0.01, 0.0, 4), clean=TRUE);
cmd = O_PGR_R2b <- PSPMdemo("Tritrophicmod5_demography2",c(16, 0.5, -0.01, 0.0, 4), clean=TRUE);
O_PGR_R2<-rbind(O_PGR_R2b$curvepoints[dim(O_PGR_R2b$curvepoints)[1]:1,], O_PGR_R2a$curvepoints[2:dim(O_PGR_R2a$curvepoints)[1],]);
### Selection for R=K=3E-4
cmd = O_PGR_R3a <- PSPMdemo("Tritrophicmod5_demography2",c(16, 0.5, 0.01, 0.0, 5), clean=TRUE);
cmd = O_PGR_R3b <- PSPMdemo("Tritrophicmod5_demography2",c(16, 0.5, -0.01, 0.0, 5), clean=TRUE);
O_PGR_R3<-rbind(O_PGR_R3b$curvepoints[dim(O_PGR_R3b$curvepoints)[1]:1,], O_PGR_R3a$curvepoints[2:dim(O_PGR_R3a$curvepoints)[1],]);

################## Individual development time ====
Beta<-seq(0,4,0.001);     # Relative body size of consumers
R<-3E-4;                  # Fixed resource productivity
### parameters settings
RH <-1.5E-5;  
NU <-0.006;   
RM <-0.003;
lb <-7.0*Beta;                #Size at birth  
lj <-110*Beta;                #Size at maturation  
lv <-rep(27.0,length(Beta));  #maximum size of predation predation range
lm <-300*Beta;                #Asymptotic size  
### Initializing
t_j<-NULL;                    #Time required to maturate
t_e<-NULL;                    #Time required to escape predation      
t_max<-NULL;                  #Time required to reach the maximal length      
t_j[1]<-NA; t_e[1]<-NA; t_max[1]<-NA;

#Loop
for(i in 2:length(Beta)){    # Iteration t=0, birth; for Beta > 0
  Length<-NULL;              # Storage vector of individual length
  Length[1]<-lb[i];          # State at birth
  l1 <- Length[1];           # Iteration t=1,     
  k<-2;                      # Step t+1
  dev <- NU*(lm[i]*R/(R+RH) - Length[k-1]);  #Calculus of individual growth for the first iteration
  l2 <- Length[k-1]+dev;     # Calculus of individual length at step k (t+1)
  
  if(l2>lm[i]){
    t_e[i]<-NA; t_j[i]<-NA;
  }else{
    while(dev>0.0099){
      if(l2>lv[i] && l1<lv[i]){t_e[i]<- k-1};
      if(l2>lj[i] && l1<lj[i]){t_j[i]<- k-1};
      l1<-l2;
      Length<-c(Length,l1);
      dev<-NU*(lm[i]*R/(R+RH) - l1);  #Calculus of individual growth for each iteration -t1
      l2<-l1+dev;
      k<-k+1;
    }
  }
  t_max[i]<-k-1;
  i<-i+1;
}
#Data merging
TableDR<-cbind(Beta,t_e,t_j,t_max);
TableDR<-as.data.frame(TableDR,colnames=TRUE)
TableDR<-cbind(Beta,t_e,t_j,t_max);
TableDR<-as.data.frame(TableDR,colnames=TRUE)
################## Critical resource density ====
cmd = output1.2 <- PSPMequi("Tritrophicmod5_demography3", "EQ", c(1.0E-08, 1.0E-08, 1.0E-08), 0.5, c(11, 0, 9e-2), NULL,c("popZE", "0", "popZE", "1", "envZE", "2", "envZE", "3", "envZE", "4"), clean=FALSE)
### Mup = 0; Absence of predator ====
cmd = o4 <- PSPMequi("Tritrophicmod5_demography3", "BP", c(output1.2$bifpoints[1:3], 1.2), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4), NULL,
                     c("popZE", "0", "envZE", "2", "envZE", "3", "envZE", "4", "popBP", "1"), clean=FALSE)
cmd = o4b <- PSPMequi("Tritrophicmod5_demography3", "BP", c(output1.2$bifpoints[1:3], 1.2), 0.1, c(11, 0, 9.5E-01, 20, 0.1, 4), NULL,
                      c("popZE", "0", "envZE", "2", "envZE", "3", "envZE", "4", "popBP", "1"), clean=FALSE)
o_4<-rbind(o4$curvepoints[dim(o4$curvepoints)[1]:1,], o4b$curvepoints[2:dim(o4b$curvepoints)[1],]);
### Mup = 0.1; mimicking additional predator-induced mortality rate ====
cmd = o5 <- PSPMequi("Tritrophicmod5_demography3", "BP", c(output1.2$bifpoints[1:3], 1.2), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4), NULL,
                     c("popZE", "0", "envZE", "2", "envZE", "3", "envZE", "4", "popBP", "1"), clean=FALSE)
cmd = o5b <- PSPMequi("Tritrophicmod5_demography3", "BP", c(output1.2$bifpoints[1:3], 1.2), 0.1, c(11, 0, 9.5E-01, 20, 0.1, 4), NULL,
                      c("popZE", "0", "envZE", "2", "envZE", "3", "envZE", "4", "popBP", "1"), clean=FALSE)
o_5<-rbind(o5$curvepoints[dim(o5$curvepoints)[1]:1,], o5b$curvepoints[2:dim(o5b$curvepoints)[1],]);
### Mup = 0.5; mimicking additional predator-induced mortality rate ====
cmd = o6 <- PSPMequi("Tritrophicmod5_demography3", "BP", c(output1.2$bifpoints[1:3], 1.2), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4), NULL,
                     c("popZE", "0", "envZE", "2", "envZE", "3", "envZE", "4", "popBP", "1"), clean=FALSE)
cmd = o6b <- PSPMequi("Tritrophicmod5_demography3", "BP", c(output1.2$bifpoints[1:3], 1.2), 0.1, c(11, 0, 9.5E-01, 20, 0.1, 4), NULL,
                      c("popZE", "0", "envZE", "2", "envZE", "3", "envZE", "4", "popBP", "1"), clean=FALSE)
o_6<-rbind(o6$curvepoints[dim(o6$curvepoints)[1]:1,], o6b$curvepoints[2:dim(o6b$curvepoints)[1],]);
### Mup = 0.9; mimicking additional predator-induced mortality rate ====
cmd = o7 <- PSPMequi("Tritrophicmod5_demography3", "BP", c(output1.2$bifpoints[1:3], 1.2), -0.1, c(11, 0, 9.5E-01, 20, 0.1, 4), NULL,
                     c("popZE", "0", "envZE", "2", "envZE", "3", "envZE", "4", "popBP", "1"), clean=FALSE)
cmd = o7b <- PSPMequi("Tritrophicmod5_demography3", "BP", c(output1.2$bifpoints[1:3], 1.2), 0.1, c(11, 0, 9.5E-01, 20, 0.1, 4), NULL,
                      c("popZE", "0", "envZE", "2", "envZE", "3", "envZE", "4", "popBP", "1"), clean=FALSE)
o_7<-rbind(o7$curvepoints[dim(o7$curvepoints)[1]:1,], o7b$curvepoints[2:dim(o7b$curvepoints)[1],])