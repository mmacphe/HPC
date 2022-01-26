require(diversitree)
require(devtools)
install_version("phangorn", version="2.7.1", upgrade = "never", repos="https://cloud.r-project.org")
library(phangorn)
require(phytools)
require(geiger)		# More functions and models of evolution
require(paleoTS)    # Useful functions for Akaike Weights calculations
require(mailR) #to send emails
require(dplyr)

write.csv(c(0,1,1,2), file="packages_loaded.csv")

setwd("./")

write.csv(c(0,1,1,3), file="directory_set_to_current.csv")
getwd()

phy<-read.tree("/home/maggie/hpc/MacPherson_Tyrannidae_subspecies.tre") #my subspecies tree (ntips=998)

write.csv(c(0,1,1,4), file="phylogeny_loaded.csv")

#load migration states for the subspecies
Migration_status=read.csv("Tyrannidae_subspecies_migration_status.csv",row.names=5) #We need a Migratory_Status_Table from R file 2a with status for each taxon (n=998)
Migration_status$Status=ifelse(Migration_status$Resid_full==1, "Resident",
                                ifelse(Migration_status$Resid_partial==1, "Partial",
                                       ifelse(Migration_status$Migr_dir_full==1, "Migratory",
                                              NA)))

#Change states from words to numbers using dplyr
Migration_status$Status=recode(Migration_status$Status, "Resident"=1, "Partial"=2, "Migratory"=3)

states=Migration_status$Status
names(states)=row.names(Migration_status)

#Match migration state data to names in phylogeny 
names(states)[is.na(names(states))]=""
names(states)=trimws(names(states))#removes leading/trailing white space
names(states)=gsub(" ","_",paste(names(states)))

write.csv(c(0,1,1,5), file="Migration_states_curated.csv")

#Resolves polytomies into bifurcating ###
phy=multi2di(phy)
phy=force.ultrametric(phy,method="extend")

#build likelihood function for characters with more than two states
diversitree:::default.argnames.musse(3) #shows you the default argnames
lik=make.musse(phy, states, 3) #this is the likelihood function for 3 states
#get the parameters for lik as well

#constrain the model (here making all transitions equal)
lik.base=constrain(lik, lambda2~lambda1, lambda3~lambda1, 
                    mu2~mu1, mu3~mu1,
                    q13~q31, q21~q12, q23~q12, q32~q12)
argnames(lik.base)

#for the 3 trait character
p=starting.point.musse(phy,3)
start_time=Sys.time()
fit.base=find.mle(lik.base, p[argnames(lik.base)])
end_time=Sys.time()
time_lik.base=end_time-start_time
time_lik.base 

write.csv(c(0,1,1,2), file="base_model_fit.csv")

#email notifications
sender = "maggie.macpherson@gmail.com"
recipients = c("maggie.macpherson@gmail.com")
PASS= as.character("Tyrannu$2018")

send.mail(from = sender,
          to = recipients,
          subject = paste0("Base Model complete"),
          body = paste0(time_lik.base),
          smtp = list(host.name = "smtp.gmail.com", port = 465, #port options 25, 465, 587, 2525
                      user.name="maggie.macpherson", 
                      passwd = PASS, ssl = TRUE),
          authenticate = TRUE,
          send = TRUE)
write.csv(c(0,1,1,2), file="first_email_sent.csv")

#Remember "Resident"=1, "Partial"=2, "Migratory"=3.

## Allow transition rates to vary
#Script to fit several models of discrete-trait evolution to a tree.
#Model A
A=constrain(lik,lambda2~lambda1, lambda3~lambda1, mu2~mu1, mu3~mu1)
start_time=Sys.time()
AFit=find.mle(A,p[argnames(A)])
end_time=Sys.time()
time_ModelA=end_time-start_time
time_ModelA

send.mail(from = sender,
          to = recipients,
          subject = paste0("ModelA complete"),
          body = paste0(time_ModelA),
          smtp = list(host.name = "smtp.gmail.com", port = 465, #port options 25, 465, 587, 2525
                      user.name="maggie.macpherson", 
                      passwd = PASS, ssl = TRUE),
          authenticate = TRUE,
          send = TRUE)

#Model B
B=constrain(lik,lambda2~lambda1, lambda3~lambda1)
start_time=Sys.time()
BFit=find.mle(B,p[argnames(B)])
end_time=Sys.time()
time_ModelB=end_time-start_time
time_ModelB

send.mail(from = sender,
          to = recipients,
          subject = paste0("ModelB complete"),
          body = paste0(time_ModelB),
          smtp = list(host.name = "smtp.gmail.com", port = 465, #port options 25, 465, 587, 2525
                      user.name="maggie.macpherson", 
                      passwd = PASS, ssl = TRUE),
          authenticate = TRUE,
          send = TRUE)
#Model C
C=constrain(lik,mu2~mu1, mu3~mu1)
start_time=Sys.time()
CFit=find.mle(C,p[argnames(C)])
end_time=Sys.time()
time_ModelC=end_time-start_time
time_ModelC

send.mail(from = sender,
          to = recipients,
          subject = paste0("ModelC complete"),
          body = paste0(time_ModelC),
          smtp = list(host.name = "smtp.gmail.com", port = 465, #port options 25, 465, 587, 2525
                      user.name="maggie.macpherson", 
                      passwd = PASS, ssl = TRUE),
          authenticate = TRUE,
          send = TRUE)
#Model D
D=constrain(lik,q31~0, q13~0)
start_time=Sys.time()
DFit=find.mle(D,p[argnames(D)])
end_time=Sys.time()
time_ModelD=end_time-start_time
time_ModelD  

send.mail(from = sender,
          to = recipients,
          subject = paste0("ModelD complete"),
          body = paste0(time_ModelD),
          smtp = list(host.name = "smtp.gmail.com", port = 465, #port options 25, 465, 587, 2525
                      user.name="maggie.macpherson", 
                      passwd = PASS, ssl = TRUE),
          authenticate = TRUE,
          send = TRUE)

#Model 2
Trans2=constrain(lik,q13~0, q31~0)
start_time=Sys.time()
Trans2Fit=find.mle(Trans2,p[argnames(Trans2)])
end_time=Sys.time()
time_Model2=end_time-start_time
time_Model2

send.mail(from = sender,
          to = recipients,
          subject = paste0("Model2 complete"),
          body = paste0(time_Model2),
          smtp = list(host.name = "smtp.gmail.com", port = 465, #port options 25, 465, 587, 2525
                      user.name="maggie.macpherson", 
                      passwd = PASS, ssl = TRUE),
          authenticate = TRUE,
          send = TRUE)

#Model 3
Trans3=constrain(lik,q13~0) 
start_time=Sys.time()
Trans3Fit=find.mle(Trans3,p[argnames(Trans3)])
end_time=Sys.time()
time_Model3=end_time-start_time
time_Model3

send.mail(from = sender,
          to = recipients,
          subject = paste0("Model3 complete"),
          body = paste0(time_Model3),
          smtp = list(host.name = "smtp.gmail.com", port = 465, #port options 25, 465, 587, 2525
                      user.name="maggie.macpherson", 
                      passwd = PASS, ssl = TRUE),
          authenticate = TRUE,
          send = TRUE)

#Model 4
Trans4=constrain(lik,q31~0)
start_time=Sys.time()
Trans4Fit=find.mle(Trans4,p[argnames(Trans4)])
end_time=Sys.time()
time_Model4=end_time-start_time
time_Model4

send.mail(from = sender,
          to = recipients,
          subject = paste0("Model4 complete"),
          body = paste0(time_Model4),
          smtp = list(host.name = "smtp.gmail.com", port = 465, #port options 25, 465, 587, 2525
                      user.name="maggie.macpherson", 
                      passwd = PASS, ssl = TRUE),
          authenticate = TRUE,
          send = TRUE)

#Model 5
Trans5=constrain(lik, q23~0, q32~0)
start_time=Sys.time()
Trans5Fit=find.mle(Trans5,p[argnames(Trans5)])
end_time=Sys.time()
time_Model5=end_time-start_time
time_Model5

send.mail(from = sender,
          to = recipients,
          subject = paste0("Model5 complete"),
          body = paste0(time_Model5),
          smtp = list(host.name = "smtp.gmail.com", port = 465, #port options 25, 465, 587, 2525
                      user.name="maggie.macpherson", 
                      passwd = PASS, ssl = TRUE),
          authenticate = TRUE,
          send = TRUE)

#Model 6
Trans6=constrain(lik,q12~0)
start_time=Sys.time()
Trans6Fit=find.mle(Trans6,p[argnames(Trans6)])
end_time=Sys.time()
time_Model6=end_time-start_time
time_Model6  

send.mail(from = sender,
          to = recipients,
          subject = paste0("Model6 complete"),
          body = paste0(time_Model6),
          smtp = list(host.name = "smtp.gmail.com", port = 465, #port options 25, 465, 587, 2525
                      user.name="maggie.macpherson", 
                      passwd = PASS, ssl = TRUE),
          authenticate = TRUE,
          send = TRUE)

#Model 7
Trans7=constrain(lik,q21~0)
start_time=Sys.time()
Trans7Fit=find.mle(Trans7,p[argnames(Trans7)])
end_time=Sys.time()
time_Model7=end_time-start_time
time_Model7

send.mail(from = sender,
          to = recipients,
          subject = paste0("Model7 complete"),
          body = paste0(time_Model7),
          smtp = list(host.name = "smtp.gmail.com", port = 465, #port options 25, 465, 587, 2525
                      user.name="maggie.macpherson", 
                      passwd = PASS, ssl = TRUE),
          authenticate = TRUE,
          send = TRUE)

#Model 8
Trans8=constrain(lik,q21~0, q12~0)
start_time=Sys.time()
Trans8Fit=find.mle(Trans8,p[argnames(Trans8)])
end_time=Sys.time()
time_Model8=end_time-start_time
time_Model8

send.mail(from = sender,
          to = recipients,
          subject = paste0("Model8 complete"),
          body = paste0(time_Model8),
          smtp = list(host.name = "smtp.gmail.com", port = 465, #port options 25, 465, 587, 2525
                      user.name="maggie.macpherson", 
                      passwd = PASS, ssl = TRUE),
          authenticate = TRUE,
          send = TRUE)

#Model 9
Trans9=constrain(lik, q23~0)
start_time=Sys.time()
Trans9Fit=find.mle(Trans9,p[argnames(Trans9)])
end_time=Sys.time()
time_Model9=end_time-start_time
time_Model9 

send.mail(from = sender,
          to = recipients,
          subject = paste0("Model9 complete"),
          body = paste0(time_Model9),
          smtp = list(host.name = "smtp.gmail.com", port = 465, #port options 25, 465, 587, 2525
                      user.name="maggie.macpherson", 
                      passwd = PASS, ssl = TRUE),
          authenticate = TRUE,
          send = TRUE)

#Model 10
Trans10=constrain(lik, q32~0)
start_time=Sys.time()
Trans10Fit=find.mle(Trans10,p[argnames(Trans10)])
end_time=Sys.time()
time_Model10=end_time-start_time
time_Model10
  
send.mail(from = sender,
          to = recipients,
          subject = paste0("Model10 complete"),
          body = paste0(time_Model10),
          smtp = list(host.name = "smtp.gmail.com", port = 465, #port options 25, 465, 587, 2525
                      user.name="maggie.macpherson", 
                      passwd = PASS, ssl = TRUE),
          authenticate = TRUE,
          send = TRUE)

rates=fit.base$par.full
write.csv(rates,file="BaseModelParameters.csv")

LnLfit.base=fit.base$lnLik 
LnLA=AFit$lnLik
LnLB=BFit$lnLik
LnLC=CFit$lnLik
LnLD=DFit$lnLik 
LnLTrans2=Trans2Fit$lnLik 
LnLTrans3=Trans3Fit$lnLik
LnLTrans4=Trans4Fit$lnLik 
LnLTrans5=Trans5Fit$lnLik
LnLTrans6=Trans6Fit$lnLik 
LnLTrans7=Trans7Fit$lnLik 
LnLTrans8=Trans8Fit$lnLik 
LnLTrans9=Trans9Fit$lnLik
LnLTrans10=Trans10Fit$lnLik

results=as.data.frame(c(LnLfit.base, LnLA, LnLB, LnLC, LnLD, LnLTrans2, LnLTrans3,
                         LnLTrans4, LnLTrans5, LnLTrans6, LnLTrans7, LnLTrans8, LnLTrans9,
                         LnLTrans10))  
    
AICcfit.base=IC(fit.base$lnLik,12,method="AICc",n=phy$Nnode)
AICcA=IC(AFit$lnLik,6,method="AICc",n=phy$Nnode)
AICcB=IC(BFit$lnLik,8,method="AICc",n=phy$Nnode)
AICcC=IC(CFit$lnLik,8,method="AICc",n=phy$Nnode)
AICcD=IC(DFit$lnLik,10,method="AICc",n=phy$Nnode)
AICcTrans2=IC(Trans2Fit$lnLik,10,method="AICc",n=phy$Nnode) 
AICcTrans3=IC(Trans3Fit$lnLik,11,method="AICc",n=phy$Nnode)
AICcTrans4=IC(Trans4Fit$lnLik,11,method="AICc",n=phy$Nnode) 
AICcTrans5=IC(Trans5Fit$lnLik,10,method="AICc",n=phy$Nnode)
AICcTrans6=IC(Trans6Fit$lnLik,11,method="AICc",n=phy$Nnode) 
AICcTrans7=IC(Trans7Fit$lnLik,11,method="AICc",n=phy$Nnode)
AICcTrans8=IC(Trans8Fit$lnLik,10,method="AICc",n=phy$Nnode) 
AICcTrans9=IC(Trans9Fit$lnLik,11,method="AICc",n=phy$Nnode)
AICcTrans10=IC(Trans10Fit$lnLik,11,method="AICc",n=phy$Nnode)

results$AIC=c(AICcfit.base, AICcA, AICcB, AICcC, AICcD, AICcTrans2, AICcTrans3, 
               AICcTrans4, AICcTrans5, AICcTrans6, AICcTrans7, AICcTrans8, 
               AICcTrans9, AICcTrans10)
      
BICfit.base=IC(fit.base$lnLik,12,method="BIC",n=phy$Nnode)
BICA=IC(AFit$lnLik,6,method="BIC",n=phy$Nnode)
BICB=IC(BFit$lnLik,8,method="BIC",n=phy$Nnode)
BICC=IC(CFit$lnLik,8,method="BIC",n=phy$Nnode)
BICD=IC(DFit$lnLik,10,method="BIC",n=phy$Nnode)
BICTrans2=IC(Trans2Fit$lnLik,10,method="BIC",n=phy$Nnode) 
BICTrans3=IC(Trans3Fit$lnLik,11,method="BIC",n=phy$Nnode)
BICTrans4=IC(Trans4Fit$lnLik,11,method="BIC",n=phy$Nnode) 
BICTrans5=IC(Trans5Fit$lnLik,10,method="BIC",n=phy$Nnode) 
BICTrans6=IC(Trans6Fit$lnLik,11,method="BIC",n=phy$Nnode) 
BICTrans7=IC(Trans7Fit$lnLik,11,method="BIC",n=phy$Nnode)
BICTrans8=IC(Trans8Fit$lnLik,10,method="BIC",n=phy$Nnode) 
BICTrans9=IC(Trans9Fit$lnLik,11,method="BIC",n=phy$Nnode)
BICTrans10=IC(Trans10Fit$lnLik,11,method="BIC",n=phy$Nnode)
  
results$BIC=c(BICfit.base, BICA, BICB, BICC, BICD, BICTrans2, BICTrans3, 
               BIDTrans4, BICTrans5, BICTrans6, BICTrans7, BICTrans8,
               BICTrans9, BICTrans10)
results$Model=c("base", "A", "B", "C", "D", "Trans2", "Trans3", "Trans4", 
                 "Trans5", "Trans6", "Trans7", "Trans8", "Trans9", "Trans10")

results=results[,c(4,1,2,3)]
results=rename(results, LnLik=`c(LnLfit.base, LnLA, LnLB, LnLC, LnLTrans3, LnLTrans5, LnLTrans7, LnLTrans9, LnLTrans10)`)
write.csv(results,file="Tyrannidae_MuSSEresults.csv",row.names=FALSE)

write.csv(AFit$par.full, file="Tyrannidae_ModelAParameters.csv", row.names=TRUE)
write.csv(BFit$par.full, file="Tyrannidae_ModelBParameters.csv", row.names=TRUE)
write.csv(CFit$par.full, file="Tyrannidae_ModelCParameters.csv", row.names=TRUE)
write.csv(DFit$par.full, file="Tyrannidae_ModelDParameters.csv", row.names=TRUE)
write.csv(Trans2Fit$par.full, file="Tyrannidae_Model2Parameters.csv", row.names=TRUE)
write.csv(Trans3Fit$par.full, file="Tyrannidae_Model3Parameters.csv", row.names=TRUE)
write.csv(Trans4Fit$par.full, file="Tyrannidae_Model4Parameters.csv", row.names=TRUE)
write.csv(Trans5Fit$par.full, file="Tyrannidae_Model5Parameters.csv", row.names=TRUE)
write.csv(Trans6Fit$par.full, file="Tyrannidae_Model6Parameters.csv", row.names=TRUE)
write.csv(Trans7Fit$par.full, file="Tyrannidae_Model7Parameters.csv", row.names=TRUE)
write.csv(Trans8Fit$par.full, file="Tyrannidae_Model8Parameters.csv", row.names=TRUE)
write.csv(Trans9Fit$par.full, file="Tyrannidae_Model9Parameters.csv", row.names=TRUE)
write.csv(Trans10Fit$par.full, file="Tyrannidae_Model10Parameters.csv", row.names=TRUE)

send.mail(from = sender,
          to = recipients,
          subject = paste0("Model4 complete"),
          body = paste0(time_Model4),
          smtp = list(host.name = "smtp.gmail.com", port = 465, #port options 25, 465, 587, 2525
                      user.name="maggie.macpherson", 
                      passwd = PASS, ssl = TRUE),
          authenticate = TRUE,
          send = TRUE,
          attach.files=c("Tyrannidae_MuSSEresults.csv","Tyrannidae_ModelAParameters.csv",
                         "Tyrannidae_ModelBParameters.csv","Tyrannidae_ModelCParameters.csv",
                         "Tyrannidae_ModelDParameters.csv","Tyrannidae_Model2Parameters.csv",
                         "Tyrannidae_Model3Parameters.csv","Tyrannidae_Model4Parameters.csv",
                         "Tyrannidae_Model5Parameters.csv","Tyrannidae_Model6Parameters.csv",
                         "Tyrannidae_Model7Parameters.csv","Tyrannidae_Model8Parameters.csv",
                         "Tyrannidae_Model9Parameters.csv","Tyrannidae_Model10Parameters.csv"))