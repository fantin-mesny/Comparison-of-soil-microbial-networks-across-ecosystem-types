library(indicspecies)

for (org in c('bacteria','fungi')){ 
    dat<-read.csv(paste('allEnvs.Genus.pa.',org,'.csv',sep=''),row.names='X') # load tables with presence absence of microbial genus, generated with file otuTableToGenusPresenceAbsence.py
    env<-read.csv('allEnvs.env.csv',row.names='group') # links each sample in dat to an ecosystem type 
    env<-env16$treatment

    phi1 <- multipatt(dat, env, duleg=TRUE, func = "r.g", control = how(nperm=999)) # compute the indicator genus analysis 
    write.csv(data.frame(phi1$sign),paste(org,'_indicGenus.csv',sep=''))
}

#save.image(file = "sessionPA_genusLevel.RData")