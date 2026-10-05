library(mgcv)

for (f in c('Cropland.csv','Woodland.csv','Grassland.csv')) { # Tables containing for each site the bacterial+fungal keystone richness, the soil pH in CaCl2 and the latitude+longitude coordinates 
    df<-read.csv(f)
    print(f)
    print('BACTERIA:')
    fit_spatial <- gam(
      Nbac_keystones ~ s(pH_CaCl2, k = 5) + s(LONGITUDE, LATITUDE, bs = "gp", k = 50),
      data = df,
      method = "REML",
      family = nb() # since keystone richness is count data
    )
    print(summary(fit_spatial))
    print('FUNGI:')
    fit_spatial <- gam(
      Nfun_keystones ~ s(pH_CaCl2, k = 5) + s(LONGITUDE, LATITUDE, bs = "gp", k = 50),
      data = df,
      method = "REML",
      family = nb() # since keystone richness is count data
    )
    print(summary(fit_spatial))
}