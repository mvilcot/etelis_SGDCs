library(vegan)
environment <- data.frame(  Ph = Ph ,Salinity =Salinity ,Temperature = Temperature,
                            MarineEcosDependency = MarineEcosDependency, HDI2019 = HDI2019, Gravity = Gravity)
dbRDA = capscale(q_matrix ~ Ph + Salinity + Temperature  + MarineEcosDependency + HDI2019 + Gravity , environment, dist=“bray”)
plot(dbRDA)