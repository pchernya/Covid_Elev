## R code for:
# Impact of altitude on COVID-19 infection and death in the United States: 
# a modeling and observational study" submitted to PLOS One. 
# Authors: Kenton E Stephens, Pavel Chernyavskiy, and Danielle R Bruns
##################################################################################################
library(sf)
library(sp)

# Code to compute elevation in meters of continental US county centroids
#* find county centroids, find their elevation in meters and append
cnty_cent<-st_centroid(st_geometry(cnty_cases))
prj_dd <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

library(elevatr)
library(rgdal)
cnty_elev<-get_elev_point(data.frame(st_coordinates(cnty_cent)),
                          prj=prj_dd,src="aws")
cnty_cases$elev_m<-cnty_elev$elevation
# note that SD of elevation is 495 meters
# elevation_m is included as a column inside cnty_cases


############################################################################################
# Statistical modeling #
#* All models below use cnty_cases as the analysis dataset, which is provided on GitHub *#

#* make a list of neighbors for each continental US county *#
library(spdep)
nb<-poly2nb(cnty_cases,snap=1.5*sqrt(.Machine$double.eps))
names(nb)<-as.character(cnty_cases$UID)
# Below must evaluate to TRUE; we now use UID in model formula
all.equal(sort(names(nb)), sort(levels(as.factor(cnty_cases$UID))))


#* 120 day incidence *#
m_tw120 <- gam(Confirmed_120d ~ offset(log(POP10)) + 
                 as.factor(RUCC_2013)*scale(POP10/HHD10) + 
                 scale(elev_m) + scale(median_age) +
                 s(as.factor(cnty_cases$Province_State),bs='re') +
                 s(as.factor(UID), bs = 'mrf', 
                   xt = list(nb = nb), k=1000),
               data = cnty_cases,
               method = 'REML', 
               family = tw()) #Tweedie
m_tw120$outer.info$conv # check convergence
summary(m_tw120)
###
m_tw120s <- gam(Confirmed_120d ~ offset(log(POP10)) + 
                  as.factor(RUCC_2013)*scale(POP10/HHD10) + 
                  s(scale(median_age),bs='ds') + s(scale(elev_m),bs='ds') +
                  s(as.factor(cnty_cases$Province_State),bs='re') +
                  s(as.factor(UID), bs = 'mrf', 
                    xt = list(nb = nb), k=1000),
                data = cnty_cases,
                method = 'REML', 
                family = tw()) #Tweedie
m_tw120s$outer.info$conv # check convergence
summary(m_tw120s)

###
m_nb120 <- gam(Confirmed_120d ~ offset(log(POP10)) + 
                 as.factor(RUCC_2013)*scale(POP10/HHD10) + scale(elev_m) +
                 s(as.factor(cnty_cases$Province_State),bs='re') +
                 s(as.factor(UID), bs = 'mrf', 
                   xt = list(nb = nb), k=1000),
               data = cnty_cases,
               method = 'REML', 
               family = nb()) #Negative Binomial
m_nb120$outer.info$conv # check convergence
summary(m_nb120)
###
m_nb120s <- gam(Confirmed_120d ~ offset(log(POP10)) + 
                  as.factor(RUCC_2013)*scale(POP10/HHD10) + s(scale(elev_m),bs='ds') +
                  s(as.factor(cnty_cases$Province_State),bs='re') +
                  s(as.factor(UID), bs = 'mrf', 
                    xt = list(nb = nb), k=1000),
                data = cnty_cases,
                method = 'REML', 
                family = nb()) #Negative Binomial
m_nb120s$outer.info$conv # check convergence
summary(m_nb120s)
###
AIC(m_tw120,m_tw120s,m_nb120,m_nb120s)


###############################################################################################
#* 90 day incidence *#
m_tw90 <- gam(Confirmed_90d ~ offset(log(POP10)) + 
                as.factor(RUCC_2013)*scale(POP10/HHD10) + 
                scale(elev_m) + scale(median_age) +
                s(as.factor(cnty_cases$Province_State),bs='re') +
                s(as.factor(UID), bs = 'mrf', 
                  xt = list(nb = nb), k=1000),
              data = cnty_cases,
              method = 'REML', 
              family = tw()) #Tweedie
m_tw90$outer.info$conv # check convergence
summary(m_tw90)
###
m_tw90s <- gam(Confirmed_90d ~ offset(log(POP10)) + 
                 as.factor(RUCC_2013)*scale(POP10/HHD10) + s(scale(elev_m),bs='ds') +
                 s(as.factor(cnty_cases$Province_State),bs='re') +
                 s(as.factor(UID), bs = 'mrf', 
                   xt = list(nb = nb), k=1000),
               data = cnty_cases,
               method = 'REML', 
               family = tw()) #Tweedie
m_tw90s$outer.info$conv # check convergence
summary(m_tw90s)

###
m_nb90 <- gam(Confirmed_90d ~ offset(log(POP10)) + 
                as.factor(RUCC_2013)*scale(POP10/HHD10) + scale(elev_m) +
                s(as.factor(cnty_cases$Province_State),bs='re') +
                s(as.factor(UID), bs = 'mrf', 
                  xt = list(nb = nb), k=1000),
              data = cnty_cases,
              method = 'REML', 
              family = nb()) #Negative Binomial
m_nb90$outer.info$conv # check convergence
summary(m_nb90)
###
m_nb90s <- gam(Confirmed_90d ~ offset(log(POP10)) + 
                 as.factor(RUCC_2013)*scale(POP10/HHD10) + s(scale(elev_m),bs='ds') +
                 s(as.factor(cnty_cases$Province_State),bs='re') +
                 s(as.factor(UID), bs = 'mrf', 
                   xt = list(nb = nb), k=1000),
               data = cnty_cases,
               method = 'REML', 
               family = nb()) #Negative Binomial
m_nb90s$outer.info$conv # check convergence
summary(m_nb90s)
###
AIC(m_tw90,m_tw90s,m_nb90,m_nb90s)


###############################################################################################
#* 30 day incidence *#
m_tw30 <- gam(Confirmed_30d ~ offset(log(POP10)) + 
                as.factor(RUCC_2013)*scale(POP10/HHD10) + 
                scale(elev_m) + scale(median_age) +
                s(as.factor(cnty_cases$Province_State),bs='re') +
                s(as.factor(UID), bs = 'mrf', 
                  xt = list(nb = nb), k=1000),
              data = cnty_cases,
              method = 'REML', 
              family = tw()) #Tweedie
m_tw30$outer.info$conv # check convergence
summary(m_tw30)
###
m_tw30s <- gam(Confirmed_30d ~ offset(log(POP10)) + 
                 as.factor(RUCC_2013)*scale(POP10/HHD10) + s(scale(elev_m),bs='ds') +
                 s(as.factor(cnty_cases$Province_State),bs='re') +
                 s(as.factor(UID), bs = 'mrf', 
                   xt = list(nb = nb), k=1000),
               data = cnty_cases,
               method = 'REML', 
               family = tw()) #Tweedie
m_tw30s$outer.info$conv # check convergence
summary(m_tw30s)

####
m_nb30 <- gam(Confirmed_30d ~ offset(log(POP10)) + 
                as.factor(RUCC_2013)*scale(POP10/HHD10) + scale(elev_m) +
                s(as.factor(cnty_cases$Province_State),bs='re') +
                s(as.factor(UID), bs = 'mrf', 
                  xt = list(nb = nb), k=1000),
              data = cnty_cases,
              method = 'REML', 
              family = nb()) #Negative Binomial
m_nb30$outer.info$conv # check convergence
summary(m_nb30)
###
m_nb30s <- gam(Confirmed_30d ~ offset(log(POP10)) + 
                 as.factor(RUCC_2013)*scale(POP10/HHD10) + s(scale(elev_m),bs='ds') +
                 s(as.factor(cnty_cases$Province_State),bs='re') +
                 s(as.factor(UID), bs = 'mrf', 
                   xt = list(nb = nb), k=1000),
               data = cnty_cases,
               method = 'REML', 
               family = nb()) #Negative Binomial
m_nb30s$outer.info$conv # check convergence
summary(m_nb30s)
AIC(m_tw30,m_tw30s,m_nb30,m_nb30s)


#############################################################################################
#* Visualization of non-linear elevation effects *#
library(mgcViz)
# choose a model (one model per figure panel)
vis<-getViz(m_tw120s) 
# vis<-getViz(m_tw90s) 
# vis<-getViz(m_tw30) 

p1<-plot(sm(vis, 2))    #elevation
p1+ l_fitLine(colour = "red") + l_rug(alpha=0.5) +
    l_ciLine(mul = 5, colour = "blue", linetype = 2) +
    xlab('Centered and scaled elevation (SD=495m)') + 
    ylab('Coefficient (contrib. to Incidence)') + #formally: basis function coefficient
    ylim(-1.0, 1.0) +
    ggtitle('120-day Incidence') + 
#    ggtitle('90-day Incidence') + 
#    ggtitle('30-day Incidence') + 
    theme_bw()
