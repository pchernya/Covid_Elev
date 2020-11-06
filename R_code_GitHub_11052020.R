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
# Map of centroid elevation for continental US #
# county centroids map #
library(tmap)
tm_shape(us_states) + 
  tm_borders(lwd=3) +
tm_shape(st_as_sf(cnty_elev,crs=4269)) + 
  tm_symbols(col="elevation",title.col = 'Centroid elevation\n(meters)',
             palette='-Spectral',style='cont', midpoint=NA,
             size=0.3,alpha=0.8) +
tm_legend(position=c("right","bottom"),outside=FALSE)

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
                 as.factor(RUCC_2013)*scale(POP10/HHD10) + scale(elev_m) + 
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
                  as.factor(RUCC_2013)*scale(POP10/HHD10) + s(scale(elev_m),bs='ds') +
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
                as.factor(RUCC_2013)*scale(POP10/HHD10) + scale(elev_m) + 
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
                as.factor(RUCC_2013)*scale(POP10/HHD10) + scale(elev_m) + 
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
vis<-getViz(m_tw120s) 
p1_<-plot(sm(vis,1),trans=exp)    #elevation
p1<-p1_ + l_fitLine(colour = "red",size=2) + l_rug(alpha=0.5) +
     l_ciLine(mul = 5, colour = "blue", linetype = 2, size=1.5) +
     xlab('Centered and scaled centroid elevation (SD=495m)') + 
     ylab('Incidence Relative Risk (exp(Basis coefficient))') + 
     ylim(0,2.0) + geom_hline(yintercept=1,colour='gray',linetype=3,size=2) +
     ggtitle('120-day Incidence') + theme_bw() +
     theme(panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank())
p1
###
vis<-getViz(m_tw90s) 
p2_<-plot(sm(vis,1),trans=exp)    #elevation
p2<-p2_ + l_fitLine(colour = "red",size=2) + l_rug(alpha=0.5) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2, size=1.5) +
  xlab('Centered and scaled centroid elevation (SD=495m)') + 
  ylab('Incidence Relative Risk (exp(Basis coefficient))') + 
  ylim(0,2.0) + geom_hline(yintercept=1,colour='gray',linetype=3,size=2) +
  ggtitle('90-day Incidence') + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p2
###
vis<-getViz(m_tw30s) 
p3_<-plot(sm(vis,1),trans=exp)    #elevation
p3<-p3_ + l_fitLine(colour = "red",size=2) + l_rug(alpha=0.5) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2, size=1.5) +
  xlab('Centered and scaled centroid elevation (SD=495m)') + 
  ylab('Incidence Relative Risk (exp(Basis coefficient))') + 
  ylim(0,2.0) + geom_hline(yintercept=1,colour='gray',linetype=3,size=2) +
  ggtitle('30-day Incidence') + theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p3
### final figure
library(ggpubr)
ggarrange(p1$ggObj,p2$ggObj,p3$ggObj,nrow=1,ncol=3)
