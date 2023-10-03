###################### Seasonal operating model for sprat#############
#####################Adapted from annual model (N.Walker)###########
#################### Author: Rosana Ourens

required_pckgs <- c("FLasher", "FLAssess", "ggplotFL", "FLBRP", "data.table")
source("C:/Users/ro03/Documents/sprat/MSE/season_MSE/season_functions.R")

load("C:/Users/ro03/Documents/sprat/MSE/IBP_OMs/VB4_LW1_s0.65_sR0.78_FH1.RData")


sprat4=expand(stk,season=1:4)
#Divide natural mortality equally across seasons
m( sprat4)=m(sprat4)/dim(sprat4)[4]

# Assume that the fishery is seasonal,and most of the activity takes place in quarter 3 and 4 


harvest(sprat4)[,,,1,,]<- harvest(stk)*0.06
harvest(sprat4)[,,,2,,]<- harvest(stk)*0.04
harvest(sprat4)[,,,c(3,4),,]<- harvest(stk)*0.45

#To ensure that the initial stock vector is consistent with the new structure, update the stock numbers values
#in the 1st year
stock.n(sprat4)[,1,,2]=stock.n(sprat4)[,1,,1]*exp(-m(sprat4)[,1,,1]-harvest(sprat4)[,1,,1])
stock.n(sprat4)[,1,,3]=stock.n(sprat4)[,1,,2]*exp(-m(sprat4)[,1,,2]-harvest(sprat4)[,1,,2])
stock.n(sprat4)[,1,,4]=stock.n(sprat4)[,1,,3]*exp(-m(sprat4)[,1,,3]-harvest(sprat4)[,1,,3])


## Seasonal growth                                     ###
stock.wt( sprat4)[,-dim(stock.wt(sprat4))[2]]=wtInterp(stock.wt( sprat4))
catch.wt(  sprat4)[,-dim(catch.wt( sprat4))[2]]=wtInterp(catch.wt( sprat4))
landings.wt(sprat4)[,-dim(landings.wt(sprat4))[2]]=wtInterp(landings.wt(sprat4))
discards.wt(sprat4)[,-dim(discards.wt(sprat4))[2]]=wtInterp(discards.wt(sprat4))

###
m.spwn(sprat4)[,,,c(1,2,3)]<-1
m.spwn(sprat4)[,,,4]<-0.5
harvest.spwn(sprat4)[,,,c(1,2,3)]<-1
harvest.spwn(sprat4)[,,,4]<-0.5
#mat(sprat4)[,,,-3]=NA

#S-R takes places in q4
sr.par=FLPar(NA, dimnames=list(params=c("a","b"), season=1:4, iter=1))
sr.par[1,4]=445
sr.par[2,4]=308

sprat_sr <- predictModel(params=sr.par, model = "segreg")

########### modify the residuals of the SR model
############### no adding interaction at this stage for a quick run
#sr.res =FLQuant(dimnames=list(quant="all",unit="unique",year=2:26,season=1:4,area="unique",iter=1:500))
#sr.res@.Data[,,,4,,]<-sr.om.res@.Data[,c(2:26),,,,]
sr.res =FLQuant(dimnames=list(quant="all",unit="unique",year=2:26,season=1:4,area="unique",iter=1))
sr.res@.Data[,,,c(1:2,3),,]<-NA
sr.res@.Data[,,,4,,]<-sr.om.res@.Data[,c(2:26),,,,1]
  

## forecast for the first 25 years

control=as(FLQuants(fbar=fbar(sprat4)[,c(2:25),,,,1]), 'fwdControl')

#sprat4=fwd(sprat4[,c(1:25),,,,], control=control, sr=sprat_sr, deviances=sr.res)[,-1]
sprat4=fwd(sprat4[,c(1:25),,,,1], control=control, sr=sprat_sr, deviances=sr.res)[,-1]


## compare with the annual model
plot(FLStocks("Annual"=stk,"Seasonal"=sprat4))+
  scale_colour_manual("Model",values=c("blue","purple"))+
  theme_bw()+
  theme(legend.position="bottom")

# same pattern in catch
ggplot(as.data.frame(FLQuants("annual"=catch(stk)[,c(1:25),,,,1],"Seasonal"=apply(catch(sprat4)[,,,,,1],2,sum,na.rm=T))))+
  geom_line(aes(year,data,col=qname))+
  theme_bw()+
  xlab("Year")+ylab("Catch")
#same fbar
ggplot(as.data.frame(FLQuants("annual"=fbar(stk)[,c(1:25),,,,1],"Seasonal"=apply(fbar(sprat4)[,,,,,1],2,sum, na.rm=T))))+
  geom_line(aes(year,data,col=qname))+
  theme_bw()+ylab("fbar")

#same pattern in ssb
ggplot(as.data.frame(FLQuants("annual"=ssb(stk)[,c(1:25),,,,1],"Seasonal"=apply(ssb(sprat4)[,,,,,1],2,mean, na.rm=T))))+
  geom_line(aes(year,data,col=qname))+
  theme_bw()+ylab("SSB")

#rec in annual model should be = rec in q4
ggplot(as.data.frame(FLQuants("annual"=rec(stk)[,c(1:25),,,,1],"Seasonal"=rec(sprat4)[,,,4,,1])))+
  geom_line(aes(year,data,col=qname))+
  theme_bw()+ylab("rec")
# mean seasonal rec should be lower than in annual model
ggplot(as.data.frame(FLQuants("annual"=rec(stk)[,c(1:25),,,,1],"Seasonal"=apply(rec(sprat4)[,,,,,1],2,mean, na.rm=T))))+
  geom_line(aes(year,data,col=qname))+
  theme_bw()+ylab("rec")


### ------------------------------------------------------------------------ ###
### CHECK OM IS BEHAVING
### ------------------------------------------------------------------------ ###


discards.n(sprat4)=0

### Attempt to crash the stock


ctrl <- fwdControl(data.frame(season=c(rep(seq(1:4),25)),
                              year = c(rep(26:50,each=4)),
                              quant = "f",
                              value = c(rep(1,60),rep(0,40))))


spr_stf <- stf(sprat4, nyears = 25)
pstk <- fwd(spr_stf, control = ctrl, sr = sprat_sr)

plot(pstk)




##### F=0
ctrl <- fwdControl(data.frame(season=c(rep(seq(1:4),25)),
                              year = c(rep(26:50,each=4)),
                              quant = "catch",
                              value = rep(0,100)))

pstk <- fwd(spr_stf, control = ctrl, sr = sprat_sr)

plot(pstk)




## recruitment check


rec(sprat4)[,ac(25),,4,,1] 
sprat_sr@params[1,4] * sprat_sr@params[2,4] * sr.res[,ac(25),,4,,1]

