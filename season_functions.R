library(FLCore)
library(ggplotFL)
library(FLBRP)
library(FLash)
#library(mydas)
#library(FLife)



seasonalise<-function(object, season=1:4){
  
  ## Stock and recruit                                   ###
  ## set expected to 1 and model variability as deviates ###
  sr=as.FLSR(object,model="segreg")
  params(sr)=FLPar(NA,dimnames=list(params=c("a","b"),iter=1))

    recs=FLCore:::expand(rec(object),season=season)
  
  ## Add seasons                                         ###
  object=FLCore:::expand(object,season=season)
  
  ## Divide up mortality by season                       ###
  m(      object)=m(object)/dim(object)[4]
  harvest(object)=harvest(object)/dim(object)[4]
  # change this for a seasonal fishery?
 # harvest(object)[,,,1,,]=harvest(object)*0.1
 # harvest(object)[,,,2,,]=harvest(object)*0
 # harvest(object)[,,,c(3,4),,]=harvest(object)*0.45
  ## Seasonal growth                                     ###
  stock.wt(   object)[,-dim(stock.wt(   object))[2]]=wtInterp(stock.wt(   object))
  catch.wt(   object)[,-dim(catch.wt(   object))[2]]=wtInterp(catch.wt(   object))
  landings.wt(object)[,-dim(landings.wt(object))[2]]=wtInterp(landings.wt(object))
  discards.wt(object)[,-dim(discards.wt(object))[2]]=wtInterp(discards.wt(object))
  
  object=adjust(object)
  
  ## Project for historic F                              ###
  #fbar=as(FLQuants("fbar"=fbar(object)[,-1]),"fwdControl")
  #object=fwd(object,control=fbar,sr=sr,residuals=recs)
  
  return(object,sr)}

wtInterp<-function(wt){
  
  mlt=wt[,-dim(wt)[2]]
  mlt=FLQuant(rep(seq(0,(dim(mlt)[4]-1)/dim(mlt)[4],1/dim(mlt)[4]),each=max(cumprod(dim(mlt)[1:3]))),
              dimnames=dimnames(mlt))[-dim(mlt)[1]]
  
  incmt=wt[,,,1]
  incmt=-incmt[-dim(incmt)[1],-dim(incmt)[2]]+incmt[-1,-1]
  
  wt[dimnames(incmt)$age,dimnames(incmt)$year]=
    wt[dimnames(incmt)$age,dimnames(incmt)$year]+
    incmt%+%(mlt%*%incmt)
  
  wt[,-dim(wt)[2]]}

annualise<-function(x) {
  
  # ADD slots (catch/landings.discards, m)
  res <- qapply(x, function(s) unitSums(seasonSums(s)))
  
  # MEAN wts
  catch.wt(res)[] <- unitSums(seasonSums(catch.wt(x) * catch.n(x))) / 
    unitSums(seasonSums(catch.n(x)))
  
  landings.wt(res)[] <- unitSums(seasonSums(landings.wt(x) * landings.n(x))) / 
    unitSums(seasonSums(landings.n(x)))
  
  discards.wt(res)[] <- unitSums(seasonSums(discards.wt(x) * discards.n(x))) / 
    unitSums(seasonSums(discards.n(x)))
  discards.wt(res)[is.na(discards.wt(res))]=seasonMeans(discards.wt(x))[is.na(discards.wt(res))]
  
  stock.wt(res)[] <- unitSums(seasonSums(stock.wt(x) * stock.n(x))) / 
    unitSums(seasonSums(stock.n(x)))
  
  
  # RECONSTRUCT N: N0 = N1 / exp(-M - F)
  
  stkn <- unitSums(stock.n(x)[,,,4])
  m(res) <- unitMeans(seasonSums(m(x)))
  harvest(res) <- unitMeans(seasonSums(harvest(x)))
  
  stock.n(res)[] <- stkn / exp(-m(res) - harvest(res))
  
  # mat
  
  mat(res)[] <- unitSums(seasonSums(mat(x) * stock.n(x))) / 
    unitSums(seasonSums(stock.n(x)))
  mat(res) <- mat(res) %/% apply(mat(res), c(1,3:6), max)
  mat(res)[is.na(mat(res))] <- 0
  
  # spwn
  
  m.spwn(res) <- 0.5
  harvest.spwn(res) <- 0.5
  
  # totals
  
  catch(res) <- computeCatch(res)
  landings(res) <- computeLandings(res)
  discards(res) <- computeDiscards(res)
  stock(res) <- computeStock(res)
  
  return(res)}


adjust<-function(object){
  dim=dim(object)
  
  un =units(catch.n( object))
  uwt=units(catch.wt(object))
  
  n  =stock.n(object)
  m  =m(object)
  f  =harvest(object)
  pg=stock.n(object)[dim[1],,,dim[4]]*exp(-f[dim[1],,,dim[4]]-m[dim[1],,,dim[4]])
  
  for (i in seq(dim(object)[2]-1))
    for (j in seq(dim(object)[4])){
      if (j!=dim(object)[4])
        stock.n(object)[,i,,j+1]=stock.n(object)[,i,,j]*exp(-f[,i,,j]-m[,i,,j])
      else{
        stock.n(object)[-1,i+1,,1]=stock.n(object)[-dim[1],i,,j]*exp(-f[-dim[1],i,,j]-m[-dim[1],i,,j])
        stock.n(object)[dim[1],i+1,,1]=stock.n(object)[dim[1],i+1,,1]+pg[,i,,1]}
    }
  
  catch.n(object)=stock.n(object)*f/(m+f)*(1-exp(-f-m))
  landings.n(object)[is.na(landings.n(object))]=0
  discards.n(object)[is.na(discards.n(object))]=0
  
  landings.n(object)=catch.n(object)*discards.n(object)/(discards.n(object)+landings.n(object))
  discards.n(object)=catch.n(object)-landings.n(object)
  
  units(catch.n(object))   =un
  units(landings.n(object))=un
  units(discards.n(object))=un
  
  units(catch.wt(object))   =uwt
  units(landings.wt(object))=uwt
  units(discards.wt(object))=uwt
  
  catch(object)=computeCatch(object,"all")  
  
  object}
