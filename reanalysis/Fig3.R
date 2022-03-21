# Our reanalysis (Fig. 3b) only included patients who were <60 years old, reached CR1, and hence were considered eligible for allogeneic transplant (n = 995)
# Data was reprocessed from the start 
# By commenting out code lines 13-15 and rerunning this script, Fig. 3a can then be generated, which is adapted from Figure 5 (a) in Gerstung et al.'s paper

# 1.3.1 Libraries 
library(CoxHD)
library(mg14)
set1 <- brewer.pal(9, "Set1")

# 1.3.2.1 Clinical data
load("./data/AMLSG_Clinical_Anon.RData")

s <- clinicalData$AOD < 60 & !is.na(clinicalData$CR_date) &! clinicalData$TPL_Phase %in% c("RD1","PR1")
clinicalData<-clinicalData[s,]
clinicalData$PDID<-droplevels(clinicalData$PDID)

# 1.3.2.2 Mutation data
mutationData = read.table("./data/AMLSG_Genetic.txt", sep="\t", header=TRUE, strip.white = TRUE)
mutationData$SAMPLE_NAME <- factor(as.character(mutationData$SAMPLE_NAME), levels = levels(clinicalData$PDID)) ##Refactor
mutationTable <- (table(mutationData[mutationData$Result %in% c("ONCOGENIC","POSSIBLE") & mutationData$FINAL_CALL== "OK" ,c("SAMPLE_NAME","GENE")]) > 0)+0
dim(mutationTable)

all(rownames(mutationTable)==clinicalData$PDID)

# 1.3.3 Survival data
os <- Surv(clinicalData$OS, clinicalData$Status) #OS
t <- clinicalData$Time_Diag_TPL
t[is.na(t) | !clinicalData$TPL_Phase %in% "CR1" | !clinicalData$TPL_type %in% c("ALLO","FREMD") ] <- Inf ## Only allografts in CR1
o <- clinicalData$OS
tplIndexOs <- t < o
osTD <- Surv(time = rep(0, nrow(clinicalData)), time2=pmin(o, t), event=ifelse(tplIndexOs, 0, clinicalData$Status) )
osTD <- rbind(osTD,
              Surv(time=t[which(tplIndexOs)],
                   time2=o[which(tplIndexOs)],
                   event=clinicalData$Status[which(tplIndexOs)])
)
osTD = Surv(osTD[,1],osTD[,2],osTD[,3])
rm(o,t)
tplSplitOs <- c(1:nrow(clinicalData), which(tplIndexOs))
osYr <- os
osYr[,1] <- osYr[,1]/365
osYrTD <- osTD
osYrTD[,1] <- osYrTD[,1]/365

# 1.3.4 Covariates
dataList <-list(Genetics = data.frame(mutationTable[,colSums(mutationTable)>0]),
                Cytogenetics = clinicalData[,grep("^(t_)|(inv)|(abn)|(plus)|(minus)|(mono)|(complex)",colnames(clinicalData))],
                Nuisance = data.frame(MakeInteger(clinicalData$Study)[,1:2], Date=scale(as.numeric(clinicalData$ERDate),scale=FALSE), MissingCyto=is.na(clinicalData$t_15_17)+0),
                Treatment = data.frame(ATRA = clinicalData$ATRA_arm, VPA=clinicalData$VPA, TPL_os=tplIndexOs),
                Demographics = clinicalData[,c("AOD","gender")],
                Clinical = cbind(clinicalData[, c("Performance_ECOG","BM_Blasts","PB_Blasts","wbc","LDH","HB","platelet",
                                                  "Splenomegaly")], MakeInteger(clinicalData$TypeAML)[,-1]))#,
#MolRisk = makeInteger(clinicalData$M_Risk))
#dataList$Genetics$CEBPA <- clinicalData$CEBPA # encoded as 0,1,2
dataList$Genetics$CEBPA_mono <- clinicalData$CEBPA == 1 # encoded as 0,1,2
dataList$Genetics$CEBPA_bi <- clinicalData$CEBPA == 2 # encoded as 0,1,2
dataList$Genetics$CEBPA <- NULL
dataList$Genetics$FLT3 <- NULL
dataList$Genetics$FLT3_ITD <- clinicalData$FLT3_ITD != "0"
dataList$Genetics$FLT3_TKD <- clinicalData$FLT3_TKD != "0"
dataList$Genetics$FLT3_other <- clinicalData$FLT3_other != "0"
dataList$Genetics$IDH2_p172 <- table(mutationData$SAMPLE_NAME[mutationData$GENE=='IDH2' & grepl("172", mutationData$AA_CHANGE)])[]
dataList$Genetics$IDH2_p140 <- table(mutationData$SAMPLE_NAME[mutationData$GENE=='IDH2' & grepl("140", mutationData$AA_CHANGE)])[]
dataList$Genetics$IDH2 <- NULL
dataList$Genetics$NPM1 <- clinicalData$NPM1
dataList$Cytogenetics$MLL_PTD <- NULL
#dataList$Genetics = dataList$Genetics + 0
dataList$Genetics[,52:58]<-as.data.frame(lapply(dataList$Genetics[,52:58],as.numeric))
dataList$GeneGene <- MakeInteractions(data.frame(dataList$Genetics), data.frame(dataList$Genetics))[,as.vector(upper.tri(matrix(0,ncol=ncol(dataList$Genetics), nrow=ncol(dataList$Genetics))))]
dataList$GeneGene <- dataList$GeneGene[,colSums(dataList$GeneGene, na.rm=TRUE)>0]
dataList$GeneGene$`NPM1:FLT3_ITD:DNMT3A` <- (rowSums(dataList$Genetics[c('NPM1',"FLT3_ITD","DNMT3A")])==3)+0 ## Add NPM1:FLT3_ITD:DNMT3A product term as well
dataList$CytoCyto <- MakeInteractions(dataList$Cytogenetics, dataList$Cytogenetics)[,sapply(1:ncol(dataList$Cytogenetics), `<`, 1:ncol(dataList$Cytogenetics))]
dataList$CytoCyto <- dataList$CytoCyto[, colSums(dataList$CytoCyto, na.rm=TRUE) > 0]
dataList$GeneCyto <- MakeInteractions(dataList$Genetics, dataList$Cytogenetics)
dataList$GeneCyto <- dataList$GeneCyto[,colSums(dataList$GeneCyto, na.rm=TRUE) > 0]
dataList$GeneTreat <- MakeInteractions(dataList$Genetics, dataList$Treatment)
dataList$GeneTreat <- dataList$GeneTreat[,colSums(dataList$GeneTreat, na.rm=TRUE) > 0]
dataList$CytoTreat <- MakeInteractions(dataList$Cytogenetics, dataList$Treatment)
dataList$CytoTreat <- dataList$CytoTreat[,colSums(dataList$CytoTreat, na.rm=TRUE) > 0]

# Condensing to a data.frame
dataRaw <- do.call(cbind,dataList)
names(dataRaw) <- unlist(sapply(dataList, names))
dataFrame <- StandardizeMagnitude(dataRaw)
dim(dataFrame)

groups <- unlist(sapply(names(dataList), function(x) rep(x, ncol(dataList[[x]]))))
groups[grepl("^(t_)|(inv)", colnames(dataFrame)) &! grepl(":", colnames(dataFrame))] <- "Fusions"
groups[groups=="Cytogenetics"] <- "CNA"
groups <- factor(groups)
names(groups) <- colnames(dataFrame)
table(groups)

poorMansImpute <- function(x) {x[is.na(x)] <- mean(x, na.rm=TRUE); return(x)}
dataFrame <- as.data.frame(sapply(dataFrame, poorMansImpute))
rownames(dataFrame) <- clinicalData$PDID

# 1.3.5 Subclonal mutations
copyNumbers = cbind(dataList$Cytogenetics[grep(c("minus|plus|mono"), colnames(dataList$Cytogenetics))], clinicalData$gender)
copyNumbers$minus7 <- (copyNumbers$minus7 | copyNumbers$minus7q) +0
copyNumbers$minus7q <- NULL
for(i in 1:ncol(copyNumbers)){
  if(grepl("plus", colnames(copyNumbers)[i]))
    copyNumbers[,i] = copyNumbers[,i] * 3
}
copyNumbers[copyNumbers==0 | is.na(copyNumbers)] = 2
colnames(copyNumbers) = c(5,7,8,9,12,13,17,18,20,21,22,"Y",11,4,"X")
rownames(copyNumbers) <- clinicalData$PDID
copyNumbers$Y <- c(1:0)[clinicalData$gender] - mg14::na.zero(dataList$Cytogenetics$minusY)

cn = sapply(1:nrow(mutationData), function(i) {c=copyNumbers[mutationData$SAMPLE_NAME[i],match(mutationData$CHR[i], colnames(copyNumbers))]; if(length(c)==0) 2 else c})
vaf <- as.numeric(as.character(mutationData$X._MUT_IN_TUM))

depth <- as.numeric(as.character(mutationData$TUM_DEPTH))

dataFLT3_ITD <- read.table("./data/AMLSG_FLT3ITD.txt", sep="\t", header=TRUE)
dataFLT3_ITD$Sample <- sub("WGA_","", dataFLT3_ITD$Sample)
dataFLT3_ITD<-dataFLT3_ITD[dataFLT3_ITD$Sample%in%clinicalData$PDID,]

mcf <- vaf/100*cn ## Approx mutant cell fraction, assuming mutations on only one copy
mcf[which(mcf > 1.25)] <- vaf[which(mcf > 1.25)] ## Probably over adjusted
mcf[mcf > 1] <- 1 ## Random fluctuations
genesClonal <- dataFrame[groups=="Genetics"]
precedence <- matrix(0, nrow=ncol(genesClonal), ncol = ncol(genesClonal) , dimnames=list(colnames(genesClonal), colnames(genesClonal)))
plist <- list()
lesions <- as.character(mutationData$GENE)
lesions[mutationData$GENE=='IDH2' & grepl("172", mutationData$AA_CHANGE)] <- "IDH2_p172"
lesions[mutationData$GENE=='IDH2' & grepl("140", mutationData$AA_CHANGE)] <- "IDH2_p140"
lesions[mutationData$GENE=='FLT3' & grepl(paste(835:841, collapse="|"), mutationData$AA_CHANGE)] <- "FLT3_TKD"
lesions[mutationData$GENE=='FLT3' & grepl("ITD", mutationData$AA_CHANGE)] <- "FLT3_ITD"
lesions[lesions=="FLT3"] <- "FLT3_other"

# Add FLT3_ITD VAF, not the most accurate presumably, due to mapping problems for ITDs..
i <- lesions == "FLT3_ITD"
m <- match(mutationData$SAMPLE_NAME[i], dataFLT3_ITD$Sample)
mcf[i] <- as.numeric(as.character(dataFLT3_ITD$Read_count[m]))/dataFLT3_ITD$Coverage[m]

depth[i] <- dataFLT3_ITD$Coverage[m]

ix= lesions %in% colnames(precedence) & mutationData$Result %in% c("ONCOGENIC","POSSIBLE")
for(s in clinicalData$PDID){
  l <- list()
  for(i in which(mutationData$SAMPLE_NAME==s & ix))
    for(j in which(mutationData$SAMPLE_NAME==s & ix)){
      if(!is.na(cn[i]) & !is.na(cn[j]) & i!=j){
        m <- round(matrix(c(
          mcf[i]*depth[i],
          depth[i]-mcf[i]*depth[i],
          mcf[j]*depth[j],
          depth[j]-mcf[j]*depth[j]),
          ncol=2))
        f <- try(fisher.test(m, alternative="greater")$p.value< 0.01 , silent=TRUE) ## Fisher test
        if(class(f)!="try-error")
          if(f & mcf[i] >= 1 - mcf[j]){ ## Pidgeonhole
            precedence[as.character(lesions[i]),as.character(lesions[j])] <- precedence[as.character(lesions[i]),as.character(lesions[j])] + 1
            l <- c(l, list(c(as.character(lesions[i]),as.character(lesions[j]))))
            genesClonal[s, as.character(lesions[i])] <- 2
            genesClonal[s, as.character(lesions[j])] <- 3
          }
      }
    }
  plist[[s]] <- l
}

# 2 Models for overall survival
# 2.3.1 Number of oncogenic mutations
dataFrameOsTD <- dataFrame[tplSplitOs,]
dataFrameOsTD[which(tplIndexOs), grep("TPL", colnames(dataFrameOsTD), value=TRUE)] <- 0 ## Set pre-tpl variables to zero

# Define some indexes relating to subsets of variables used by the random effects model
mainGroups <- grep("[A-Z][a-z]+[A-Z]",levels(groups), invert=TRUE, value=TRUE)
mainGroups

mainIdx <- groups %in% mainGroups
osIdx <- !grepl("TPL", colnames(dataFrame)) ## Exclude TPL from OS analyses..
whichRFXOs <- which((colSums(dataFrame)>=8 | mainIdx) & osIdx) # ie, > 0.5%
mainIdxOs <- mainIdx & osIdx
osTDIdx <- !grepl("TPL_efs", colnames(dataFrame))
whichRFXOsTD <- which((colSums(dataFrame)>=8 | mainIdx) & osTDIdx) # ie, > 0.5%
mainIdxOsTD <- mainIdx & osTDIdx
whichRFXOsGG <- which((colSums(dataFrame)>=8 | mainIdxOs) & osIdx & groups %in% c(mainGroups,"GeneGene")) # ie, >0.5%

# 3.6.2 Prepare covariates
alloIdx <- clinicalData$TPL_type %in% c("ALLO","FREMD") # only allografts
alloTimeCR1 <- clinicalData$Time_1CR_TPL + .5 # +.5 to make > 0
alloTimeCR1[!alloIdx | (clinicalData$TPL_date < clinicalData$Recurrence_date & !clinicalData$TPL_Phase %in% c("CR1","RD1"))] <- NA

whichRFXOsTDGG <- which((colSums(dataFrame)>=8 | mainIdxOsTD) & osTDIdx & groups %in% c(mainGroups,"GeneGene")) #ie, > 0.5%
whichRFXRel <- whichRFXOsTDGG[grep("TPL",names(whichRFXOsTDGG), invert=TRUE)] #mainIdx & !grepl("TPL", names(dataFrame)) & groups!="Nuisance"
t <- clinicalData$Recurrence_date
t[is.na(t)] <- as.Date(1e6, origin="2000-01-01")

relData <- MakeTimeDependent(dataFrame[whichRFXRel], timeEvent=alloTimeCR1, timeStop=as.numeric(pmin(t, clinicalData$Date_LF) - clinicalData$CR_date), status=!is.na(clinicalData$Recurrence_date)+0)
relData$transplantCR1 <- relData$event
relData$event <- NULL
relData$transplantRel <- 0

nrdData <- MakeTimeDependent(dataFrame[whichRFXRel], timeEvent=alloTimeCR1, timeStop=as.numeric(pmin(t, clinicalData$Date_LF) - clinicalData$CR_date), status=is.na(clinicalData$Recurrence_date) & clinicalData$Status)
nrdData$transplantCR1 <- nrdData$event
nrdData$event <- NULL
nrdData$transplantRel <- 0

alloTimeRel <- clinicalData$TPL_date - clinicalData$Recurrence_date + .5 # +.5 to make > 0
alloTimeRel[!alloIdx | (clinicalData$TPL_date < clinicalData$Recurrence_date & !clinicalData$TPL_Phase %in% c("CR1","RD1"))] <- NA
i <- !is.na(clinicalData$Recurrence_date)
prdData <- MakeTimeDependent(dataFrame[i,whichRFXRel], timeEvent=alloTimeRel[i], timeStop=as.numeric(clinicalData$Date_LF- clinicalData$Recurrence_date)[i], status=clinicalData$Status[i])
prdData$transplantCR1 <- rep(0,nrow(prdData))
w <- sub("\\.1","",rownames(relData))[relData$status==1 & relData$transplantCR1==1]
prdData$transplantCR1[sub("\\.1","",rownames(prdData)) %in% w] <- 1
prdData$transplantRel <- prdData$event
prdData$event <- NULL
w <- which(prdData$time1 == prdData$time2) ## 5 cases with LF=Rec
prdData$time2[w] <- prdData$time2[w] + .5
prdData$time0 <- as.numeric(clinicalData$Recurrence_date[i]-clinicalData$CR_date[i])[prdData$index]

# 3.6.3 RFX fit of transitions
crGroups <- c(as.character(groups[whichRFXRel]), "Treatment","Treatment")
names(crGroups) <- c(names(dataFrame)[whichRFXRel],"transplantCR1","transplantRel")
coxRFXNrdTD <- CoxRFX(nrdData[names(crGroups)], Surv(nrdData$time1, nrdData$time2, nrdData$status), groups=crGroups, which.mu = intersect(mainGroups, unique(crGroups)))
coxRFXNrdTD$coefficients["transplantRel"] <- 0
#prsData$time1[!is.na(prsData$time1)] <- 0
coxRFXPrdTD <- CoxRFX(prdData[names(crGroups)], Surv(prdData$time1, prdData$time2, prdData$status), groups=crGroups, nu=1, which.mu = intersect(mainGroups, unique(crGroups)))
coxRFXRelTD <- CoxRFX(relData[names(crGroups)], Surv(relData$time1, relData$time2, relData$status), groups=crGroups, which.mu = intersect(mainGroups, unique(crGroups)))
coxRFXRelTD$coefficients["transplantRel"] <- 0

# 3.6.3.1 OS
osData <- MakeTimeDependent(dataFrame[whichRFXRel], timeEvent=alloTimeCR1, timeStop=as.numeric(clinicalData$Date_LF- clinicalData$CR_date), status=clinicalData$Status)
osData$transplantCR1 <- osData$event
osData$transplantRel <- osData$event
w <- which(clinicalData$TPL_date > clinicalData$Recurrence_date | clinicalData$TPL_Phase != "CR1")
osData$transplantCR1[osData$index %in% w] <- 0
osData$transplantRel[!osData$index %in% w] <- 0

data <- osData[rev(!duplicated(rev(osData$index))),colnames(coxRFXRelTD$Z)]
osData$transplantRel <- 0 # Note: confounded by relapse
rownames(data) <- sub("\\.1$","", rownames(data))
data <- data[rownames(dataFrame),]

# 3.6.6 Predicting outcome after CR
library(Rcpp)
cppFunction('NumericVector computeTotalPrsC(NumericVector x, NumericVector diffCir, NumericVector prsP, NumericVector tdPrmBaseline, double risk) {
				int xLen = x.size();
				double hj;
				double r = exp(risk);
				NumericVector rs(xLen);
				for(int i = 0; i < xLen; ++i) rs[i] = 1;
				for(int j = 1; j < xLen; ++j){ 
				hj = tdPrmBaseline[j-1] * r;
				for(int i = j; i < xLen; ++i){
				rs[i] += diffCir[j-1] * (1-pow(prsP[i-j], hj));
				}
				}
				return rs;
				}', rebuild=TRUE)

# Function to predict OS from CR1
MultiRFX3 <- function(coxRFXNrdTD, coxRFXRelTD, coxRFXPrdTD, data, x =365, ciType="analytical", prdData){
  ## Step 1: Compute KM survival curves and log hazard
  getS <- function(coxRFX, data, max.x=5000) {		
    if(!is.null(coxRFX$na.action)) coxRFX$Z <- coxRFX$Z[-coxRFX$na.action,]
    data <- as.matrix(data[,match(colnames(coxRFX$Z),colnames(data))])
    r <- PredictRiskMissing(coxRFX, data, var="var2")
    H0 <- basehaz(coxRFX, centered = FALSE)
    hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
    x <- c(0:max.x,max.x)
    S <- exp(-hazardDist(x))
    return(list(S=S, r=r, x=x, hazardDist=hazardDist, r0 = coxRFX$means %*% coef(coxRFX)))
  }
  kmRel <- getS(coxRFX = coxRFXRelTD, data = data, max.x=max(x))
  kmNrd <- getS(coxRFX = coxRFXNrdTD, data = data, max.x=max(x))
  kmPrd <- getS(coxRFX = coxRFXPrdTD, data = data, max.x=max(x))
  
  ## Step 2: Adjust CIR and NRM curve for competing risks, accounting for hazard
  kmRel$Sadj <- sapply(1:nrow(data), function(i) cumsum(c(1,diff(kmRel$S^exp(kmRel$r[i,1]))) * kmNrd$S ^ exp(kmNrd$r[i,1])))
  kmNrd$Sadj <- sapply(1:nrow(data), function(i) cumsum(c(1,diff(kmNrd$S^exp(kmNrd$r[i,1]))) * kmRel$S ^ exp(kmRel$r[i,1]))) ## array times x nrow(data)
  
  stopifnot(length(x)==1 | length(x) == nrow(data))
  if(length(x)==nrow(data))
    w <- match(x,kmRel$x)
  else if(length(x)==1)
    w <- rep(match(x, kmRel$x), nrow(data))
  y <- mapply(function(i,j) kmNrd$Sadj[i,j], w,1:length(w) ) # select time for each sample
  nrs <- y
  nrsUp <- y^exp(2*sqrt(kmNrd$r[,2]))
  nrsLo <- y^exp(- 2*sqrt(kmNrd$r[,2]))
  
  y <- mapply(function(i,j) kmRel$Sadj[i,j], w,1:length(w) ) # select time for each sample
  cir <- y
  cirLo <- y^exp( 2*sqrt(kmRel$r[,2]))
  cirUp <- y^exp( - 2*sqrt(kmRel$r[,2]))
  
  ## Step 3: Compute post-relapse survival
  survPredict <- function(surv){
    s <- survfit(surv~1)
    splinefun(s$time, s$surv, method="monoH.FC")
  }
  xx <- 0:max(x)
  # Baseline Prs (measured from relapse)
  kmPrs0 <- survPredict(Surv(prdData$time1, prdData$time2, prdData$status))(xx) 
  # PRS baseline with spline-based dep on CR length)
  
  coxphPrs <- coxph(Surv(time1, time2, status)~ pspline(time0, df=10), data=prdData ) 
  tdPrmBaseline <- exp(predict(coxphPrs, newdata=data.frame(time0=xx[-1])))						
  rs <- sapply(1:nrow(data), function(i){
    ### Different approach				
    xLen <- 1+floor(x)
    cir <- kmRel$Sadj[1:xLen,i]
    rs <- computeTotalPrsC(x = xx, diffCir = diff(cir), prsP = kmPrs0, tdPrmBaseline = tdPrmBaseline, risk = kmPrd$r[i,1]-kmPrd$r0)
    rs[xLen]
  })
  
  ## Step 4: Combine into overall survival
  if(any(1-(1-rs)-(1-nrs)<0)) warning("OS < 0 occured.")	
  os <- pmax(pmin(1-(1-rs)-(1-nrs),1),0)
  
  ## Step 5: Confidence intervals for OS
  osCi <- sapply(1:nrow(data), function(i){
    if("analytical" == ciType){
      ## Confidence intervals
      PlogP2 <- function(x) {(x * log(x))^2}
      errOs <- kmNrd$r[i,2] * PlogP2(kmNrd$S[w[i]]) * (1-kmRel$S[w[i]] * kmPrd$S[w[i]])^2 + kmRel$r[i,2]  * (1-kmNrd$S[w[i]])^2* kmPrd$S[w[i]]^2 * PlogP2(kmRel$S[w[i]]) +  kmPrd$r[i,2]  * (1-kmNrd$S[w[i]])^2* kmRel$S[w[i]]^2 * PlogP2(kmPrd$S[w[i]])
      errOs <- errOs / PlogP2(1-(1-kmNrd$S[w[i]])*(1-kmRel$S[w[i]]*kmPrd$S[w[i]]))
      return(c(osUp=os[i] ^ exp(-2* errOs), osLo= os[i] ^ exp(+2*errOs)))
    } else if("simulated" == ciType){
      ## Simulate CI
      nSim <- 200
      osCiMc <- sapply(1:nSim, function(foo){
        H <- exp(rnorm(3,c(kmRel$r[i,1],kmNrd$r[i,1],kmPrd$r[i,1]),sqrt(c(kmRel$r[i,2],kmNrd$r[i,2],kmPrd$r[i,2]))))
        nrs <- cumsum(c(1,diff(kmNrd$S^H[2]) * kmRel$S[-1]^H[1])) ## Correct KM estimate for competing risk
        diffCir <- diff(kmRel$inc^H[1]) * kmNrd$inc[-1]^H[2] ## Correct KM estimate for competing risk							
        rs <- computeTotalPrsC(x = x, diffCir = diffCir, prsP = kmPrs0, tdPrmBaseline = tdPrmBaseline, risk = -kmPrd$r0+log(H[3]))
        return((1-(1-nrs)-(1-rs))[w[i]])
      })
      osCiMcQ <- quantile(osCiMc, c(0.025,0.975))
      return(c(osUp = osCiMcQ[2], osLo = osCiMcQ[1]))
    }
  })
  
  return(data.frame(os=os, osLo = osCi[2,], osUp = osCi[1,],  cir=cir, cirLo=cirLo, cirUp=cirUp, nrs=nrs, nrsLo=nrsLo, nrsUp=nrsUp, rs=rs ))
}

# Create a data.frame with all data in cr
allData <- MakeTimeDependent(dataFrame[whichRFXRel], timeEvent=alloTimeCR1, timeStop=as.numeric(clinicalData$Date_LF- clinicalData$CR_date), status=clinicalData$Status)
allData$transplantCR1 <- allData$event
allData$transplantRel <- allData$event
w <- which(clinicalData$TPL_date > clinicalData$Recurrence_date)  
allData$transplantCR1[allData$index %in% w] <- 0
allData$transplantRel[!allData$index %in% w] <- 0

allDataTpl <- osData[rep(1:nrow(dataFrame), each=3),]
allDataTpl$transplantCR1 <- rep(c(0,1,0), nrow(dataFrame))
allDataTpl$transplantRel <- rep(c(0,0,1), nrow(dataFrame))

# 3.6.6.3 Allogeneic hematopoietic stem cell transplants
multiRFX3Tpl <- MultiRFX3(coxRFXNrdTD, coxRFXRelTD, coxRFXPrdTD, data=allDataTpl, x=3*365, prdData=prdData)
multiRFX3Tpl <- as.data.frame(matrix(multiRFX3Tpl$os, ncol=3, byrow=TRUE, dimnames=list(NULL, c("None","CR1","Relapse"))), row.names=rownames(dataFrame))

MultiRFX3TplCi <- function(coxRFXNrdTD, coxRFXRelTD, coxRFXPrdTD, data, x =365, prdData, ciType="simulated", nSim = 200, mc.cores=10){
  ## Step 1: Compute KM survival curves and log hazard
  getS <- function(coxRFX, data, max.x=5000) {		
    if(!is.null(coxRFX$na.action)) coxRFX$Z <- coxRFX$Z[-coxRFX$na.action,]
    data <- as.matrix(data[,match(colnames(coxRFX$Z),colnames(data))])
    r <- PredictRiskMissing(coxRFX, data, var="var2")
    H0 <- basehaz(coxRFX, centered = FALSE)
    hazardDist <- splinefun(H0$time, H0$hazard, method="monoH.FC")
    x <- c(0:max.x,max.x)
    S <- exp(-hazardDist(x))
    return(list(S=S, r=r, x=x, hazardDist=hazardDist, r0 = coxRFX$means %*% coef(coxRFX)))
  }
  
  data$transplantCR1 <- 0
  data$transplantRel <- 0
  
  kmCir <- getS(coxRFX = coxRFXRelTD, data = data, max.x=max(x))
  kmNrs <- getS(coxRFX = coxRFXNrdTD, data = data, max.x=max(x))
  kmPrs <- getS(coxRFX = coxRFXPrdTD, data = data, max.x=max(x))
  
  survPredict <- function(surv){
    s <- survfit(surv~1)
    splinefun(s$time, s$surv, method="monoH.FC")
  }
  xx <- 0:max(x)
  
  # Baseline Prs (measured from relapse)
  kmPrs0 <- survPredict(Surv(prdData$time1, prdData$time2, prdData$status))(xx) 
  
  # PRS baseline with spline-based dep on CR length)
  coxphPrs <- coxph(Surv(time1, time2, status)~ pspline(time0, df=10), data=prdData) 
  tdPrmBaseline <- exp(predict(coxphPrs, newdata=data.frame(time0=xx[-1])))						
  
  stopifnot(length(x)==1 | length(x) == nrow(data))
  if(length(x)==nrow(data))
    w <- match(x,kmCir$x)
  else if(length(x)==1)
    w <- rep(match(x, kmCir$x), nrow(data))
  
  survival <- sapply(c("None","Rel","CR1"), function(type){
    if(type=="None"){
      data$transplantCR1 <- 0
      data$transplantRel <- 0
    }else if(type=="Rel"){
      data$transplantCR1 <- 0
      data$transplantRel <- 1					
    }else if(type=="CR1"){
      data$transplantCR1 <- 1
      data$transplantRel <- 0
    }
    
    
    kmCir <- getS(coxRFX = coxRFXRelTD, data = data, max.x=max(x))
    kmNrm <- getS(coxRFX = coxRFXNrdTD, data = data, max.x=max(x))
    kmPrs <- getS(coxRFX = coxRFXPrdTD, data = data, max.x=max(x))
    
    ## Step 2: Adjust CIR and NRM curve for competing risks, accounting for hazard
    kmCir$Sadj <- sapply(1:nrow(data), function(i) cumsum(c(1,diff(kmCir$S^exp(kmCir$r[i,1]))) * kmNrm$S ^ exp(kmNrm$r[i,1])))
    kmNrm$Sadj <- sapply(1:nrow(data), function(i) cumsum(c(1,diff(kmNrm$S^exp(kmNrm$r[i,1]))) * kmCir$S ^ exp(kmCir$r[i,1]))) ## array times x nrow(data)
    
    y <- mapply(function(i,j) kmNrm$Sadj[i,j], w,1:length(w) ) # select time for each sample
    nrs <- y
    nrsUp <- y^exp(2*sqrt(kmNrm$r[,2]))
    nrsLo <- y^exp(- 2*sqrt(kmNrm$r[,2]))
    
    y <- mapply(function(i,j) kmCir$Sadj[i,j], w,1:length(w) ) # select time for each sample
    cir <- y
    cirLo <- y^exp( 2*sqrt(kmCir$r[,2]))
    cirUp <- y^exp( - 2*sqrt(kmCir$r[,2]))
    
    ## Step 3: Compute post-relapse survival
    rs <- sapply(1:nrow(data), function(i){
      ### Different approach				
      xLen <- 1+floor(x)
      cir <- kmCir$Sadj[1:xLen,i]
      rs <- computeTotalPrsC(x = xx, diffCir = diff(cir), prsP = kmPrs0, tdPrmBaseline = tdPrmBaseline, risk = kmPrs$r[i,1]-kmPrs$r0)
      rs[xLen]
    })
    
    ## Step 4: Combine into overall survival
    if(any(1-(1-rs)-(1-nrs)<0)) warning("OS < 0 occured.")	
    os <- pmax(pmin(1-(1-rs)-(1-nrs),1),0)
    cbind(os, rs, nrs, aar=rs-cir)
  }, simplify='array')
  
  ## Step 5: Confidence intervals for OS
  osCi <- sapply(mclapply(1:nrow(data), function(i){
    {
      ## Simulate CI
      osCiMc <- sapply(1:nSim, function(foo){
        r0 <- rnorm(3,c(kmCir$r[i,1],kmNrs$r[i,1],kmPrs$r[i,1]),sqrt(c(kmCir$r[i,2],kmNrs$r[i,2],kmPrs$r[i,2])))
        H0 <- exp(r0)
        nrs0 <- cumsum(c(1,diff(kmNrs$S^H0[2])) * kmCir$S^H0[1]) ## Correct KM estimate for competing risk
        diffCir <- diff(c(1,kmCir$S^H0[1])) * kmNrs$S^H0[2] ## Correct KM estimate for competing risk			
        cir0 <- 1+cumsum(diffCir)
        rs0 <- computeTotalPrsC(x = xx, diffCir = diffCir, prsP = kmPrs0, tdPrmBaseline = tdPrmBaseline, risk = -kmPrs$r0+log(H0[3]))
        aar0 <- rs0[1:w[i]]-cir0[1:w[i]]
        
        Hcr1 <- exp(r0 + rnorm(3,c(coxRFXRelTD$coefficients["transplantCR1"],coxRFXNrdTD$coefficients["transplantCR1"],coxRFXPrdTD$coefficients["transplantCR1"]), 
                               sqrt(c(coxRFXRelTD$var2["transplantCR1","transplantCR1"],coxRFXNrdTD$var2["transplantCR1","transplantCR1"],coxRFXPrdTD$var2["transplantCR1","transplantCR1"])))) 
        nrsCr1 <- cumsum(c(1,diff(kmNrs$S^Hcr1[2])) * kmCir$S^Hcr1[1]) ## Correct KM estimate for competing risk
        diffCir <- diff(c(1,kmCir$S^Hcr1[1])) * kmNrs$S^Hcr1[2] ## Correct KM estimate for competing risk	
        cirCr1 <- 1+cumsum(diffCir)
        rsCr1 <- computeTotalPrsC(x = xx, diffCir = diffCir, prsP = kmPrs0, tdPrmBaseline = tdPrmBaseline, risk = -kmPrs$r0+log(Hcr1[3]))
        aarCr1 <- rsCr1[1:w[i]]-cirCr1[1:w[i]]
        
        
        Hrel <- exp(r0 + rnorm(3,c(coxRFXRelTD$coefficients["transplantRel"],coxRFXNrdTD$coefficients["transplantRel"],coxRFXPrdTD$coefficients["transplantRel"]), 
                               sqrt(c(coxRFXRelTD$var2["transplantRel","transplantRel"],coxRFXNrdTD$var2["transplantRel","transplantRel"],coxRFXPrdTD$var2["transplantRel","transplantRel"])))) 
        nrsRel <- cumsum(c(1,diff(kmNrs$S^Hrel[2])) * kmCir$S^Hrel[1]) ## Correct KM estimate for competing risk
        diffCir <- diff(c(1,kmCir$S^Hrel[1])) * kmNrs$S^Hrel[2] ## Correct KM estimate for competing risk							
        cirRel <- 1+cumsum(diffCir)
        rsRel <- computeTotalPrsC(x = xx, diffCir = diffCir, prsP = kmPrs0, tdPrmBaseline = tdPrmBaseline, risk = -kmPrs$r0+log(Hrel[3]))
        aarRel <- rsRel[1:w[i]]-cirRel[1:w[i]]
        
        
        os0 <- (1-(1-nrs0[1:w[i]])-(1-rs0))[w[i]]
        osCr1 <- (1-(1-nrsCr1[1:w[i]])-(1-rsCr1))[w[i]]
        osRel <- (1-(1-nrsRel[1:w[i]])-(1-rsRel))[w[i]]
        return(cbind(os=c(none=os0, cr1=osCr1, rel=osRel, dCr1=osCr1-os0, dRel=osRel-os0, dCr1Rel=osCr1-osRel),
                     rs=c(none=rs0[w[i]], cr1=rsCr1[w[i]], rel=rsRel[w[i]], dCr1=rsCr1[w[i]]-rs0[w[i]], dRel=rsRel[w[i]]-rs0[w[i]], dCr1Rel=rsCr1[w[i]]-rsRel[w[i]]),
                     nrs=c(none=nrs0[w[i]], cr1=nrsCr1[w[i]], rel=nrsRel[w[i]], dCr1=nrsCr1[w[i]]-nrs0[w[i]], dRel=nrsRel[w[i]]-nrs0[w[i]], dCr1Rel=nrsCr1[w[i]]-nrsRel[w[i]]),
                     aar=c(none=aar0[w[i]], cr1=aarCr1[w[i]], rel=aarRel[w[i]], dCr1=aarCr1[w[i]]-aar0[w[i]], dRel=aarRel[w[i]]-aar0[w[i]], dCr1Rel=aarCr1[w[i]]-aarRel[w[i]])))
      }, simplify='array')
      osCiMcQ <- apply(osCiMc,1:2,quantile, c(0.025,0.5,0.975))
      return(sapply(c("os","rs","nrs","aar"), function(t) 
        cbind(hat = c(survival[i,t,1], survival[i,t,3], survival[i,t,2], survival[i,t,3]-survival[i,t,1], survival[i,t,2]-survival[i,t,1], survival[i,t,3]-survival[i,t,2]), 
              median = osCiMcQ[2,,t], lower = osCiMcQ[1,,t], upper = osCiMcQ[3,,t]), simplify="array"))
    }
  }, mc.cores=mc.cores), I, simplify="array")
  #cat(os, "\n")
  return(osCi)
}


# 3.6.6.4 Leave one out cross-validation
# 3.6.6.4.1 Three state model

multiRFX3TplCiLoo <- sapply(mclapply(rownames(dataFrame), function(pd){
  i <- which(rownames(dataFrame)==pd)
  whichTrain <<- which(rownames(dataFrame)!=pd)
  rfxNrs <- CoxRFX(nrdData[nrdData$index %in% whichTrain, names(crGroups)], Surv(nrdData$time1, nrdData$time2, nrdData$status)[nrdData$index %in% whichTrain], groups=crGroups, which.mu = intersect(mainGroups, unique(crGroups)))
  rfxNrs$coefficients["transplantRel"] <- 0
  #prsData$time1[!is.na(prsData$time1)] <- 0
  rfxPrs <-  CoxRFX(prdData[prdData$index %in% whichTrain, names(crGroups)], Surv(prdData$time1, prdData$time2, prdData$status)[prdData$index %in% whichTrain], groups=crGroups, nu=1, which.mu = intersect(mainGroups, unique(crGroups)))
  rfxRel <-  CoxRFX(relData[relData$index %in% whichTrain, names(crGroups)], Surv(relData$time1, relData$time2, relData$status)[relData$index %in% whichTrain], groups=crGroups, which.mu = intersect(mainGroups, unique(crGroups)))
  rfxRel$coefficients["transplantRel"] <- 0
  multiRFX3TplCi <- MultiRFX3TplCi(rfxNrs, rfxRel, rfxPrs, data=data[i,colnames(rfxPrs$Z), drop=FALSE], x=3*365, nSim=200, prdData=prdData[prdData$index!=i,], mc.cores=1)
}, mc.cores=10), I, simplify="array")[,,,1,]

multiRFX3TplLoo <- t(multiRFX3TplCiLoo[1:3,"hat","os",])

pdf(file = "./Figure3.pdf")

# Fig.3
riskCol=set1[c(1,3,4,2)]
names(riskCol) <- levels(clinicalData$M_Risk)
nSim <- 200
par(mar=c(3,3,1,3), las=2, mgp=c(2,.5,0), bty="n")
benefit <- multiRFX3TplLoo[,2]-multiRFX3TplLoo[,3]
benefitGroup <- factor(benefit > 0.1, labels=c("Low","High"))
absrisk <- multiRFX3TplLoo[,1]
names(absrisk) <- names(benefit) <- rownames(dataFrame)
s <- clinicalData$AOD < 60 & !is.na(clinicalData$CR_date) &! clinicalData$TPL_Phase %in% c("RD1","PR1")
x <- 1-absrisk
y <- benefit
plot(x[s], y[s], pch=NA, ylab="Mortality reduction from CR1 allograft compared to salvage allograft", xlab="3yr mortality with standard chemo only (995 patients)", col=riskCol[clinicalData$M_Risk], cex=.8, las=1, ylim=range(benefit))
abline(h=seq(-.1,.3,.1), col='grey', lty=3)
abline(v=seq(.2,.9,0.2), col='grey', lty=3)
points(x[s], y[s], pch=16,  col=riskCol[clinicalData$M_Risk[s]], cex=.8)
# Add loess fit, accounting for correlations of errors
xn <- seq(0.01,0.99,0.01)
fit <- sapply(1:nSim, function(i){
  #benefit <- multiRFX3TplCiCorLoo[2,"os",i,]-multiRFX3TplCiCorLoo[3,"os",i,]
  #absrisk <- multiRFX3TplCiCorLoo[1,"os",i,]
  s <- clinicalData$AOD < 60 & !is.na(clinicalData$CR_date) &! clinicalData$TPL_Phase %in% c("RD1","PR1")
  x <- 1-absrisk
  y <- benefit
  p <- predict(loess(y~x, data=data.frame(x=x[s], y=y[s])), newdata=data.frame(x=xn), se=TRUE)
  yn <- c(p$fit + 2*p$se.fit,rev(p$fit - 2*p$se.fit))
  p$fit
})
q <- apply(fit, 1, quantile, c(0.025, 0.5, 0.975), na.rm=TRUE)
#polygon(c(xn, rev(xn)), c(q[1,], rev(q[3,])), border=NA, col="#00000044", lwd=1)
lines(xn, rowMeans(fit, na.rm=TRUE))
legend("topleft", pch=c(16,16,16,16,NA),lty=c(NA,NA,NA,NA,1), col=c(riskCol[c(2,4,3,1)],1),fill=c(NA,NA,NA,NA,NA), border=NA, c(names(riskCol)[c(2,4,3,1)],"loess average"), box.lty=0)
n <- c(100,50,20,10,5,4,3)
axis(side=4, at=1/n, labels=n, las=1)
mtext("Number needed to treat", side=4, at=.2, line=2, las=0)
axis(side=4, at=-1/n, labels=n, las=1)
mtext("Number needed to harm", side=4, at=-.1, line=2, las=0)

dev.off()
