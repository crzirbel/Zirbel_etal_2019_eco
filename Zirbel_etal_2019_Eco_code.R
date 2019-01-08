#=================================================================#
#Code for Zirbel et al. (2019) "Landscape context                 #
#explains ecosystem multifunctionality in restored grasslands"    #
#better than plant diversity                                      #
#Created by: Chad Zirbel                                          #
#11/4/16                                                          #
#=================================================================#

###load required libraries###
library(vegan)
library(reshape2)
library(FD)
library(PhyloMeasures)
library(pez)
library(lme4)
library(car)
library(ggplot2)
library(MuMIn)
library(r2glmm)
library(data.table)
library(mgcv)

##load data##
#Environmental data
site_data<-read.csv("http://datadryad.org/bitstream/handle/10255/dryad.135949/site_environmental_Zirbel_2017_dryad.csv?sequence=1")

#Species composition
species_abund<-read.csv("http://datadryad.org/bitstream/handle/10255/dryad.135948/plant_community_Zirbel_2017_dryad.csv?sequence=1")

#species composition vs. phylogeny name table
tax_phylo<-read.csv("https://datadryad.org/bitstream/handle/10255/dryad.135952/tax_phylo_names_Zirbel_2017_dryad.csv?sequence=1")

#Functional traits
trait_con<-read.csv("http://datadryad.org/bitstream/handle/10255/dryad.135950/trait_data_field_Zirbel_2017_dryad.csv?sequence=1")
trait_lit<-read.csv("http://datadryad.org/bitstream/handle/10255/dryad.135951/trait_data_literature_Zirbel_2017_dryad.csv?sequence=1")

#Phylogeny
phylo<-read.tree("https://datadryad.org/bitstream/handle/10255/dryad.135953/phylogeny_Zirbel_2017_dryad.txt?sequence=1")

#Ecosystem functions
functions.byplot<-read.csv("http://datadryad.org/bitstream/handle/10255/dryad.135621/ecosystem_function_Zirbel_2017_dryad.csv?sequence=1")

#####Environmental data#####
#site age
site_data$age.2013<-2013-site_data$year

#Fire return interval
fire.return.2013<-site_data$no.burns/(site_data$age.2013)
site_data$fire.return<-fire.return.2013

#Soil PCA
soildata<-site_data
soil<-soildata[8:15]
pca.soil<-prcomp(soil, center=TRUE, scale=TRUE)#scale=TRUE, different scales across variables 
soilPC1<-pca.soil$x[,1]#56% variance; high values=sandy sites with low WHC (P loads high)

#invert the soil PCA axis (increasing WHC) and add to the data frame
site_data$soilPC1<-soilPC1*-1

#landscape PCA
land <- site_data[3:7]
pca.land<-prcomp(land, center=TRUE, scale=FALSE)
landPC1<-pca.land$x[,1]#59% of variance; high values = high forest/grass, low ag

site_data$landscape<-landPC1

#cut to relevant variables
site_data<-site_data[c(1,17:21)]

#####Taxonomic Diversity#####
rownames(species_abund)<-paste(species_abund$site.id,species_abund$plot,sep=".")
species_abund<-species_abund[-c(1:2)]

sim.div<-diversity(species_abund, index = "simpson")
div<-sim.div

#####Functional Diversity#####
#calculate average trait values for each species
trait_con<-trait_con[,1:4] #drop seed mass collection column
trait.melt<-melt(trait_con, id="species",na.rm=T)
trait.cast=dcast(trait.melt, species~variable,mean,drop=F)

#merge field collected and literature traits
trait_data<-merge(trait_lit,trait.cast, by="species",all.x=T)
rownames(trait_data)<-trait_data[,1]
trait_data<-trait_data[,-1]
trait_data[ is.na(trait_data) ] <- NA

#caluculate functional diversity
fd<-dbFD(gowdis(trait_data),species_abund,cor= "lingoes",stand.FRic = T)

#Add to data.frame
div<-cbind(div,fd$RaoQ)
colnames(div)<-c("sim.div","FRaoQ")

#####Phylogenetic Diversity#####
#build a three speices tree that adds in Equisetum and Botyrchium which were not in the Zanne et al. (2014) tree
#branch lengths are from timetree.org
pter.tree<-read.tree(text="(Botrychium_sp: 372.6,(Onoclea_sensibilis: 352.3 , Equisetum_sp: 352.3):20.3);")

#plot the tree to see if it looks right
plot.phylo(pter.tree,show.tip.label = TRUE, edge.width = .8, cex = 1)
add.scale.bar(length=372.6)

# replace the Ono_sen tip with the pter.tree
phylo.add<-bind.replace(phylo, pter.tree, replacing.tip.label="Onoclea_sensibilis", donor.length = (372.6*2))

#plot the full tree
plot.phylo(phylo.add,type="fan",show.tip.label = TRUE, cex = .5)

#format species abundance matrix to have same species names as the phylogeny
species_abund.p<-species_abund
colnames(species_abund.p)<-tax_phylo$phylo.name

#calculate phylogenetic diversity
ses.pd<-pd.query(phylo.add,species_abund.p,standardize = T)

div<-cbind(div,ses.pd)
colnames(div)<-c("sim.div","FRaoQ","ses.pd")

#add plot and site.id columns to div df
div<-as.data.frame(div)
div$site.id<-sub('\\..*', '', rownames(div))
div$plot<-sub('.*\\.', '', rownames(div))

#####Ecosystem Functions#####
#averaging floral cover and SR across three measurements:
functions.byplot$floral.cover<-(functions.byplot$floral.cover.june + functions.byplot$floral.cover.july + functions.byplot$floral.cover.sept)/3
functions.byplot$floral.plotSR<-(functions.byplot$floral.plotSR.june + functions.byplot$floral.plotSR.july + functions.byplot$floral.plotSR.sept)/3
functions.byplot$floral.siteSR<-(functions.byplot$floral.siteSR.june + functions.byplot$floral.siteSR.july + functions.byplot$floral.siteSR.sept)/3

#wide for analysis
melt.functions.byplot<-melt(functions.byplot, id=c("site.id", "plot"))
mean.functions<-dcast(melt.functions.byplot, site.id~variable, mean, na.rm=TRUE)
var.functions<-dcast(melt.functions.byplot, site.id~variable, var, na.rm=TRUE)

#long for plotting with facets
mean.functions.l<-dcast(melt.functions.byplot, site.id+variable~., mean, na.rm=TRUE)
names(mean.functions.l)<-c("site.id", "funct", "mean.function")
var.functions.l<-dcast(melt.functions.byplot, site.id+variable~., var, na.rm=TRUE)
names(var.functions.l)<-c("site.id", "funct", "var.function")

#standardizing functions by the plot max (as recommended by byrnes et al)
relevant.fun<-melt.functions.byplot[melt.functions.byplot$variable %in% c("biomass", "mass.loss", "roots.per.g.sand", "total.seeds.removed", "aug.worms.gone", "floral.cover","floral.plotSR"),]
max.fun<-dcast(relevant.fun, variable~., fun=max, na.rm=T) #finds highest plot value from any site (universal) for each variable;
names(max.fun)<-c("variable", "universal.max")
stand.functions<-merge(relevant.fun, max.fun, by="variable")
stand.functions$prop.max<-stand.functions$value/stand.functions$universal.max

##Calculate average multifunctionality
#averaging across functions in a plot:
stand.functions.byplot<-dcast(stand.functions, site.id+plot~., value.var="prop.max", fun=mean, na.rm=T)
names(stand.functions.byplot)<-c("site.id", "plot", "multifunct.a")

##THRESHOLD APPROACH
melt.stand.functions.byplot<-dcast(stand.functions, site.id+plot+variable~., value.var="prop.max", fun=mean, na.rm=T) #finds site mean of each standardized function
names(melt.stand.functions.byplot)=c("site.id","plot", "funct", "site.prop.max")

melt.stand.functions.byplot$multifunct.t25<-1*(melt.stand.functions.byplot$site.prop.max>0.25) #is the function performing at >25% of max?
melt.stand.functions.byplot$multifunct.t50<-1*(melt.stand.functions.byplot$site.prop.max>0.50)
melt.stand.functions.byplot$multifunct.t75<-1*(melt.stand.functions.byplot$site.prop.max>0.75)
multifunct.t<-aggregate(.~site.id+plot,data=melt.stand.functions.byplot,sum) #counts number of functions performing above the threshold

#clean threshold data
multifunct.t<-multifunct.t[,-3:-4]

####Combine data.frames####
full.data<-merge(stand.functions.byplot,div,by= c("site.id","plot"),all=T)
full.data<-merge(full.data,multifunct.t, by=c("site.id","plot"),all=T)
full.data<-merge(full.data,site_data,by="site.id",all.x=T)

#add in relevant ecosystem functions
functions.byplot.rel<-functions.byplot[,c(1:5,11:12,22:23)]
full.data<-merge(full.data,functions.byplot.rel,by=c("site.id","plot"),all.x=T)

#logit transform floral cover data
min(full.data$floral.cover[full.data$floral.cover!=0])#minimum non-zero values to be added before transformation
full.data$fc.logit<-logit(full.data$floral.cover,adjust=0.003333333)

#add andger abundance
andger<-species_abund["Andropogon.gerardii"]
andger$site.id<-sub('\\..*', '', rownames(andger))
andger$plot<-sub('.*\\.', '', rownames(andger))
full.data<-merge(full.data,andger,by= c("site.id","plot"),all=T)

#scale variables
data.scale<-scale(full.data[3:23],center=T,scale=T)
data.scale<-cbind(full.data[1:2],data.scale)

#Site average data.frame
full.data.av<-aggregate(full.data,by=list(full.data$site),mean,na.rm=T)

data.scale.av<-scale(full.data.av[3:23],center=T,scale=T)
data.scale.av<-cbind(full.data.av[1:2],data.scale)

####Model Selection with AICc####
#REML=F in all models

##Multifunctionality
#Averaging approach (Maestre et al. 2012)

#Phylogenetic
mod.mult.phylo.aic<-lmer(multifunct.a~ses.pd+soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii+
                (1|site.id),data=data.scale,REML=F)

#Functional
mod.mult.func.aic<-lmer(multifunct.a~FRaoQ+soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii+
               (1|site.id),data=data.scale,REML=F)

#Taxonomic
mod.mult.tax.aic<-lmer(multifunct.a~sim.div+soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii+
              (1|site.id),data=data.scale,REML=F)

#Null model
#not included in results but good reassurance
mod.mult.null.aic<-lmer(multifunct.a~1+(1|site.id),data=data.scale,REML=F)

#model comparison
AICc(mod.mult.phylo.aic,mod.mult.func.aic,mod.mult.tax.aic,mod.mult.null.aic)

##Threshold models (Gamfeldt et al. 2008)
#25% Threshold

#Phylogenetic
mod.t25.phylo.aic<-lmer(multifunct.t25~ses.pd+soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii+
                           (1|site.id),data=data.scale,REML=F)

#Functional
mod.t25.func.aic<-lmer(multifunct.t25~FRaoQ+soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii+
                          (1|site.id),data=data.scale,REML=F)

#Taxonomic
mod.t25.tax.aic<-lmer(multifunct.t25~sim.div+soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii+
                         (1|site.id),data=data.scale,REML=F)

#Null model
#not included in results but good reassurance
mod.t25.null.aic<-lmer(multifunct.t25~1+(1|site.id),data=data.scale,REML=F)

#model comparison
AICc(mod.t25.phylo.aic,mod.t25.func.aic,mod.t25.tax.aic,mod.t25.null.aic)

#50% Threshold
#Phylogenetic
mod.t50.phylo.aic<-lmer(multifunct.t50~ses.pd+soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii+
                          (1|site.id),data=data.scale,REML=F)

#Functional
mod.t50.func.aic<-lmer(multifunct.t50~FRaoQ+soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii+
                         (1|site.id),data=data.scale,REML=F)

#Taxonomic
mod.t50.tax.aic<-lmer(multifunct.t50~sim.div+soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii+
                        (1|site.id),data=data.scale,REML=F)

#Null model
#not included in results but good reassurance
mod.t50.null.aic<-lmer(multifunct.t50~1+(1|site.id),data=data.scale,REML=F)

#model comparison
AICc(mod.t50.phylo.aic,mod.t50.func.aic,mod.t50.tax.aic,mod.t50.null.aic)

#75% Threshold
#Phylogenetic
mod.t75.phylo.aic<-lmer(multifunct.t75~ses.pd+soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii+
                          (1|site.id),data=data.scale,REML=F)

#Functional
mod.t75.func.aic<-lmer(multifunct.t75~FRaoQ+soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii+
                         (1|site.id),data=data.scale,REML=F)

#Taxonomic
mod.t75.tax.aic<-lmer(multifunct.t75~sim.div+soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii+
                        (1|site.id),data=data.scale,REML=F)

#Null model
#not included in results but good reassurance
mod.t75.null.aic<-lmer(multifunct.t75~1+(1|site.id),data=data.scale,REML=F)

#model comparison
AICc(mod.t75.phylo.aic,mod.t75.func.aic,mod.t75.tax.aic,mod.t75.null.aic)

###Individual ecosystem function models###

##Biomass
#Phylogenetic
mod.bio.phylo.aic<-lmer(biomass~ses.pd+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                      (1|site.id),data=data.scale,REML=F)

#Functional
mod.bio.func.aic<-lmer(biomass~FRaoQ+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                     (1|site.id),data=data.scale,REML=F)

#Taxonomic
mod.bio.tax.aic<-lmer(biomass~sim.div+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                    (1|site.id),data=data.scale,REML=F)

#Null model
mod.bio.null<-lmer(biomass~1+(1|site.id),data=data.scale,REML=F)

AICc(mod.bio.phylo.aic,mod.bio.func.aic,mod.bio.tax.aic,mod.bio.null)

##decomposition
#Phylo
mod.decomp.phylo.aic<-lmer(mass.loss~ses.pd+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                         (1|site.id),data=data.scale,REML=F)

#Functional
mod.decomp.func.aic<-lmer(mass.loss~FRaoQ+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                        (1|site.id),data=data.scale,REML=F)

#Taxonomic
mod.decomp.tax.aic<-lmer(mass.loss~sim.div+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                       (1|site.id),data=data.scale,REML=F)

#Null model
mod.decomp.null<-lmer(mass.loss~1+(1|site.id),data=data.scale,REML=F)

AICc(mod.decomp.phylo.aic,mod.decomp.func.aic,mod.decomp.tax.aic,mod.decomp.null)

##Belowground biomass
#Phylo
mod.root.phylo.aic<-lmer(roots.per.g.sand~ses.pd+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                       (1|site.id),data=data.scale,REML=F)

#Functional
mod.root.func.aic<-lmer(roots.per.g.sand~FRaoQ+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                      (1|site.id),data=data.scale,REML=F)

#Taxonomic
mod.root.tax.aic<-lmer(roots.per.g.sand~sim.div+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                     (1|site.id),data=data.scale,REML=F)

#Null model
mod.root.null<-lmer(roots.per.g.sand~1+(1|site.id),data=data.scale,REML=F)

AICc(mod.root.phylo.aic,mod.root.func.aic,mod.root.tax.aic,mod.root.null)

#Seed predation
#Phylo
mod.seed.phylo.aic<-lmer(total.seeds.removed~ses.pd+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                       (1|site.id),data=data.scale,REML=F)

#Functional
mod.seed.func.aic<-lmer(total.seeds.removed~FRaoQ+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                      (1|site.id),data=data.scale,REML=F)

#Taxonomic
mod.seed.tax.aic<-lmer(total.seeds.removed~sim.div+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                     (1|site.id),data=data.scale,REML=F)

#Null model
mod.seed.null<-lmer(total.seeds.removed~1+(1|site.id),data=data.scale,REML=F)

AICc(mod.seed.phylo.aic,mod.seed.func.aic,mod.seed.tax.aic,mod.seed.null)

##Arthropod predation
#Phylo
mod.arthro.phylo.aic<-lmer(aug.worms.gone~ses.pd+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                         (1|site.id),data=data.scale,REML=F)

#Functional
mod.arthro.func.aic<-lmer(aug.worms.gone~FRaoQ+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                        (1|site.id),data=data.scale,REML=F)

#Taxonomic
mod.arthro.tax.aic<-lmer(aug.worms.gone~sim.div+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                       (1|site.id),data=data.scale,REML=F)

#Null model
mod.arthro.null<-lmer(aug.worms.gone~1+(1|site.id),data=data.scale,REML=F)

AICc(mod.arthro.phylo.aic,mod.arthro.func.aic,mod.arthro.tax.aic,mod.arthro.null)


##Floral cover
#Phylo
mod.fc.phylo.aic<-lmer(fc.logit~ses.pd+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                   (1|site.id),data=data.scale,REML=F)

#Functional
mod.fc.func.aic<-lmer(fc.logit~FRaoQ+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                    (1|site.id),data=data.scale,REML=F)

#Taxonomic
mod.fc.tax.aic<-lmer(fc.logit~sim.div+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                   (1|site.id),data=data.scale,REML=F)

#Null model
mod.fc.null<-lmer(fc.logit~1+(1|site.id),data=data.scale,REML=F)

AICc(mod.fc.phylo.aic,mod.fc.func.aic,mod.fc.tax.aic,mod.fc.null)


#Floral richness
#Phylo
mod.fr.phylo.aic<-lmer(floral.plotSR~ses.pd+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                     (1|site.id),data=data.scale,REML=F)


#Functional
mod.fr.func.aic<-lmer(floral.plotSR~FRaoQ+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                    (1|site.id),data=data.scale,REML=F)

#Taxonomic
mod.fr.tax.aic<-lmer(floral.plotSR~sim.div+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                   (1|site.id),data=data.scale,REML=F)

#Null model
mod.fr.null<-lmer(floral.plotSR~1+(1|site.id),data=data.scale,REML=F)

AICc(mod.fr.phylo.aic,mod.fr.func.aic,mod.fr.tax.aic,mod.fr.null)

####Calculating regression coefs and p-values for multifunctionality and individual function models####
#Calculate full and partial R2 for best models
#REML=T for all models

#Variance inflation factor (VIF) function for glmms
vif.glmm <- function (fit) {
  # adapted from rms::vif
  v <- vcov(fit)
  nam <- names(fixef(fit))
  # exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}
##Multifunctionality models##
#averaging approach

#Phylogenetic
mod.mult.phylo<-lmer(multifunct.a~ses.pd+soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii+
                (1|site.id),data=data.scale,REML=T)
summary(mod.mult.phylo)
Anova(mod.mult.phylo,type=3)
vif.glmm(mod.mult.phylo)

#Functional
mod.mult.func<-lmer(multifunct.a~FRaoQ+soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii+
               (1|site.id),data=data.scale,REML=T)
summary(mod.mult.func)
Anova(mod.mult.func)
vif.glmm(mod.mult.func)

#Taxonomic
mod.mult.tax<-lmer(multifunct.a~sim.div+soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii+
              (1|site.id),data=data.scale,REML=T)
summary(mod.mult.tax)
Anova(mod.mult.tax)
vif.glmm(mod.mult.tax)

#model R^2 and eta^2 for fixed effects
r2beta(mod.mult.phylo)
r2beta(mod.mult.func)

##Threshold approach
#25% threshold

#Phylogenetic
mod.t25.phylo<-lmer(multifunct.t25~ses.pd+soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii+
                (1|site.id),data=data.scale,REML=T)
summary(mod.t25.phylo)
Anova(mod.t25.phylo,type=3)
vif.glmm(mod.t25.phylo)

#Functional
mod.t25.func<-lmer(multifunct.t25~FRaoQ+soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii+
              (1|site.id),data=data.scale,REML=T)
summary(mod.t25.func)
Anova(mod.t25.func, type=3)
vif.glmm(mod.t25.func)

#Taxonomic
mod.t25.tax<-lmer(multifunct.t25~sim.div+soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii+
             (1|site.id),data=data.scale,REML=T)
summary(mod.t25.tax)
Anova(mod.t25.tax, type=3)
vif.glmm(mod.t25.tax)

#50% threshold
#Phylogenetic
mod.t50.phylo<-lmer(multifunct.t50~ses.pd+soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii+
               (1|site.id),data=data.scale,REML=T)
summary(mod.t50.phylo)
Anova(mod.t50.phylo,type=3)
vif.glmm(mod.t50.phylo)

#Functional
mod.t50.func<-lmer(multifunct.t50~FRaoQ+soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii+
              (1|site.id),data=data.scale,REML=T)
summary(mod.t50.func)
Anova(mod.t50.func, type=3)
vif.glmm(mod.t50.func)

#Taxonomic
mod.t50.tax<-lmer(multifunct.t50~sim.div+soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii+
             (1|site.id),data=data.scale,REML=T)
summary(mod.t50.tax)
Anova(mod.t50.tax, type=3)
vif.glmm(mod.t50.tax)

#75% threshold
#Phylo
mod.t75.phylo<-lmer(multifunct.t75~ses.pd+soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii+
               (1|site.id),data=data.scale,REML=T)
summary(mod.t75.phylo)
Anova(mod.t75.phylo,type=3)
vif.glmm(mod.t75.phylo)

#Functional
mod.t75.func<-lmer(multifunct.t75~FRaoQ+soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii+
              (1|site.id),data=data.scale,REML=T)
summary(mod.t75.func)
Anova(mod.t75.func, type=3)
vif.glmm(mod.t75.func)

#Taxonomic
mod.t75.tax<-lmer(multifunct.t75~sim.div+soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii+
             (1|site.id),data=data.scale,REML=T)
summary(mod.t75.tax)
Anova(mod.t75.tax, type=3)
vif.glmm(mod.t75.tax)

####Individual ecosystem function models####

##Biomass
#Phylogenetic
mod.bio.phylo<-lmer(biomass~ses.pd+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
               (1|site.id),data=data.scale,REML=T)
summary(mod.bio.phylo)
Anova(mod.bio.phylo)

#Functional
mod.bio.func<-lmer(biomass~FRaoQ+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
              (1|site.id),data=data.scale,REML=T)
summary(mod.bio.func)
Anova(mod.bio.func)

#Taxonomic
mod.bio.tax<-lmer(biomass~sim.div+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
              (1|site.id),data=data.scale,REML=T)
summary(mod.bio.tax)
Anova(mod.bio.tax)

r2beta(mod.bio.func)

##decomposition
#Phylo
mod.decomp.phylo<-lmer(mass.loss~ses.pd+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                    (1|site.id),data=data.scale,REML=T)
summary(mod.decomp.phylo)
Anova(mod.decomp.phylo,type=3)

#Functional
mod.decomp.func<-lmer(mass.loss~FRaoQ+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                  (1|site.id),data=data.scale,REML=T)
summary(mod.decomp.func)
Anova(mod.decomp.func)

#Taxonomic
mod.decomp.tax<-lmer(mass.loss~sim.div+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                  (1|site.id),data=data.scale,REML=T)
summary(mod.decomp.tax)
Anova(mod.decomp.tax)

r2beta(mod.decomp.tax)

##Belowground biomass
#Phylo
mod.root.phylo<-lmer(roots.per.g.sand~ses.pd+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                  (1|site.id),data=data.scale,REML=T)
summary(mod.root.phylo)
Anova(mod.root.phylo)

#Functional
mod.root.func<-lmer(roots.per.g.sand~FRaoQ+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                (1|site.id),data=data.scale,REML=T)
summary(mod.root.func)
Anova(mod.root.func)

#Taxonomic
mod.root.tax<-lmer(roots.per.g.sand~sim.div+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                (1|site.id),data=data.scale,REML=T)
summary(mod.root.tax)
Anova(mod.root.tax)

r2beta(mod.root.tax)

#Seed predation
#Phylo
mod.seed.phylo<-lmer(total.seeds.removed~ses.pd+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                  (1|site.id),data=data.scale,REML=T)
summary(mod.seed.phylo)
Anova(mod.seed.phylo)

#Functional
mod.seed.func<-lmer(total.seeds.removed~FRaoQ+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                (1|site.id),data=data.scale,REML=T)
summary(mod.seed.func)
Anova(mod.seed.func)

#Taxonomic
mod.seed.tax<-lmer(total.seeds.removed~sim.div+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                (1|site.id),data=data.scale,REML=T)
summary(mod.seed.tax)
Anova(mod.seed.tax)

r2beta(mod.seed.phylo)

##Arthropod predation
#Phylo
mod.arthro.phylo<-lmer(aug.worms.gone~ses.pd+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                    (1|site.id),data=data.scale,REML=T)
summary(mod.arthro.phylo)
Anova(mod.arthro.phylo,type=3)

#Functional
mod.arthro.func<-lmer(aug.worms.gone~FRaoQ+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                  (1|site.id),data=data.scale,REML=T)
summary(mod.arthro.func)
Anova(mod.arthro.func)

#Taxonomic
mod.arthro.tax<-lmer(aug.worms.gone~sim.div+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                  (1|site.id),data=data.scale,REML=T)
summary(mod.arthro.tax)
Anova(mod.arthro.tax)

r2beta(mod.arthro.phylo)

##Floral cover
#Phylo
mod.fc.phylo<-lmer(fc.logit~ses.pd+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
              (1|site.id),data=data.scale,REML=T)
summary(mod.fc.phylo)
Anova(mod.fc.phylo,type=3)

#Functional
mod.fc.func<-lmer(fc.logit~FRaoQ+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
              (1|site.id),data=data.scale,REML=T)
summary(mod.fc.func)
Anova(mod.fc.func)

#Taxonomic
mod.fc.tax<-lmer(fc.logit~sim.div+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                (1|site.id),data=data.scale,REML=T)
summary(mod.fc.tax)
Anova(mod.fc.tax,type=3)

r2beta(mod.fc.func)

#Floral richness
#Phylo
mod.fr.phylo<-lmer(floral.plotSR~ses.pd+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                (1|site.id),data=data.scale,REML=T)
summary(mod.fr.phylo)
Anova(mod.fr.phylo)

#Functional
mod.fr.func<-lmer(floral.plotSR~FRaoQ+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
              (1|site.id),data=data.scale,REML=T)
summary(mod.fr.func)
Anova(mod.fr.func,type=3)

#Taxonomic
mod.fr.tax<-lmer(floral.plotSR~sim.div+soilPC1+age.2013+fire.return+landscape+Andropogon.gerardii+
                (1|site.id),data=data.scale,REML=T)
summary(mod.fr.tax)
Anova(mod.fr.tax)

r2beta(mod.fr.func)

####Create Table S3####
#betas, partial R^2, and full model R^2 for all models

#create a list of all models
modlist<-list(mod.mult.phylo,mod.mult.func,mod.mult.tax,
              mod.bio.phylo,mod.bio.func,mod.bio.tax,
              mod.decomp.phylo,mod.decomp.func,mod.decomp.tax,
              mod.root.phylo,mod.root.func,mod.root.tax
              ,mod.seed.phylo,mod.seed.func,mod.seed.tax
              ,mod.arthro.phylo,mod.arthro.func,mod.arthro.tax
              ,mod.fc.phylo,mod.fc.func,mod.fc.tax
              ,mod.fr.phylo,mod.fr.func,mod.fr.tax)

#extract the beta for each term in each model
effect.table<-matrix()
effect.table2<-NULL
for(i in modlist){
  mod<-rep(all.vars(formula(i))[1],6)
  func<-rep(names(fixef(i))[2],6)
  var<-names(fixef(i)[2:7])
  effect<-as.vector(fixef(i)[2:7])
  pval<-Anova(i,type=3)$Pr[2:7]
  
  effect.table<-data.frame(mod,func,var,effect,pval)
  effect.table2<-rbind(effect.table2,effect.table)
}

#extract full model and partial R^2
rsq.table<-matrix()
rsq.table2<-NULL
for(i in modlist){
  mod<-rep(all.vars(formula(i))[1],7)
  func<-rep(names(fixef(i))[2],7)
  var<-r2beta(i)[1:7,1]
  rsq<-as.vector(r2beta(i)[1:7,6])
  
  rsq.table<-data.frame(mod,func,var,rsq)
  rsq.table2<-rbind(rsq.table2,rsq.table)
}

#Calculate confidence intervals for each coef
confint.table<-matrix()
confint.table2<-NULL
for(i in modlist){
  mod<-rep(all.vars(formula(i))[1],6)
  func<-rep(names(fixef(i))[2],6)
  var<-names(fixef(i)[2:7])
  lower<-as.vector(confint(i)[4:9,1])
  upper<-as.vector(confint(i)[4:9,2])
  
  confint.table<-data.frame(mod,func,var,lower,upper)
  confint.table2<-rbind(confint.table2,confint.table)
}

#Combine the above into a single table
results.table<-merge(effect.table2,rsq.table2,by=c("mod","func","var"),all=T)
results.cast<-dcast(mod+var~func,data=setDT(results.table),value.var=c("effect","pval","rsq"),sum)

#reorder the levels for making table S3
results.cast$var<-factor(results.cast$var,
                  levels = (c("Model","ses.pd","FRaoQ","sim.div","Andropogon.gerardii","age.2013","fire.return","soilPC1","landscape")))
results.cast<-results.cast[order(results.cast$mod, results.cast$var),]
#create a .csv (uncommet code to create .csv file)
#write.csv(results.cast,"results_table2.csv")

####Create Table S2####
##AICc table

#make a list of all AICc models
modlist.aic<-list(mod.mult.phylo.aic,mod.mult.func.aic,mod.mult.tax.aic,
              mod.bio.phylo.aic,mod.bio.func.aic,mod.bio.tax.aic,
              mod.decomp.phylo.aic,mod.decomp.func.aic,mod.decomp.tax.aic,
              mod.root.phylo.aic,mod.root.func.aic,mod.root.tax.aic
              ,mod.seed.phylo.aic,mod.seed.func.aic,mod.seed.tax.aic
              ,mod.arthro.phylo.aic,mod.arthro.func.aic,mod.arthro.tax.aic
              ,mod.fc.phylo.aic,mod.fc.func.aic,mod.fc.tax.aic
              ,mod.fr.phylo.aic,mod.fr.func.aic,mod.fr.tax.aic)

#extract AICc from each model
aic.table<-matrix()
aic.table2<-NULL
for(i in modlist.aic){
  mod<-all.vars(formula(i))[1]
  div<-names(fixef(i))[2]
  AICc<-AICc(i)
  
  aic.table<-data.frame(mod,div,AICc)
  aic.table2<-rbind(aic.table2,aic.table)
}

#output AICc table .csv file (uncomment code to create .csv file)
#write.csv(aic.table2,"aic.table.csv")

####Create Figure 2####

#Create models with term of interest removed to plot against residuals 
mod.no.func<-lm(multifunct.a~soilPC1+landscape+age.2013+fire.return+Andropogon.gerardii,data=full.data)
mod.no.lan<-lm(multifunct.a~FRaoQ+soilPC1+age.2013+fire.return+Andropogon.gerardii,data=full.data)
mod.no.andger<-lm(multifunct.a~FRaoQ+soilPC1+age.2013+fire.return+landscape,data=full.data)

#Functional diversity
p1<-ggplot(full.data, aes(x = FRaoQ, y = resid(mod.no.func)))+
  geom_point(size=1.5)+
  geom_smooth(method="lm",size=2,color="black",se=F)+
  annotate("text",label="p=0.002",x=0.005,y=.15,size=3)+
  annotate("text",label="paste(R ^ 2, \" = 0.05\")",x=0.005,y=.1225,parse = TRUE,size=3)+
  labs(x = "Functional diversity", y = "Ecosystem\nmultifunctionality")+
  theme(text = element_text(size=24),axis.text=element_text(colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(size=.7, color="black"),legend.position="none")

#Landscape context
p2<-ggplot(full.data, aes(x = landscape, y = resid(mod.no.lan)))+
  geom_point(size=1.5)+
  geom_smooth(method="lm",size=2,color="black",se=F)+
  annotate("text",label="p>0.0001",x=-4.75e+05,y=.17,size=3)+
  annotate("text",label="paste(R ^ 2, \" = 0.19\")",x=-4.75e+05,y=.1425,parse = TRUE,size=3)+
  labs(x = "Landscape Context", y = "Ecosystem\nmultifunctionality")+
  theme(text = element_text(size=24),axis.text=element_text(colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(size=.7, color="black"),legend.position="none")

#A. gerardii abundance
p3<-ggplot(full.data, aes(x = Andropogon.gerardii, y = resid(mod.no.andger)))+
  geom_point(size=1.5)+
  geom_smooth(method="lm",size=2,color="black",se=F)+
  annotate("text",label="p=0.007",x=12,y=.17,size=3)+
  annotate("text",label="paste(R ^ 2, \" = 0.04\")",x=12,y=.1425,parse = TRUE,size=3)+
  labs(x = "A. gerardii abundance", y = "Ecosystem\nmultifunctionality")+
  theme(text = element_text(size=24),axis.text=element_text(colour="black"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line = element_line(size=.7, color="black"),legend.position="none")

#combine ggplot outputs into a single pannel
multiplot <- function(..., plotlist=NULL, cols) {
  require(grid)
  
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  
  plotCols = cols                         
  plotRows = ceiling(numPlots/plotCols) 
  
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  } 
}

#output figure 2 to PDF (uncomment code)
#pdf("figure_2.pdf",height=11,width=6)
multiplot(p1,p3,p2, cols=1)
dev.off()

####Create Figure 1#####

#list of best models based on AICc
aic.modlist<-list(mod.mult.func,
              mod.bio.phylo,
              mod.decomp.phylo,
              mod.root.phylo,mod.root.func,mod.root.tax,
              mod.seed.phylo,mod.seed.func,mod.seed.tax,
              mod.arthro.phylo,mod.arthro.func,mod.arthro.tax
              ,mod.fc.func
              ,mod.fr.func)

#extract betas
effect.table.aic<-matrix()
effect.table2.aic<-NULL
for(i in aic.modlist){
  mod<-rep(all.vars(formula(i))[1],6)
  func<-rep(names(fixef(i))[2],6)
  var<-names(fixef(i)[2:7])
  effect<-as.vector(fixef(i)[2:7])
  pval<-Anova(i,type=3)$Pr[2:7]
  
  effect.table.aic<-data.frame(mod,func,var,effect,pval)
  effect.table2.aic<-rbind(effect.table2.aic,effect.table.aic)
}

#extract R2
rsq.table.aic<-matrix()
rsq.table2.aic<-NULL
for(i in aic.modlist){
  mod<-rep(all.vars(formula(i))[1],7)
  func<-rep(names(fixef(i))[2],7)
  var<-r2beta(i)[1:7,1]
  rsq<-as.vector(r2beta(i)[1:7,6])
  
  rsq.table.aic<-data.frame(mod,func,var,rsq)
  rsq.table2.aic<-rbind(rsq.table2.aic,rsq.table.aic)
}

#calculate confidence intervals
confint.table.aic<-matrix()
confint.table2.aic<-NULL
for(i in aic.modlist){
  mod<-rep(all.vars(formula(i))[1],6)
  func<-rep(names(fixef(i))[2],6)
  var<-names(fixef(i)[2:7])
  lower<-as.vector(confint(i)[4:9,1])
  upper<-as.vector(confint(i)[4:9,2])
  
  confint.table.aic<-data.frame(mod,func,var,lower,upper)
  confint.table2.aic<-rbind(confint.table2.aic,confint.table.aic)
}

#create combined data.frame
results.table.aic<-merge(effect.table2.aic,rsq.table2.aic,by=c("mod","func","var"),all=T)
results.table.aic<-merge(results.table.aic,confint.table2.aic,by=c("mod","func","var"),all=T)

results.table.aic<-results.table.aic[results.table.aic$var!="Model",]
results.table.aic<-aggregate(.~mod+var,data=results.table.aic,FUN=mean)

#Reorder and relabel variables
results.table.aic$var <- factor(results.table.aic$var, levels = c("Model","landscape","soilPC1","fire.return","age.2013","Andropogon.gerardii","sim.div","FRaoQ","ses.pd"))
labels2<- c(multifunct.a="Ecosystem multifunctionality",biomass="Aboveground biomass", mass.loss="Decomposition rate", roots.per.g.sand="Belowground biomass",total.seeds.removed="Seed predation",aug.worms.gone="Arthropod predation",fc.logit="Pollinator resource availability",floral.plotSR="Floral richness")
var.label<-c("landscape"="Landscape context","soilPC1"="Soil moisture","fire.return"="Fire frequency","age.2013"="Site age","Andropogon.gerardii"="A. gerardii\nabundance","sim.div"="Taxonomic\ndiversity","FRaoQ"="Functional\ndiversity","ses.pd"="Phylogenetic\ndiversity")

#create the plot (uncoment code to export PDF file)
#pdf("effects.pdf",width=8.5,height=11)
ggplot(results.table.aic, aes(var,effect))+
  geom_point(size=1.5)+
  facet_wrap(~mod,ncol=2,labeller=labeller(mod = labels2),scales="free_y")+
  geom_hline(yintercept=0, linetype="dotted")+
  coord_flip()+
  geom_errorbar(width=0, aes(ymin=lower, ymax=upper)) +
  xlab(label="")+
  ylab(label="Standardized regression coefficents")+
  scale_x_discrete(labels = var.label)+
  theme(text = element_text(size=14),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")
dev.off()

####Create Figure S1####

#list of best threshold models
thresh.modlist<-list(mod.mult.func,
                  mod.t25.func,
                  mod.t50.phylo,mod.t50.func,mod.t50.tax,
                  mod.t75.phylo)

#extract betas
effect.table.thresh<-matrix()
effect.table2.thresh<-NULL
for(i in thresh.modlist){
  mod<-rep(all.vars(formula(i))[1],6)
  func<-rep(names(fixef(i))[2],6)
  var<-names(fixef(i)[2:7])
  effect<-as.vector(fixef(i)[2:7])
  pval<-Anova(i,type=3)$Pr[2:7]
  
  effect.table.thresh<-data.frame(mod,func,var,effect,pval)
  effect.table2.thresh<-rbind(effect.table2.thresh,effect.table.thresh)
}

#extract R2
rsq.table.thresh<-matrix()
rsq.table2.thresh<-NULL
for(i in thresh.modlist){
  mod<-rep(all.vars(formula(i))[1],7)
  func<-rep(names(fixef(i))[2],7)
  var<-r2beta(i)[1:7,1]
  rsq<-as.vector(r2beta(i)[1:7,6])
  
  rsq.table.thresh<-data.frame(mod,func,var,rsq)
  rsq.table2.thresh<-rbind(rsq.table2.thresh,rsq.table.thresh)
}

#calculate confidence intervals
confint.table.thresh<-matrix()
confint.table2.thresh<-NULL
for(i in thresh.modlist){
  mod<-rep(all.vars(formula(i))[1],6)
  func<-rep(names(fixef(i))[2],6)
  var<-names(fixef(i)[2:7])
  lower<-as.vector(confint(i)[4:9,1])
  upper<-as.vector(confint(i)[4:9,2])
  
  confint.table.thresh<-data.frame(mod,func,var,lower,upper)
  confint.table2.thresh<-rbind(confint.table2.thresh,confint.table.thresh)
}

#combine into data.frame
results.table.thresh<-merge(effect.table2.thresh,rsq.table2.thresh,by=c("mod","func","var"),all=T)
results.table.thresh<-merge(results.table.thresh,confint.table2.thresh,by=c("mod","func","var"),all=T)

results.table.thresh<-results.table.thresh[results.table.thresh$var!="Model",]
results.table.thresh<-aggregate(.~mod+var,data=results.table.thresh,FUN=mean)

#reorder and rename variables
results.table.thresh$var <- factor(results.table.thresh$var, levels = c("landscape","soilPC1","fire.return","age.2013","Andropogon.gerardii","sim.div","FRaoQ","ses.pd"))
labels2<- c(multifunct.a="Ecosystem multifunctionality\n(Averaging)",multifunct.t25="Ecosystem multifunctionality\n(25% Threshold)",multifunct.t50="Ecosystem multifunctionality\n(50% Threshold)",multifunct.t75="Ecosystem multifunctionality\n(75% Threshold)")
var.label<-c("landscape"="Landscape context","soilPC1"="Soil moisture","fire.return"="Fire frequency","age.2013"="Site age","Andropogon.gerardii"="A. gerardii abundance","sim.div"="Taxonomic diversity","FRaoQ"="Functional diversity","ses.pd"="Phylogenetic diversity")

#create plot (uncomment to export PDF file)
#pdf("effects_thresh.pdf",width=8.5,height=6)
ggplot(results.table.thresh, aes(var,effect))+
  geom_point(size=1.5)+
  facet_wrap(~mod,ncol=2,labeller=labeller(mod = labels2),scales="free_y")+
  geom_hline(yintercept=0, linetype="dotted")+
  coord_flip()+
  geom_errorbar(width=0, aes(ymin=lower, ymax=upper)) +
  xlab(label="")+
  ylab(label="Standardized regression coefficents")+
  scale_x_discrete(labels = var.label)+
  theme(text = element_text(size=14),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")
dev.off()

####Create Figure S4####

#This function computes both the loess and line of best fit for bivariate plots
#commenting out loess fit line for this figure
panel.fit<-function(x, y, ...) {
  ll <- loess(y~x)
  points(x,y, ...)
  nx<-seq(min(x), max(x), length.out=150)
  #lines(nx, predict(ll, nx), col="blue",lwd=2)
  abline(a = lm(y ~ x)$coefficients[1] , b = lm(y ~ x)$coefficients[2] ,col="black",lwd=2, ...)
}

#This function computes the correlation between every combination of variables
panel.cor <- function(x, y, digits=2, prefix="",cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 1/strwidth(txt)
  text(0.5, 0.5, txt, cex = 3)
}

#This function computes the density histogram and boxplot of each variable
panel.density = function(x,...)
{
  pu <- par("usr")
  d <- density(x,...)
  par("usr" = c(pu[1:2], 0, max(d$y)*1.5))
  lines(d,lwd=2)
  par("usr" = c(pu[1:2], 0, 1))
  boxplot(x, at=0.5, boxwex=0.3, horizontal=TRUE, add=TRUE, axes=F)
}

#create data frame with diversity metrics
div.cor.data<-full.data[c(4:6)]
names(div.cor.data)<-c("Simpson's diversity","Rao's Q","SES Phylogenetic diversity")

#plot (uncomment to export PDF file)
#pdf("div_cor.pdf",width=8.5,height=6)
pairs(div.cor.data,na.action = stats::na.omit,
      upper.panel=panel.cor,
      lower.panel=panel.fit,
      diag.panel=panel.density)
dev.off()

#####Create figures S2 & S3####

#Create a data frame of the ecosystem functions
eco.funct<-dcast(site.id+plot~variable,data=relevant.fun)
eco.funct$roots.per.g.sand<-eco.funct$roots.per.g.sand*1000
funct.scale<-scale(eco.funct[,3:9])
funct.scale<-cbind(eco.funct[,1:2],funct.scale)

#Create models and plots (both linear and non-linear) for each pairwise EF combination
mod.bio.root<-lmer(roots.per.g.sand~biomass+(1|site.id),data=eco.funct,REML=T)
mod.bio.root.gam<-gamm(roots.per.g.sand~biomass,random=list(site.id=~1),data=funct.scale)
mod.bio.root.gam.s<-gamm(roots.per.g.sand~s(biomass),random=list(site.id=~1),data=funct.scale)
AICc(mod.bio.root.gam,mod.bio.root.gam.s)

c1<-ggplot(eco.funct,aes(biomass,roots.per.g.sand))+
 geom_point(size=.75)+
 stat_function(fun=function(x)(fixef(mod.bio.root)[1] + fixef(mod.bio.root)[2]*x),size=1,color="black")+
  annotate("text",label="r=0.06",x=250,y=0.015,size=3)+
  labs(x = "Aboveground biomass (g)", y="Belowgorund biomass\n(mg roots/g sand)")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

d1<-ggplot(eco.funct,aes(biomass,roots.per.g.sand))+
  geom_point(size=.75)+
  geom_smooth(method=gam,formula=y~s(x),se=F,color="black")+
  annotate("text",label="paste(R ^ 2, \" = 0.01\")",parse=T,x=250,y=0.015,size=3)+
  labs(x = "Aboveground biomass (g)", y="Belowgorund biomass\n(mg roots/g sand)")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

summary(lmer(roots.per.g.sand~biomass+(1|site.id),data=funct.scale))

mod.bio.decomp<-lmer(mass.loss~biomass+(1|site.id),data=funct.scale,REML=F)
mod.bio.decomp.gam<-gamm(mass.loss~biomass,random=list(site.id=~1),data=funct.scale)
mod.bio.decomp.gam.s<-gamm(mass.loss~s(biomass),random=list(site.id=~1),data=funct.scale)
AICc(mod.bio.decomp.gam,mod.bio.decomp.gam.s)

c2<-ggplot(eco.funct,aes(biomass,mass.loss))+
  geom_point(size=.75)+
  stat_function(fun=function(x)(fixef(mod.bio.decomp)[1] + fixef(mod.bio.decomp)[2]*x),size=1,color="black")+
  annotate("text",label="r=0.12*",x=250,y=0.03,size=3)+
  labs(x = "Aboveground biomass (g)", y="Decomposition rate\n(g/day)")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

d2<-ggplot(eco.funct,aes(biomass,mass.loss))+
  geom_point(size=.75)+
  geom_smooth(method="gam",formula=y~s(x),se=F,color="black")+
  annotate("text",label="paste(R ^ 2, \" = 0.01*\")",parse=T,x=250,y=0.03,size=3)+
  labs(x = "Aboveground biomass (g)", y="Decomposition rate\n(g/day)")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

summary(lmer(mass.loss~biomass+(1|site.id),data=funct.scale))

mod.bio.seed<-lmer(total.seeds.removed~biomass+(1|site.id),data=funct.scale,REML=F)
mod.bio.seed.gam<-gamm(total.seeds.removed~biomass,random=list(site.id=~1),data=funct.scale)
mod.bio.seed.gam.s<-gamm(total.seeds.removed~s(biomass),random=list(site.id=~1),data=funct.scale)
AICc(mod.bio.seed.gam,mod.bio.seed.gam.s)

c3<-ggplot(eco.funct,aes(biomass,total.seeds.removed))+
  geom_point(size=.75)+
  stat_function(fun=function(x)(fixef(mod.bio.seed)[1] + fixef(mod.bio.seed)[2]*x),size=1,color="black")+
  annotate("text",label="r=-0.02",x=250,y=70,size=3)+
  labs(x = "Aboveground biomass (g)", y="Seed predation rate\n(seeds/day)")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

d3<-ggplot(eco.funct,aes(biomass,total.seeds.removed))+
  geom_point(size=0.75)+
  geom_smooth(method="gam",formula=y~s(x),se=F,color="black")+
  annotate("text",label="paste(R ^ 2, \" = 0.003\")",parse=T,x=250,y=70,size=3)+
  labs(x = "Aboveground biomass (g)", y="Seed predation rate\n(seeds/day)")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

summary(lmer(total.seeds.removed~biomass+(1|site.id),data=funct.scale))

mod.bio.arthro<-lmer(aug.worms.gone~biomass+(1|site.id),data=funct.scale,REML=F)
mod.bio.arthro.gam<-gamm(aug.worms.gone~biomass,random=list(site.id=~1),data=funct.scale)
mod.bio.arthro.gam.s<-gamm(aug.worms.gone~s(biomass),random=list(site.id=~1),data=funct.scale)
AICc(mod.bio.arthro.gam,mod.bio.arthro.gam.s)

c4<-ggplot(eco.funct,aes(biomass,aug.worms.gone))+
  geom_point(size=.75)+
  stat_function(fun=function(x)(fixef(mod.bio.arthro)[1] + fixef(mod.bio.arthro)[2]*x),size=1,color="black")+
  annotate("text",label="r=-0.05",x=250,y=4.5,size=3)+
  labs(x = "Aboveground biomass (g)", y="Arthropod predation rate\n(g/day)")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

d4<-ggplot(eco.funct,aes(biomass,aug.worms.gone))+
  geom_point(size=0.75)+
  geom_smooth(method="gam",formula=y~s(x),se=F,color="black")+
  annotate("text",label="paste(R ^ 2, \" = 0.001\")",parse=T,x=250,y=3.5,size=3)+
  labs(x = "Aboveground biomass (g)", y="Arthropod predation rate\n(g/day)")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

summary(lmer(aug.worms.gone~biomass+(1|site.id),data=funct.scale))

mod.bio.fc<-lmer(floral.cover~biomass+(1|site.id),data=funct.scale,REML=F)
mod.bio.fc.gam<-gamm(floral.cover~biomass,random=list(site.id=~1),data=funct.scale)
mod.bio.fc.gam.s<-gamm(floral.cover~s(biomass),random=list(site.id=~1),data=funct.scale)
AICc(mod.bio.fc.gam,mod.bio.fc.gam.s)

c5<-ggplot(eco.funct,aes(biomass,floral.cover))+
  geom_point(size=.75)+
  stat_function(fun=function(x)(fixef(mod.bio.fc)[1] + fixef(mod.bio.fc)[2]*x),size=1,color="black")+
  annotate("text",label="r=-0.07",x=325,y=45,size=3)+
  labs(x = "Aboveground biomass (g)", y="Pollinator resource\navailability")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

d5<-ggplot(eco.funct,aes(biomass,floral.cover))+
  geom_point(size=0.75)+
  geom_smooth(method="gam",formula=y~s(x),se=F,color="black")+
  annotate("text",label="paste(R ^ 2, \" = 0.10*\")",parse=T,x=325,y=45,size=3)+
  labs(x = "Aboveground biomass (g)", y="Pollinator resource\navailability")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

summary(lmer(floral.cover~biomass+(1|site.id),data=funct.scale))

mod.bio.fr<-lmer(floral.plotSR~biomass+(1|site.id),data=funct.scale,REML=F)
mod.bio.fr.gam<-gamm(floral.plotSR~biomass,random=list(site.id=~1),data=funct.scale)
mod.bio.fr.gam.s<-gamm(floral.plotSR~s(biomass),random=list(site.id=~1),data=funct.scale)
AICc(mod.bio.fr.gam,mod.bio.fr.gam.s)

c6<-ggplot(eco.funct,aes(biomass,floral.plotSR))+
  geom_point(size=.75)+
  stat_function(fun=function(x)(fixef(mod.bio.fr)[1] + fixef(mod.bio.fr)[2]*x),size=1,color="black")+
  annotate("text",label="r=-0.03",x=250,y=3,size=3)+
  labs(x = "Aboveground biomass (g)", y="Floral richness")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

d6<-ggplot(eco.funct,aes(biomass,floral.plotSR))+
  geom_point(size=0.75)+
  geom_smooth(method="gam",formula=y~s(x),se=F,color="black")+
  annotate("text",label="paste(R ^ 2, \" = 0.002\")",parse=T,x=250,y=3,size=3)+
  labs(x = "Aboveground biomass (g)", y="Floral richness")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

summary(lmer(floral.plotSR~biomass+(1|site.id),data=funct.scale))

mod.decomp.roots<-lmer(roots.per.g.sand~mass.loss+(1|site.id),data=funct.scale,REML=F)
mod.decomp.roots.gam<-gamm(roots.per.g.sand~mass.loss,random=list(site.id=~1),data=funct.scale)
mod.decomp.roots.gam.s<-gamm(roots.per.g.sand~s(mass.loss),random=list(site.id=~1),data=funct.scale)
AICc(mod.decomp.roots.gam,mod.decomp.roots.gam.s)

c7<-ggplot(eco.funct,aes(mass.loss,roots.per.g.sand))+
  geom_point(size=.75)+
  stat_function(fun=function(x)(fixef(mod.decomp.roots)[1] + fixef(mod.decomp.roots)[2]*x),size=1,color="black")+
  annotate("text",label="r=-0.09",x=0.01,y=0.015,size=3)+
  labs(x = "Decomposition rate (g/day)", y="Belowground biomass\n(mg roots/g sand)")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

d7<-ggplot(eco.funct,aes(mass.loss,roots.per.g.sand))+
  geom_point(size=0.75)+
  geom_smooth(method="gam",formula=y~s(x),se=F,color="black")+
  annotate("text",label="paste(R ^ 2, \" = 0.02\")",parse=T,x=0.01,y=0.015,size=3)+
  labs(x = "Decomposition rate (g/day)", y="Belowground biomass\n(mg roots/g sand)")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

summary(lmer(roots.per.g.sand~mass.loss+(1|site.id),data=funct.scale))

mod.decomp.seed<-lmer(total.seeds.removed~mass.loss+(1|site.id),data=funct.scale,REML=F)
mod.decomp.seed.gam<-gamm(total.seeds.removed~mass.loss,random=list(site.id=~1),data=funct.scale)
mod.decomp.seed.gam.s<-gamm(total.seeds.removed~s(mass.loss),random=list(site.id=~1),data=funct.scale)
AICc(mod.decomp.seed.gam,mod.decomp.seed.gam.s)

c8<-ggplot(eco.funct,aes(mass.loss,total.seeds.removed))+
  geom_point(size=.75)+
  stat_function(fun=function(x)(fixef(mod.decomp.seed)[1] + fixef(mod.decomp.seed)[2]*x),size=1,color="black")+
  annotate("text",label="r=0.01",x=0.01,y=70,size=3)+
  labs(x = "Decomposition rate (g/day)", y="Seed predation rate\n(seeds/day)")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

d8<-ggplot(eco.funct,aes(mass.loss,total.seeds.removed))+
  geom_point(size=0.75)+
  geom_smooth(method="gam",formula=y~s(x),se=F,color="black")+
  annotate("text",label="paste(R ^ 2, \" = 0.004\")",parse=T,x=0.01,y=70,size=3)+
  labs(x = "Decomposition rate (g/day)", y="Seed predation rate\n(seeds/day)")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

summary(lmer(total.seeds.removed~mass.loss+(1|site.id),data=funct.scale))

mod.decomp.arthro<-lmer(aug.worms.gone~mass.loss+(1|site.id),data=funct.scale,REML=F)
mod.decomp.arthro.gam<-gamm(aug.worms.gone~mass.loss,random=list(site.id=~1),data=funct.scale)
mod.decomp.arthro.gam.s<-gamm(aug.worms.gone~s(mass.loss),random=list(site.id=~1),data=funct.scale)
AICc(mod.decomp.arthro.gam,mod.decomp.arthro.gam.s)

c9<-ggplot(eco.funct,aes(mass.loss,aug.worms.gone))+
  geom_point(size=.75)+
  stat_function(fun=function(x)(fixef(mod.decomp.arthro)[1] + fixef(mod.decomp.arthro)[2]*x),size=1,color="black")+
  annotate("text",label="r=0.03",x=0.01,y=3.5,size=3)+
  labs(x = "Decomposition rate (g/day)", y="Arthopod Predation rate\n(worms/day)")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

d9<-ggplot(eco.funct,aes(mass.loss,aug.worms.gone))+
  geom_point(size=.75)+
  geom_smooth(method="gam",formula=y~s(x),se=F,color="black")+
  annotate("text",label="paste(R ^ 2, \" = 0.003\")",parse=T,x=0.01,y=3.5,size=3)+
  labs(x = "Decomposition rate (g/day)", y="Arthopod Predation rate\n(worms/day)")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

summary(lmer(aug.worms.gone~mass.loss+(1|site.id),data=funct.scale))

mod.decomp.fc<-lmer(floral.cover~mass.loss+(1|site.id),data=funct.scale,REML=F)
mod.decomp.fc.gam<-gamm(floral.cover~mass.loss,random=list(site.id=~1),data=funct.scale)
mod.decomp.fc.gam.s<-gamm(floral.cover~s(mass.loss),random=list(site.id=~1),data=funct.scale)
AICc(mod.decomp.fc.gam,mod.decomp.fc.gam.s)

c10<-ggplot(eco.funct,aes(mass.loss,floral.cover))+
  geom_point(size=.75)+
  stat_function(fun=function(x)(fixef(mod.decomp.fc)[1] + fixef(mod.decomp.fc)[2]*x),size=1,color="black")+
  annotate("text",label="r=0.02",x=0.01,y=45,size=3)+
  labs(x = "Decomposition rate (g/day)", y="Pollinator resource\navailability")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

d10<-ggplot(eco.funct,aes(mass.loss,floral.cover))+
  geom_point(size=.75)+
  geom_smooth(method="gam",formula=y~s(x),se=F,color="black")+
  annotate("text",label="paste(R ^ 2, \" = 0.006\")",parse=T,x=0.01,y=45,size=3)+
  labs(x = "Decomposition rate (g/day)", y="Pollinator resource\navailability")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

summary(lmer(floral.cover~mass.loss+(1|site.id),data=funct.scale))

mod.decomp.fr<-lmer(floral.plotSR~mass.loss+(1|site.id),data=eco.funct)
mod.decomp.fr.gam<-gamm(floral.plotSR~mass.loss,random=list(site.id=~1),data=funct.scale)
mod.decomp.fr.gam.s<-gamm(floral.plotSR~s(mass.loss),random=list(site.id=~1),data=funct.scale)
AICc(mod.decomp.fr.gam,mod.decomp.fr.gam.s)

c11<-ggplot(eco.funct,aes(mass.loss,floral.plotSR))+
  geom_point(size=.75)+
  stat_function(fun=function(x)(fixef(mod.decomp.fr)[1] + fixef(mod.decomp.fr)[2]*x),size=1,color="black")+
  annotate("text",label="r=-0.13*",x=0.002,y=3.5,size=3)+
  labs(x = "Decomposition rate (g/day)", y="Floral richness")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

d11<-ggplot(eco.funct,aes(mass.loss,floral.plotSR))+
  geom_point(size=.75)+
  geom_smooth(method="gam",formula=y~s(x),se=F,color="black")+
  annotate("text",label="paste(R ^ 2, \" = 0.04*\")",parse=T,x=0.025,y=3,size=3)+
  labs(x = "Decomposition rate (g/day)", y="Floral richness")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

summary(lmer(floral.plotSR~mass.loss+(1|site.id),data=funct.scale))

mod.roots.seed<-lmer(total.seeds.removed~roots.per.g.sand+(1|site.id),data=eco.funct)
mod.roots.seed.gam<-gamm(total.seeds.removed~roots.per.g.sand,random=list(site.id=~1),data=funct.scale)
mod.roots.seed.gam.s<-gamm(total.seeds.removed~s(roots.per.g.sand),random=list(site.id=~1),data=funct.scale)
AICc(mod.roots.seed.gam,mod.roots.seed.gam.s)

c12<-ggplot(eco.funct,aes(roots.per.g.sand,total.seeds.removed))+
  geom_point(size=.75)+
  stat_function(fun=function(x)(fixef(mod.roots.seed)[1] + fixef(mod.roots.seed)[2]*x),size=1,color="black")+
  annotate("text",label="r=0.01",x=0.001,y=70,size=3)+
  labs(x = "Belowground biomass\n(mg roots/g sand)", y="Seed predation rate\n(seeds/day)")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

d12<-ggplot(eco.funct,aes(roots.per.g.sand,total.seeds.removed))+
  geom_point(size=.75)+
  geom_smooth(method="gam",formula=y~s(x),se=F,color="black")+
  annotate("text",label="paste(R ^ 2, \" = 0.001\")",parse=T,x=0.009,y=70,size=3)+
  labs(x = "Belowground biomass\n(mg roots/g sand)", y="Seed predation rate\n(seeds/day)")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

summary(lmer(total.seeds.removed~roots.per.g.sand+(1|site.id),data=funct.scale))

mod.roots.arthro<-lmer(aug.worms.gone~roots.per.g.sand+(1|site.id),data=eco.funct)
mod.roots.arthro.gam<-gamm(aug.worms.gone~roots.per.g.sand,random=list(site.id=~1),data=funct.scale)
mod.roots.arthro.gam.s<-gamm(aug.worms.gone~s(roots.per.g.sand),random=list(site.id=~1),data=funct.scale)
AICc(mod.roots.arthro.gam,mod.roots.arthro.gam.s)

c13<-ggplot(eco.funct,aes(roots.per.g.sand,aug.worms.gone))+
  geom_point(size=.75)+
  stat_function(fun=function(x)(fixef(mod.roots.arthro)[1] + fixef(mod.roots.arthro)[2]*x),size=1,color="black")+
  annotate("text",label="r=-0.09",x=0.001,y=3.5,size=3)+
  labs(x = "Belowground biomass\n(mg roots/g sand)", y="Arthropod predation rate\n(seeds/day)")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

d13<-ggplot(eco.funct,aes(roots.per.g.sand,aug.worms.gone))+
  geom_point(size=.75)+
  geom_smooth(method="gam",formula=y~s(x),se=F,color="black")+
  annotate("text",label="paste(R ^ 2, \" = 0.006\")",parse=T,x=0.0025,y=3.5,size=3)+
  labs(x = "Belowground biomass\n(mg roots/g sand)", y="Arthropod predation rate\n(seeds/day)")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

summary(lmer(aug.worms.gone~roots.per.g.sand+(1|site.id),data=funct.scale))

mod.roots.fc<-lmer(floral.cover~roots.per.g.sand+(1|site.id),data=eco.funct)
mod.roots.fc.gam<-gamm(floral.cover~roots.per.g.sand,random=list(site.id=~1),data=funct.scale)
mod.roots.fc.gam.s<-gamm(floral.cover~s(roots.per.g.sand),random=list(site.id=~1),data=funct.scale)
AICc(mod.roots.fc.gam,mod.roots.fc.gam.s)

c14<-ggplot(eco.funct,aes(roots.per.g.sand,floral.cover))+
  geom_point(size=.75)+
  stat_function(fun=function(x)(fixef(mod.roots.fc)[1] + fixef(mod.roots.fc)[2]*x),size=1,color="black")+
  annotate("text",label="r=-0.05",x=0.001,y=45,size=3)+
  labs(x = "Belowground biomass\n(mg roots/g sand)", y="Pollinator resource\navailability")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

d14<-ggplot(eco.funct,aes(roots.per.g.sand,floral.cover))+
  geom_point(size=.75)+
  geom_smooth(method="gam",formula=y~s(x),se=F,color="black")+
  annotate("text",label="paste(R ^ 2, \" = 0.001\")",parse=T,x=0.005,y=45,size=3)+
  labs(x = "Belowground biomass\n(mg roots/g sand)", y="Pollinator resource\navailability")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

summary(lmer(floral.cover~roots.per.g.sand+(1|site.id),data=funct.scale))

mod.roots.fr<-lmer(floral.plotSR~roots.per.g.sand+(1|site.id),data=eco.funct)
mod.roots.fr.gam<-gamm(floral.plotSR~roots.per.g.sand,random=list(site.id=~1),data=funct.scale)
mod.roots.fr.gam.s<-gamm(floral.plotSR~s(roots.per.g.sand),random=list(site.id=~1),data=funct.scale)
AICc(mod.roots.fr.gam,mod.roots.fr.gam.s)

c15<-ggplot(eco.funct,aes(roots.per.g.sand,floral.plotSR))+
  geom_point(size=.75)+
  stat_function(fun=function(x)(fixef(mod.roots.fr)[1] + fixef(mod.roots.fr)[2]*x),size=1,color="black")+
  annotate("text",label="r=-0.0008",x=0.001,y=3.5,size=3)+
  labs(x = "Belowground biomass\n(mg roots/g sand)", y="Floral richness")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

d15<-ggplot(eco.funct,aes(roots.per.g.sand,floral.plotSR))+
  geom_point(size=.75)+
  geom_smooth(method="gam",formula=y~s(x),se=F,color="black")+
  annotate("text",label="paste(R ^ 2, \" = 0.003\")",parse=T,x=0.0075,y=3,size=3)+
  labs(x = "Belowground biomass\n(mg roots/g sand)", y="Floral richness")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

summary(lmer(floral.plotSR~roots.per.g.sand+(1|site.id),data=funct.scale))

mod.seed.arthro<-lmer(aug.worms.gone~total.seeds.removed+(1|site.id),data=eco.funct)
mod.seed.arthro.gam<-gamm(aug.worms.gone~total.seeds.removed,random=list(site.id=~1),data=funct.scale)
mod.seed.arthro.gam.s<-gamm(aug.worms.gone~s(total.seeds.removed),random=list(site.id=~1),data=funct.scale)
AICc(mod.seed.arthro.gam,mod.seed.arthro.gam.s)

c16<-ggplot(eco.funct,aes(total.seeds.removed,aug.worms.gone))+
  geom_point(size=.75)+
  stat_function(fun=function(x)(fixef(mod.seed.arthro)[1] + fixef(mod.seed.arthro)[2]*x),size=1,color="black")+
  annotate("text",label="r=0.03",x=10,y=3.5,size=3)+
  labs(x = "Seed predation rate\n(seeds/day)", y="Arthropod predation rate\n(worms/day)")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

d16<-ggplot(eco.funct,aes(total.seeds.removed,aug.worms.gone))+
  geom_point(size=.75)+
  geom_smooth(method="gam",formula=y~s(x),se=F,color="black")+
  annotate("text",label="paste(R ^ 2, \" = 0.01\")",parse=T,x=11,y=3.5,size=3)+
  labs(x = "Seed predation rate\n(seeds/day)", y="Arthropod predation rate\n(worms/day)")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

summary(lmer(aug.worms.gone~total.seeds.removed+(1|site.id),data=funct.scale))

mod.seed.fc<-lmer(floral.cover~total.seeds.removed+(1|site.id),data=eco.funct)
mod.seed.fc.gam<-gamm(floral.cover~total.seeds.removed,random=list(site.id=~1),data=funct.scale)
mod.seed.fc.gam.s<-gamm(floral.cover~s(total.seeds.removed),random=list(site.id=~1),data=funct.scale)
AICc(mod.seed.fc.gam,mod.seed.fc.gam.s)

c17<-ggplot(eco.funct,aes(total.seeds.removed,floral.cover))+
  geom_point(size=.75)+
  stat_function(fun=function(x)(fixef(mod.seed.fc)[1] + fixef(mod.seed.fc)[2]*x),size=1,color="black")+
  annotate("text",label="r=0.14*",x=11,y=45,size=3)+
  labs(x = "Seed predation rate\n(seeds/day)", y="Pollinator resource\navailability")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

d17<-ggplot(eco.funct,aes(total.seeds.removed,floral.cover))+
  geom_point(size=.75)+
  geom_smooth(method="gam",formula=y~s(x),se=F,color="black")+
  annotate("text",label="paste(R ^ 2, \" = 0.01*\")",parse=T,x=12,y=45,size=3)+
  labs(x = "Seed predation rate\n(seeds/day)", y="Pollinator resource\navailability")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

summary(lmer(floral.cover~total.seeds.removed+(1|site.id),data=funct.scale))

mod.seed.fr<-lmer(floral.plotSR~total.seeds.removed+(1|site.id),data=eco.funct)
mod.seed.fr.gam<-gamm(floral.plotSR~total.seeds.removed,random=list(site.id=~1),data=funct.scale)
mod.seed.fr.gam.s<-gamm(floral.plotSR~s(total.seeds.removed),random=list(site.id=~1),data=funct.scale)
AICc(mod.seed.fr.gam,mod.seed.fr.gam.s)

c18<-ggplot(eco.funct,aes(total.seeds.removed,floral.plotSR))+
  geom_point(size=.75)+
  stat_function(fun=function(x)(fixef(mod.seed.fr)[1] + fixef(mod.seed.fr)[2]*x),size=1,color="black")+
  annotate("text",label="r=0.06",x=10,y=3.5,size=3)+
  labs(x = "Seed predation rate\n(seeds/day)", y="Floral richness")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

d18<-ggplot(eco.funct,aes(total.seeds.removed,floral.plotSR))+
  geom_point(size=.75)+
  geom_smooth(method="gam",formula=y~s(x),se=F,color="black")+
  annotate("text",label="paste(R ^ 2, \" = 0.002\")",parse=T,x=65,y=3,size=3)+
  labs(x = "Seed predation rate\n(seeds/day)", y="Floral richness")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

summary(lmer(floral.plotSR~total.seeds.removed+(1|site.id),data=funct.scale))

mod.fc.arthro<-lmer(aug.worms.gone~floral.cover+(1|site.id),data=eco.funct)
mod.fc.arthro.gam<-gamm(aug.worms.gone~floral.cover,random=list(site.id=~1),data=funct.scale)
mod.fc.arthro.gam.s<-gamm(aug.worms.gone~s(floral.cover),random=list(site.id=~1),data=funct.scale)
AICc(mod.fc.arthro.gam,mod.fc.arthro.gam.s)

c19<-ggplot(eco.funct,aes(floral.cover,aug.worms.gone))+
  geom_point(size=.75)+
  stat_function(fun=function(x)(fixef(mod.fc.arthro)[1] + fixef(mod.fc.arthro)[2]*x),size=1,color="black")+
  annotate("text",label="r=-0.11",x=35,y=4,size=3)+
  labs(y = "Arthropod predation rate\n(worms/day)", x="Pollinator resource\navailability")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

d19<-ggplot(eco.funct,aes(floral.cover,aug.worms.gone))+
  geom_point(size=.75)+
  geom_smooth(method="gam",formula=y~s(x),se=F,color="black")+
  annotate("text",label="paste(R ^ 2, \" = 0.03*\")",parse=T,x=35,y=3.75,size=3)+
  labs(y = "Arthropod predation rate\n(worms/day)", x="Pollinator resource\navailability")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

summary(lmer(aug.worms.gone~floral.cover+(1|site.id),data=funct.scale))

mod.fr.arthro<-lmer(aug.worms.gone~floral.plotSR+(1|site.id),data=eco.funct)
mod.fr.arthro.gam<-gamm(aug.worms.gone~floral.plotSR,random=list(site.id=~1),data=funct.scale)
mod.fr.arthro.gam.s<-gamm(aug.worms.gone~s(floral.plotSR),random=list(site.id=~1),data=funct.scale)
AICc(mod.fr.arthro.gam,mod.fr.arthro.gam.s)

c20<-ggplot(eco.funct,aes(floral.plotSR,aug.worms.gone))+
  geom_point(size=.75)+
  stat_function(fun=function(x)(fixef(mod.fr.arthro)[1] + fixef(mod.fr.arthro)[2]*x),size=1,color="black")+
  annotate("text",label="r=0.01",x=0.5,y=3.5,size=3)+
  labs(y = "Arthropod predation rate\n(worms/day)", x="Floral richness")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

d20<-ggplot(eco.funct,aes(floral.plotSR,aug.worms.gone))+
  geom_point(size=.75)+
  geom_smooth(method="gam",formula=y~s(x),se=F,color="black")+
  annotate("text",label="paste(R ^ 2, \" = 0.002\")",parse=T,x=0.5,y=3.5,size=3)+
  labs(y = "Arthropod predation rate\n(worms/day)", x="Floral richness")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

summary(lmer(floral.plotSR~aug.worms.gone+(1|site.id),data=funct.scale))

mod.fc.fr<-lmer(floral.cover~floral.plotSR+(1|site.id),data=eco.funct)
mod.fc.fr.gam<-gamm(floral.cover~floral.plotSR,random=list(site.id=~1),data=funct.scale)
mod.fc.fr.gam.s<-gamm(floral.cover~s(floral.plotSR),random=list(site.id=~1),data=funct.scale)
AICc(mod.fc.fr.gam,mod.fc.fr.gam.s)

c21<-ggplot(eco.funct,aes(floral.plotSR,floral.cover))+
  geom_point(size=.75)+
  stat_function(fun=function(x)(fixef(mod.fc.fr)[1] + fixef(mod.fc.fr)[2]*x),size=1,color="black")+
  annotate("text",label="r=0.34*",x=0.5,y=50,size=3)+
  labs(y = "Pollinator resource\navailability", x="Floral richness")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

d21<-ggplot(eco.funct,aes(floral.plotSR,floral.cover))+
  geom_point(size=.75)+
  geom_smooth(method="gam",formula=y~s(x),se=F,color="black")+
  annotate("text",label="paste(R ^ 2, \" = 0.15*\")",parse=T,x=0.5,y=45,size=3)+
  labs(y = "Pollinator resource\navailability", x="Floral richness")+
  theme(text = element_text(size=10),axis.text=element_text(colour="black"),strip.background =element_rect(fill="white"),
        panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_line(size=1, color="black"),legend.position="none")

summary(lmer(floral.cover~floral.plotSR+(1|site.id),data=funct.scale))

#Create Fig. S2 "Linear tradeoffs" (uncomment to export PDF file)
#pdf("cor_fig.pdf",height=11,width = 8.5)
multiplot(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,cols=3)
dev.off()

#Create Fig. S3 "non-linear tradeoffs" (uncomment to export PDF file)
#pdf("cor_fig_gam.pdf",height=11,width = 8.5)
multiplot(d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12,d13,d14,d15,d16,d17,d18,d19,d20,d21,cols=3)
dev.off()

####Site level average models####
mod.mult.phylo.av<-lm(multifunct.a~ses.pd+Andropogon.gerardii+soilPC1+landscape+age.2013+fire.return,data=data.scale.av)
mod.mult.funct.av<-lm(multifunct.a~FRaoQ+Andropogon.gerardii+soilPC1+landscape+age.2013+fire.return,data=data.scale.av)
mod.mult.tax.av<-lm(multifunct.a~sim.div+Andropogon.gerardii+soilPC1+landscape+age.2013+fire.return,data=data.scale.av)
