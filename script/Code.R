setwd('~/Desktop/Project/6.OAW_interaction/2.Analysis_2021')
### packages-----
library(stringr)
library(dplyr)
library(tidyr)
library(metafor)
library(ggplot2)
library(devtools)
library(tidyverse) 
library(patchwork) 
library(R.rsp) 
library(tibble) 
library(cowplot) 
library(metaAidR) 
library(orchaRd) 
library(corpcor)
library(maps)
## note that the code for making orchard plots used the orchaRd 1.0., which may not rerun with the orchaRd 2.0 updated in 2023. Nakagawa 
dt=read.csv('OAW_interaction.csv')
new_orchard = function(a){ # for orchard plots
  b <- mod_results(a, mod = "Stressor") 
  orchard_plot(b, mod = "Stressor", xlab = "LnRR (main effect size)",
               angle = 0, cb = FALSE)+
    theme(panel.grid = element_blank())+coord_flip()+
    scale_x_continuous(limits=c(-1.6,1.6),breaks = c(-1.5,-1.0, -0.5,0,0.5,1,1.5))
}

##### 1.Generating variance-covariance matrices-------------------------------------------------
V_0.5 <- make_VCV_matrix(data = dt, V = "variance",cluster = "share_organism", 
                         obs = "Effect_size_ID", type = "vcv", rho = 0.5) 
V_0.9 <- make_VCV_matrix(data = dt, V = "variance",cluster = "share_organism", 
                         obs = "Effect_size_ID", type = "vcv", rho = 0.9)  # for sensitivity analysis
is.positive.definite(V_0.5) # TRUE
is.positive.definite(V_0.9) # TRUE
##### 2.Defining random effects structure----
# Effect_size_ID to estimate the residual heterogeneity,and StudyID to account for dependence between effect
# size from the same publication were always included.

## for whole dataset
RE1.dt <- rma.mv(yi = Main_effect, V = variance, 
                 random = list(~1 | StudyID, ~1 | Effect_size_ID,~1 | Species,~1 | Shared_first_author,
                               ~1 | Country), data = dt, method = "ML")
RE2.dt <- rma.mv(yi = Main_effect, V = variance, 
                 random = list(~1 | StudyID, ~1 | Effect_size_ID,~1 | Species,~1 | Shared_first_author), 
                 data = dt, method = "ML")
RE3.dt <- rma.mv(yi = Main_effect, V = variance, 
                 random = list(~1 | StudyID, ~1 | Effect_size_ID,~1 | Species,~1 | Country), 
                 data = dt, method = "ML")
RE4.dt <- rma.mv(yi = Main_effect, V = variance, 
                 random = list(~1 | StudyID, ~1 | Effect_size_ID,~1 | Shared_first_author,~1 | Country), 
                 data = dt, method = "ML")

RE5.dt <- rma.mv(yi = Main_effect, V = variance, 
                 random = list(~1 | StudyID, ~1 | Effect_size_ID,~1 | Shared_first_author), 
                 data = dt, method = "ML")
RE6.dt <- rma.mv(yi = Main_effect, V = variance, 
                 random = list(~1 | StudyID, ~1 | Effect_size_ID,~1 | Country), 
                 data = dt, method = "ML")
RE7.dt <- rma.mv(yi = Main_effect, V = variance, 
                 random = list(~1 | StudyID, ~1 | Effect_size_ID,~1 | Species), 
                 data = dt, method = "ML")
RE8.dt <- rma.mv(yi = Main_effect, V = variance, random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                 data = dt, method = "ML")
AIC(RE1.dt,RE2.dt,RE3.dt,RE4.dt,RE5.dt,RE6.dt,RE7.dt,RE8.dt)


##### FIG 1 creating the map-----
world_map = map_data('world')
mapplot=ggplot()+geom_polygon(data=world_map,aes(x=long,y=lat, group=group),fill='grey')+
  geom_point(data=dt, aes(x=Longitude, y=Latitude), color='blue')+theme_bw()+
  theme(panel.grid = element_blank())
##### FIG 2 overall OA,OW,and interaction effects on marine trophic levels----
dt= read.csv('OAW_interaction.csv')

# Predator unmerged using the moderator "TrophicLevel"
str_TL <- rma.mv(yi = Main_effect, V = V_0.5, mods = ~Stressor*TrophicLevel-1,random = list(~1 | StudyID, ~1 | Effect_size_ID), 
              data = dt, method = "REML")

res_RV <- mod_results(str_TL, mod = "Stressor*TrophicLevel",group = "Study_ID") 
box <- res_RV$mod_table
new_orchard(str_TL)+
  scale_x_discrete()+
  scale_color_manual(values=c(rep(c('#00BA38','#619CFF','#F8766D'),4)))+
  scale_fill_manual(values=c(rep(c('#00BA38','#619CFF','#F8766D'),4)))


unmerge_plot <- orchard_plot(res_RV, mod = "Stressor*TrophicLevel", xlab = "LnRR (main effect size)",
             angle = 0, cb = FALSE)+
  theme(panel.grid = element_blank())+coord_flip()+
  scale_x_continuous(limits=c(-1.6,1.6),breaks = c(-1.5,-1.0, -0.5,0,0.5,1,1.5))

# Predator merged using the moderator "Trophic_level"
str_tl=rma.mv(yi = Main_effect, V = V_0.5, mods = ~Stressor*Trophic_level-1,random = list(~1 | StudyID, ~1 | Effect_size_ID), 
              data = dt, method = "REML")
res_RV <- mod_results(str_tl, mod = "Stressor*Trophic_level") 
box <- res_RV$mod_table
new_orchard(str_TL)

merge_plot <- orchard_plot(res_RV, mod = "Stressor*Trophic_level", xlab = "LnRR (main effect size)",
             angle = 0, cb = FALSE)+
  theme(panel.grid = element_blank())+coord_flip()+
  scale_x_continuous(limits=c(-1.6,1.6),breaks = c(-1.5,-1.0, -0.5,0,0.5,1,1.5))

##### FIG 3A Effects on calcifiers at different trophic levels----
## calcifiers
# unmerged predators
calstr_TL=rma.mv(yi = Main_effect, V = V_0.5, mods = ~Stressor*TrophicLevel-1,
                 subset=(Calcifier=='y'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                 data = dt, method = "REML")

res_RV <- mod_results(calstr_TL, mod = "Stressor*TrophicLevel",group = "Study_ID") 
box <- res_RV$mod_table

orchard_plot(res_RV, mod = "Stressor*TrophicLevel", xlab = "LnRR (main effect size)",
             angle = 0, cb = FALSE)+
  theme(panel.grid = element_blank())+coord_flip()+
  scale_x_continuous(limits=c(-1.6,1.6),breaks = c(-1.5,-1.0, -0.5,0,0.5,1,1.5))

# Merged
calstr_tl=rma.mv(yi = Main_effect, V = V_0.5, mods = ~Stressor*Trophic_level-1,
                 subset=(Calcifier=='y'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                 data = dt, method = "REML")

res_RV <- mod_results(calstr_tl, mod = "Stressor*Trophic_level") 
box <- res_RV$mod_table

orchard_plot(res_RV, mod = "Stressor*Trophic_level", xlab = "LnRR (main effect size)",
             angle = 0, cb = FALSE)+
  theme(panel.grid = element_blank())+coord_flip()+
  scale_x_continuous(limits=c(-1.6,1.6),breaks = c(-1.5,-1.0, -0.5,0,0.5,1,1.5))

##### FIG 3B Effects on non-calcifiers at different trophic levels----
## non-calcifiers
# Unmerged
ncalstr_TL=rma.mv(yi = Main_effect, V = V_0.5, mods = ~Stressor*TrophicLevel-1,
                  subset=(Calcifier=='n'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                  data = dt, method = "REML")

res_RV <- mod_results(ncalstr_TL, mod = "Stressor*TrophicLevel") 
box <- res_RV$mod_table
orchard_plot(res_RV, mod = "Stressor*TrophicLevel", xlab = "LnRR (main effect size)",
             angle = 0, cb = FALSE)+
  theme(panel.grid = element_blank())+coord_flip()+
  scale_x_continuous(limits=c(-1.6,1.6),breaks = c(-1.5,-1.0, -0.5,0,0.5,1,1.5))

# Merged
ncalstr_tl=rma.mv(yi = Main_effect, V = V_0.5, mods = ~Stressor*Trophic_level-1,
                  subset=(Calcifier=='n'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                  data = dt, method = "REML")

res_RV <- mod_results(ncalstr_tl, mod = "Stressor*Trophic_level") 
box <- res_RV$mod_table
orchard_plot(res_RV, mod = "Stressor*Trophic_level", xlab = "LnRR (main effect size)",
             angle = 0, cb = FALSE)+
  theme(panel.grid = element_blank())+coord_flip()+
  scale_x_continuous(limits=c(-1.6,1.6),breaks = c(-1.5,-1.0, -0.5,0,0.5,1,1.5))

noncalstre <- rma.mv(yi = Main_effect, V = V_0.5, mods = ~Stressor-1,
                     subset=(Calcifier=='n'),
                     random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                     data = dt, method = "REML")

##### FIG 4 Relationship beteen latitudes and main effect sizes----
new_ggplot =function(a){
  ggplot(a, aes(x=Abs_lat, y=Main_effect))+
    geom_point(aes(size=1/sqrt(variance)))+
    geom_smooth(method = lm, alpha=0.3,color='Firebrick4',fill='Firebrick4')+
    scale_y_continuous(limits=c(-1.5,1),breaks = c(-1.5,-1.0, -0.5,0,0.5,1))+
    scale_x_continuous(limits=c(10,80),breaks = c(10,20,30,40,50,60,70,80))+
    theme_bw()+ylab("LnRR (main effect size)")+xlab('Absolute latitude')+
    annotate('text',label="Slope = 0.0011 p = 0.7846", x=60, y=-1.0)+
    theme(panel.grid = element_blank(),
          panel.border = element_rect(size=1.2),
          axis.text = element_text(size=14, color='black'),
          axis.title= element_text(size=16, color='black'),
          legend.direction = "horizontal",
          legend.position = c(0.5,0.9),
          legend.background = element_blank())+
    geom_hline(yintercept=0, linetype='dashed', alpha=0.5)
} # This function is used for plots in the main text with purple color
LatPlot = function(a,b){ # a is the dataset, b is the model
  changeLine <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
  colnames(changeLine) <- c('X','Y','Upper','Lower')
  tempX <- seq(min(a$Abs_lat),max(a$Abs_lat),length.out = 100)
  for(i in seq(1,length(tempX))){
    changeLine[i,1] <- tempX[i]
    tempResult <- predict.rma(b,newmods = tempX[i])
    changeLine[i,2] <- as.numeric(tempResult$pred)
    changeLine[i,3] <- as.numeric(tempResult$ci.ub)
    changeLine[i,4] <- as.numeric(tempResult$ci.lb)
  }
  # CI polygon
  changeCIPolygon <- as.data.frame(matrix(data = 0,nrow = (2*nrow(changeLine)+1),ncol = 2))
  colnames(changeCIPolygon) <- c('X','Y')
  count <- 1
  for(i in seq(1,nrow(changeLine))){
    changeCIPolygon$X[count] <- changeLine$X[i]
    changeCIPolygon$Y[count] <- changeLine$Upper[i]
    count <- count + 1
  }
  for(i in seq(nrow(changeLine),1,-1)){
    changeCIPolygon$X[count] <- changeLine$X[i]
    changeCIPolygon$Y[count] <- changeLine$Lower[i]
    count <- count + 1
  }
  changeCIPolygon$X[count] <- changeLine$X[1]
  changeCIPolygon$Y[count] <- changeLine$Upper[1]
  
  a %>%
    ggplot(aes(x = Abs_lat, y = Main_effect, size = sqrt(1/vi))) + 
    geom_point(shape = 21, fill="#00B9E3", color='#00B9E3',alpha=0.75) +
    geom_line(data = changeLine,mapping = aes(x = X,y = Y),size = 1,color = 'black')+
    geom_polygon(data = changeCIPolygon,mapping = aes(x = X,y = Y),size = 0.5,fill = "#00B9E3",alpha = 0.3)+ 
    labs(x = "Latitude", y = "lnRR (effect size)", size = "Precision (1/SE)")+
    geom_hline(yintercept=0, linetype='dashed', alpha=0.5)+ 
    guides(fill = "none", colour = "none")+theme_bw()+
    theme(panel.grid = element_blank(),
          panel.border = element_rect(size=1.2),
          axis.text = element_text(size=14, color='black'),
          axis.title= element_text(size=16, color='black'),
          legend.direction = "horizontal",
          legend.position = c(0.5,0.9),
          legend.background = element_blank())+
    scale_y_continuous(limits=c(-1.5,1.5),breaks = c(-1.5,-1.0, -0.5,0,0.5,1.0,1.5))+
    scale_x_continuous(limits=c(10,80),breaks = c(10,20,30,40,50,60,70,80))
} # for new plot with two arguments, one is the dataset, the other is the regression model


# OA effects on the primary producer
ppOA = filter(dt, TrophicLevel=='primary producer'& Stressor=='OA')

# meta-regression model using the absolute latitude as the moderator
pp_oa = rma.mv(yi = Main_effect, V = V_0.5, mods = ~Abs_lat,
             subset=(TrophicLevel=='primary producer' & Stressor=='OA'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
             data = dt, method = "REML")

LatPlot(ppOA, pp_oa)

# OA effects on the herbivore
herOA = filter(dt, TrophicLevel=='herbivore'& Stressor=='OA')
her_oa = rma.mv(yi = Main_effect, V = V_0.5, mods = ~Abs_lat,
                subset=(TrophicLevel=='herbivore' & Stressor=='OA'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                data = dt, method = "REML")
LatPlot(herOA, her_oa)

# OA effects on the Meso-predator
mOA = filter(dt, TrophicLevel=='meso-predator'& Stressor=='OA')

m_oa=rma.mv(yi = Main_effect, V = V_0.5, mods = ~Abs_lat,
            subset=(TrophicLevel=='meso-predator' & Stressor=='OA'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
            data = dt, method = "REML")
LatPlot(mOA, m_oa)

# OA effects on the top-predator
topOA <- filter(dt, TrophicLevel=='top-predator'& Stressor=='OA')
top_oa <- rma.mv(yi = Main_effect, V = V_0.5, mods = ~Abs_lat,
       subset=(TrophicLevel=='top-predator' & Stressor=='OA'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
       data = dt, method = "REML")
LatPlot(topOA, top_oa)


# OW effects on the primary producer
ppOW <- filter(dt, Stressor=='OW'& TrophicLevel=='primary producer')
herOW <- filter(dt, Stressor=='OW'& TrophicLevel=='herbivore')
mesoOW <- filter(dt, Stressor=='OW'& TrophicLevel=='meso-predator')
topOW <- filter(dt, Stressor=='OW'& TrophicLevel=='top-predator')

pp_OW <- rma.mv(yi = Main_effect, V = V_0.5, mods = ~Abs_lat,
       subset=(TrophicLevel=='primary producer' & Stressor=='OW'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
       data = dt, method = "REML")
LatPlot(ppOW, pp_OW)

# OW effects on the herbivore
her_OW <- rma.mv(yi = Main_effect, V = V_0.5, mods = ~Abs_lat,
       subset=(TrophicLevel=='herbivore' & Stressor=='OW'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
       data = dt, method = "REML")
LatPlot(herOW, her_OW)

# OW effects on the meso-predator
meso_OW <- rma.mv(yi = Main_effect, V = V_0.5, mods = ~Abs_lat,
       subset=(TrophicLevel=='meso-predator' & Stressor=='OW'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
       data = dt, method = "REML")
LatPlot(mesoOW, meso_OW)

# OW effects on the top-predator
top_OW <- rma.mv(yi = Main_effect, V = V_0.5, mods = ~Abs_lat,
       subset=(TrophicLevel=='top-predator' & Stressor=='OW'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
       data = dt, method = "REML")
LatPlot(topOW, top_OW)


# Interaction effects on the primary producer
int_pp <- filter(dt, Stressor=='interaction'& TrophicLevel=='primary producer')
int_her <- filter(dt, Stressor=='interaction'& TrophicLevel=='herbivore')
int_meso <- filter(dt, Stressor=='interaction'& TrophicLevel=='meso-predator')
int_top <- filter(dt, Stressor=='interaction'& TrophicLevel=='top-predator')
ppINT <- rma.mv(yi = Main_effect, V = V_0.5, mods = ~Abs_lat,
       subset=(TrophicLevel=='primary producer' & Stressor=='interaction'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
       data = dt, method = "REML")
LatPlot(int_pp, ppINT)

# Interaction effects on the herbivore
herINT <- rma.mv(yi = Main_effect, V = V_0.5, mods = ~Abs_lat,
       subset=(TrophicLevel=='herbivore' & Stressor=='interaction'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
       data = dt, method = "REML")
summary(herINT)

LatPlot(int_her, herINT)

# Interaction effects on the meso-predator
mesoINT <- rma.mv(yi = Main_effect, V = V_0.5, mods = ~Abs_lat,
       subset=(TrophicLevel=='meso-predator' & Stressor=='interaction'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
       data = dt, method = "REML")
LatPlot(int_meso, mesoINT)

# Interaction effects on the top-predator
topINT <- rma.mv(yi = Main_effect, V = V_0.5, mods = ~Abs_lat,
       subset=(TrophicLevel=='top-predator' & Stressor=='interaction'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
       data = dt, method = "REML")
LatPlot(int_top, topINT)


# OA OW and their interaction effects on the merged Predators
oa_pre <- filter(dt, Stressor=='OA'& Trophic_level=='predator')
ow_pre <- filter(dt, Stressor=='OW'& Trophic_level=='predator')
int_pre <- filter(dt, Stressor=='interaction'& Trophic_level=='predator')

preOA <- rma.mv(yi = Main_effect, V = V_0.5, mods = ~Abs_lat,
       subset=(Trophic_level=='predator' & Stressor=='OA'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
       data = dt, method = "REML")
preOW <- rma.mv(yi = Main_effect, V = V_0.5, mods = ~Abs_lat,
       subset=(Trophic_level=='predator' & Stressor=='OW'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
       data = dt, method = "REML")

preINT <- rma.mv(yi = Main_effect, V = V_0.5, mods = ~Abs_lat,
       subset=(Trophic_level=='predator' & Stressor=='interaction'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
       data = dt, method = "REML")
##plotting
LatPlot(oa_pre, preOA)
LatPlot(ow_pre, preOW)
LatPlot(int_pre, preINT)

##### FIG 5 Effects on climate regions----
tropic <- filter(dt, Region2=='Tropical')
V_0.5_tropic <- make_VCV_matrix(data = tropic, V = "variance",cluster = "share_organism",
                                obs = "Effect_size_ID", type = "vcv", rho = 0.5) 
tropic_model <- rma.mv(yi = Main_effect, V = V_0.5_tropic, mods = ~Stressor-1,random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                       data = tropic, method = "REML")

orchard_plot(tropic_model, mod = "Stressor", group = "Study_ID", xlab = "LnRR (main effect size)",transfm = "none", flip = FALSE)+
  theme(panel.grid = element_blank())+ scale_x_discrete(limits=c('OA','OW', 'Interaction'))+
  scale_y_continuous(limits=c(-2,2),breaks = c(-2.0, -1.0, 0, 1.0,2.0))


sub_tropic <- filter(dt, Region2=='Sub-tropical')
V_0.5_sub_tropical <- make_VCV_matrix(data = sub_tropic, V = "variance",cluster = "share_organism", 
                                      obs = "Effect_size_ID", type = "vcv", rho = 0.5) 
sub_t_model <- rma.mv(yi = Main_effect, V = V_0.5_sub_tropical, mods = ~Stressor-1,random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                      data = sub_tropic, method = "REML")

orchard_plot(sub_t_model, mod = "Stressor", group = "Study_ID", xlab = "LnRR (main effect size)",transfm = "none", flip = FALSE)+
  theme(panel.grid = element_blank())+ scale_x_discrete(limits=c('OA','OW', 'Interaction'))+
  scale_y_continuous(limits=c(-1,1))


temperate <- filter(dt, Region2=='Temperate')

V_0.5_temperate <- make_VCV_matrix(data = temperate, V = "variance",cluster = "share_organism", 
                                   obs = "Effect_size_ID", type = "vcv", rho = 0.5) 
temperate_model <- rma.mv(yi = Main_effect, V = V_0.5_temperate, mods = ~Stressor-1,random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                          data = temperate, method = "REML")

orchard_plot(temperate_model, mod = "Stressor", group = "Study_ID", xlab = "LnRR (main effect size)",transfm = "none", flip = FALSE)+
  theme(panel.grid = element_blank())+ scale_x_discrete(limits=c('OA','OW', 'Interaction'))+
  scale_y_continuous(limits=c(-1,1))


##### Publication bias----
#funnel plot
# selecting the best random effects combination
PB1 = rma.mv(yi = Main_effect, V = V_0.5,mods = ~Stressor-1,random = list(~1 | StudyID, ~1 | Effect_size_ID), 
             data = dt, method = "ML")
PB2 = rma.mv(yi = Main_effect, V = V_0.5,mods = ~TrophicLevel-1,random = list(~1 | StudyID, ~1 | Effect_size_ID), 
             data = dt, method = "ML")
PB3 = rma.mv(yi = Main_effect, V = V_0.5,mods = ~Calcifier-1,random = list(~1 | StudyID, ~1 | Effect_size_ID), 
             data = dt, method = "ML")
PB4 = rma.mv(yi = Main_effect, V = V_0.5,mods = ~Calcifier+TrophicLevel-1,random = list(~1 | StudyID, ~1 | Effect_size_ID), 
             data = dt, method = "ML")
PB5 = rma.mv(yi = Main_effect, V = V_0.5,mods = ~Calcifier+Stressor-1,random = list(~1 | StudyID, ~1 | Effect_size_ID), 
             data = dt, method = "ML")
PB6 = rma.mv(yi = Main_effect, V = V_0.5,mods = ~TrophicLevel+Stressor-1,random = list(~1 | StudyID, ~1 | Effect_size_ID), 
             data = dt, method = "ML")
PB7 = rma.mv(yi = Main_effect, V = V_0.5,mods = ~TrophicLevel+Stressor+Calcifier-1,random = list(~1 | StudyID, ~1 | Effect_size_ID), 
             data = dt, method = "ML")
PB8 = rma.mv(yi = Main_effect, V = V_0.5,mods = ~TrophicLevel*Stressor-1,random = list(~1 | StudyID, ~1 | Effect_size_ID), 
             data = dt, method = "ML")
PB9 = rma.mv(yi = Main_effect, V = V_0.5,mods = ~Calcifier*Stressor-1,random = list(~1 | StudyID, ~1 | Effect_size_ID), 
             data = dt, method = "ML")
AIC(PB1,PB2,PB3,PB4,PB5,PB6,PB7,PB8,PB9)


funnel <- rma.mv(yi = Main_effect, V = V_0.5,mods = ~Calcifier-1,random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                 data = dt, method = "REML")
funnel(funnel, yaxis = "seinv", level = c(90, 95, 99),
       ylim = c(0.2, 2.4), shade = c("white", "gray55", "gray75"),
       refline = 0, legend = TRUE)
egger <- rma.mv(yi = Main_effect, V = V_0.5, mods = ~Calcifier+sqrt(variance),
                random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                data = dt, method = "REML")
summary(egger) # p=0.4569 no significant funnel asymmetry

which(dt2==max(dt2$variance),arr.ind=T)
nrow(dt)
nrow(dt2)
dt2=dt2[-182,]
ggplot(dt2, aes(x=sqrt(variance), y=Main_effect))+
  geom_point(aes(size=1/sqrt(variance)))+
  geom_smooth(method = lm, alpha=0.3,color='Firebrick4',fill='Firebrick4')+
  theme_bw()+ylab("LnRR (main effect size)")+xlab('sqrt(variance)')+ylim(-5,5)+
  annotate('text',label="Slope = -0.137 p = 0.457", x=0.75, y=-5.0)+
  theme(panel.grid = element_blank(),
        panel.border = element_rect(size=1.2),
        axis.text = element_text(size=14, color='black'),
        axis.title= element_text(size=16, color='black'),
        legend.direction = "horizontal",
        legend.position = c(0.5,0.9),
        legend.background = element_blank())+
  geom_hline(yintercept=0, linetype='dashed', alpha=0.5)



##### Sensitivity analysis with setting the correlation coefficient to 0.9----
## 1.Generating variance-covariance matrices
V_0.9 <- make_VCV_matrix(data = dt, V = "variance",cluster = "share_organism", 
                         obs = "Effect_size_ID", type = "vcv", rho = 0.9) 

is.positive.definite(V_0.9) # TRUE

##### Sensitivity analysis to FIG 2
dt= read.csv('OAW_interaction.csv')

# Predator unmerged
str_TL <- rma.mv(yi = Main_effect, V = V_0.9, mods = ~Stressor*TrophicLevel-1,random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                 data = dt, method = "REML")

res_RV <- mod_results(str_TL, mod = "Stressor*TrophicLevel",group = "Study_ID") 
box <- res_RV$mod_table
new_orchard(str_TL)+
  scale_x_discrete()+
  scale_color_manual(values=c(rep(c('#00BA38','#619CFF','#F8766D'),4)))+
  scale_fill_manual(values=c(rep(c('#00BA38','#619CFF','#F8766D'),4)))


unmerge_plot <- orchard_plot(res_RV, mod = "Stressor*TrophicLevel", xlab = "LnRR (main effect size)",
                             angle = 0, cb = FALSE)+
  theme(panel.grid = element_blank())+coord_flip()+
  scale_x_continuous(limits=c(-1.6,1.6),breaks = c(-1.5,-1.0, -0.5,0,0.5,1,1.5))

# Predator merged
str_tl=rma.mv(yi = Main_effect, V = V_0.9, mods = ~Stressor*Trophic_level-1,random = list(~1 | StudyID, ~1 | Effect_size_ID), 
              data = dt, method = "REML")
res_RV <- mod_results(str_tl, mod = "Stressor*Trophic_level") 
box <- res_RV$mod_table
new_orchard(str_TL)

merge_plot <- orchard_plot(res_RV, mod = "Stressor*Trophic_level", xlab = "LnRR (main effect size)",
                           angle = 0, cb = FALSE)+
  theme(panel.grid = element_blank())+coord_flip()+
  scale_x_continuous(limits=c(-1.6,1.6),breaks = c(-1.5,-1.0, -0.5,0,0.5,1,1.5))


##### Sensitivity analysis to FIG 3A
#calcifiers
# unmerged
calstr_TL=rma.mv(yi = Main_effect, V = V_0.9, mods = ~Stressor*TrophicLevel-1,
                 subset=(Calcifier=='y'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                 data = dt, method = "REML")

res_RV <- mod_results(calstr_TL, mod = "Stressor*TrophicLevel",group = "Study_ID") 
box <- res_RV$mod_table
orchard_plot(res_RV, mod = "Stressor*TrophicLevel", xlab = "LnRR (main effect size)",
             angle = 0, cb = FALSE)+
  theme(panel.grid = element_blank())+coord_flip()+
  scale_x_continuous(limits=c(-1.6,1.6),breaks = c(-1.5,-1.0, -0.5,0,0.5,1,1.5))

# Merged
calstr_tl=rma.mv(yi = Main_effect, V = V_0.9, mods = ~Stressor*Trophic_level-1,
                 subset=(Calcifier=='y'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                 data = dt, method = "REML")

res_RV <- mod_results(calstr_tl, mod = "Stressor*Trophic_level") 
box <- res_RV$mod_table
orchard_plot(res_RV, mod = "Stressor*Trophic_level", xlab = "LnRR (main effect size)",
             angle = 0, cb = FALSE)+
  theme(panel.grid = element_blank())+coord_flip()+
  scale_x_continuous(limits=c(-1.6,1.6),breaks = c(-1.5,-1.0, -0.5,0,0.5,1,1.5))



##### Sensitivity analysis to FIG 3B
# non-calcifiers
ncalstr_TL=rma.mv(yi = Main_effect, V = V_0.9, mods = ~Stressor*TrophicLevel-1,
                  subset=(Calcifier=='n'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                  data = dt, method = "REML")

res_RV <- mod_results(ncalstr_TL, mod = "Stressor*TrophicLevel") 
box <- res_RV$mod_table
orchard_plot(res_RV, mod = "Stressor*TrophicLevel", xlab = "LnRR (main effect size)",
             angle = 0, cb = FALSE)+
  theme(panel.grid = element_blank())+coord_flip()+
  scale_x_continuous(limits=c(-1.6,1.6),breaks = c(-1.5,-1.0, -0.5,0,0.5,1,1.5))

#Merged
ncalstr_tl=rma.mv(yi = Main_effect, V = V_0.9, mods = ~Stressor*Trophic_level-1,
                  subset=(Calcifier=='n'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                  data = dt, method = "REML")

res_RV <- mod_results(ncalstr_tl, mod = "Stressor*Trophic_level") 
box <- res_RV$mod_table
orchard_plot(res_RV, mod = "Stressor*Trophic_level", xlab = "LnRR (main effect size)",
             angle = 0, cb = FALSE)+
  theme(panel.grid = element_blank())+coord_flip()+
  scale_x_continuous(limits=c(-1.6,1.6),breaks = c(-1.5,-1.0, -0.5,0,0.5,1,1.5))

noncalstre <- rma.mv(yi = Main_effect, V = V_0.9, mods = ~Stressor-1,
                     subset=(Calcifier=='n'),
                     random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                     data = dt, method = "REML")

##### Sensitivity analysis to FIG 4
new_ggplot =function(a){
  ggplot(a, aes(x=Abs_lat, y=Main_effect))+
    geom_point(aes(size=1/sqrt(variance)))+
    geom_smooth(method = lm, alpha=0.3,color='Firebrick4',fill='Firebrick4')+
    scale_y_continuous(limits=c(-1.5,1),breaks = c(-1.5,-1.0, -0.5,0,0.5,1))+
    scale_x_continuous(limits=c(10,80),breaks = c(10,20,30,40,50,60,70,80))+
    theme_bw()+ylab("LnRR (main effect size)")+xlab('Absolute latitude')+
    annotate('text',label="Slope = 0.0011 p = 0.7846", x=60, y=-1.0)+
    theme(panel.grid = element_blank(),
          panel.border = element_rect(size=1.2),
          axis.text = element_text(size=14, color='black'),
          axis.title= element_text(size=16, color='black'),
          legend.direction = "horizontal",
          legend.position = c(0.5,0.9),
          legend.background = element_blank())+
    geom_hline(yintercept=0, linetype='dashed', alpha=0.5)
} # This function is used for plots in the main text with purple color
LatPlot = function(a,b){ # a is the dataset, b is the model
  changeLine <- as.data.frame(matrix(data = 0,nrow = 100,ncol = 4))
  colnames(changeLine) <- c('X','Y','Upper','Lower')
  tempX <- seq(min(a$Abs_lat),max(a$Abs_lat),length.out = 100)
  for(i in seq(1,length(tempX))){
    changeLine[i,1] <- tempX[i]
    tempResult <- predict.rma(b,newmods = tempX[i])
    changeLine[i,2] <- as.numeric(tempResult$pred)
    changeLine[i,3] <- as.numeric(tempResult$ci.ub)
    changeLine[i,4] <- as.numeric(tempResult$ci.lb)
  }
  # CI polygon
  changeCIPolygon <- as.data.frame(matrix(data = 0,nrow = (2*nrow(changeLine)+1),ncol = 2))
  colnames(changeCIPolygon) <- c('X','Y')
  count <- 1
  for(i in seq(1,nrow(changeLine))){
    changeCIPolygon$X[count] <- changeLine$X[i]
    changeCIPolygon$Y[count] <- changeLine$Upper[i]
    count <- count + 1
  }
  for(i in seq(nrow(changeLine),1,-1)){
    changeCIPolygon$X[count] <- changeLine$X[i]
    changeCIPolygon$Y[count] <- changeLine$Lower[i]
    count <- count + 1
  }
  changeCIPolygon$X[count] <- changeLine$X[1]
  changeCIPolygon$Y[count] <- changeLine$Upper[1]
  
  a %>%
    ggplot(aes(x = Abs_lat, y = Main_effect, size = sqrt(1/vi))) + 
    geom_point(shape = 21, fill="#00B9E3", color='#00B9E3',alpha=0.75) +
    geom_line(data = changeLine,mapping = aes(x = X,y = Y),size = 1,color = 'black')+
    geom_polygon(data = changeCIPolygon,mapping = aes(x = X,y = Y),size = 0.5,fill = "#00B9E3",alpha = 0.3)+ 
    labs(x = "Latitude", y = "lnRR (effect size)", size = "Precision (1/SE)")+
    geom_hline(yintercept=0, linetype='dashed', alpha=0.5)+ 
    guides(fill = "none", colour = "none")+theme_bw()+
    theme(panel.grid = element_blank(),
          panel.border = element_rect(size=1.2),
          axis.text = element_text(size=14, color='black'),
          axis.title= element_text(size=16, color='black'),
          legend.direction = "horizontal",
          legend.position = c(0.5,0.9),
          legend.background = element_blank())+
    scale_y_continuous(limits=c(-1.5,1.5),breaks = c(-1.5,-1.0, -0.5,0,0.5,1.0,1.5))+
    scale_x_continuous(limits=c(10,80),breaks = c(10,20,30,40,50,60,70,80))
} # for new plot with two arguments, one is the dataset, the other is the regression model


# OA effects on the primary producer
ppOA = filter(dt, TrophicLevel=='primary producer'& Stressor=='OA')

# meta-regression model using the absolute latitude as the moderator
pp_oa = rma.mv(yi = Main_effect, V = V_0.9, mods = ~Abs_lat,
               subset=(TrophicLevel=='primary producer' & Stressor=='OA'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
               data = dt, method = "REML")

LatPlot(ppOA, pp_oa)

# OA effects on the herbivore
herOA = filter(dt, TrophicLevel=='herbivore'& Stressor=='OA')
her_oa = rma.mv(yi = Main_effect, V = V_0.9, mods = ~Abs_lat,
                subset=(TrophicLevel=='herbivore' & Stressor=='OA'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                data = dt, method = "REML")
LatPlot(herOA, her_oa)

# OA effects on the Meso-predator
mOA = filter(dt, TrophicLevel=='meso-predator'& Stressor=='OA')

m_oa=rma.mv(yi = Main_effect, V = V_0.9, mods = ~Abs_lat,
            subset=(TrophicLevel=='meso-predator' & Stressor=='OA'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
            data = dt, method = "REML")
LatPlot(mOA, m_oa)

# OA effects on the top-predator
topOA <- filter(dt, TrophicLevel=='top-predator'& Stressor=='OA')
top_oa <- rma.mv(yi = Main_effect, V = V_0.9, mods = ~Abs_lat,
                 subset=(TrophicLevel=='top-predator' & Stressor=='OA'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                 data = dt, method = "REML")
LatPlot(topOA, top_oa)


# OW effects on the primary producer
ppOW <- filter(dt, Stressor=='OW'& TrophicLevel=='primary producer')
herOW <- filter(dt, Stressor=='OW'& TrophicLevel=='herbivore')
mesoOW <- filter(dt, Stressor=='OW'& TrophicLevel=='meso-predator')
topOW <- filter(dt, Stressor=='OW'& TrophicLevel=='top-predator')

pp_OW <- rma.mv(yi = Main_effect, V = V_0.9, mods = ~Abs_lat,
                subset=(TrophicLevel=='primary producer' & Stressor=='OW'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                data = dt, method = "REML")
LatPlot(ppOW, pp_OW)

# OW effects on the herbivore
her_OW <- rma.mv(yi = Main_effect, V = V_0.9, mods = ~Abs_lat,
                 subset=(TrophicLevel=='herbivore' & Stressor=='OW'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                 data = dt, method = "REML")
LatPlot(herOW, her_OW)

# OW effects on the meso-predator
meso_OW <- rma.mv(yi = Main_effect, V = V_0.9, mods = ~Abs_lat,
                  subset=(TrophicLevel=='meso-predator' & Stressor=='OW'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                  data = dt, method = "REML")
LatPlot(mesoOW, meso_OW)

# OW effects on the top-predator
top_OW <- rma.mv(yi = Main_effect, V = V_0.9, mods = ~Abs_lat,
                 subset=(TrophicLevel=='top-predator' & Stressor=='OW'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                 data = dt, method = "REML")
LatPlot(topOW, top_OW)


# Interaction effects on the primary producer
int_pp <- filter(dt, Stressor=='interaction'& TrophicLevel=='primary producer')
int_her <- filter(dt, Stressor=='interaction'& TrophicLevel=='herbivore')
int_meso <- filter(dt, Stressor=='interaction'& TrophicLevel=='meso-predator')
int_top <- filter(dt, Stressor=='interaction'& TrophicLevel=='top-predator')
ppINT <- rma.mv(yi = Main_effect, V = V_0.9, mods = ~Abs_lat,
                subset=(TrophicLevel=='primary producer' & Stressor=='interaction'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                data = dt, method = "REML")
LatPlot(int_pp, ppINT)

# Interaction effects on the herbivore
herINT <- rma.mv(yi = Main_effect, V = V_0.9, mods = ~Abs_lat,
                 subset=(TrophicLevel=='herbivore' & Stressor=='interaction'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                 data = dt, method = "REML")
summary(herINT)

LatPlot(int_her, herINT)

# Interaction effects on the meso-predator
mesoINT <- rma.mv(yi = Main_effect, V = V_0.9, mods = ~Abs_lat,
                  subset=(TrophicLevel=='meso-predator' & Stressor=='interaction'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                  data = dt, method = "REML")
LatPlot(int_meso, mesoINT)

# Interaction effects on the top-predator
topINT <- rma.mv(yi = Main_effect, V = V_0.9, mods = ~Abs_lat,
                 subset=(TrophicLevel=='top-predator' & Stressor=='interaction'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                 data = dt, method = "REML")
LatPlot(int_top, topINT)


# OA OW and their interaction effects on the merged Predators
oa_pre <- filter(dt, Stressor=='OA'& Trophic_level=='predator')
ow_pre <- filter(dt, Stressor=='OW'& Trophic_level=='predator')
int_pre <- filter(dt, Stressor=='interaction'& Trophic_level=='predator')

preOA <- rma.mv(yi = Main_effect, V = V_0.9, mods = ~Abs_lat,
                subset=(Trophic_level=='predator' & Stressor=='OA'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                data = dt, method = "REML")
preOW <- rma.mv(yi = Main_effect, V = V_0.9, mods = ~Abs_lat,
                subset=(Trophic_level=='predator' & Stressor=='OW'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                data = dt, method = "REML")

preINT <- rma.mv(yi = Main_effect, V = V_0.9, mods = ~Abs_lat,
                 subset=(Trophic_level=='predator' & Stressor=='interaction'),random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                 data = dt, method = "REML")
##plotting
LatPlot(oa_pre, preOA)
LatPlot(ow_pre, preOW)
LatPlot(int_pre, preINT)

##### Sensitivity analysis to FIG 5
tropic <- filter(dt, Region2=='Tropical')
V_0.9_tropic <- make_VCV_matrix(data = tropic, V = "variance",cluster = "share_organism",
                                obs = "Effect_size_ID", type = "vcv", rho = 0.9) 
tropic_model <- rma.mv(yi = Main_effect, V = V_0.9_tropic, mods = ~Stressor-1,random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                       data = tropic, method = "REML")

orchard_plot(tropic_model, mod = "Stressor", group = "Study_ID", xlab = "LnRR (main effect size)",transfm = "none", flip = FALSE)+
  theme(panel.grid = element_blank())+ scale_x_discrete(limits=c('OA','OW', 'Interaction'))+
  scale_y_continuous(limits=c(-2,2),breaks = c(-2.0, -1.0, 0, 1.0,2.0))


sub_tropic <- filter(dt, Region2=='Sub-tropical')
V_0.9_sub_tropical <- make_VCV_matrix(data = sub_tropic, V = "variance",cluster = "share_organism", 
                                      obs = "Effect_size_ID", type = "vcv", rho = 0.9) 
sub_t_model <- rma.mv(yi = Main_effect, V = V_0.9_sub_tropical, mods = ~Stressor-1,random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                      data = sub_tropic, method = "REML")

orchard_plot(sub_t_model, mod = "Stressor", group = "Study_ID", xlab = "LnRR (main effect size)",transfm = "none", flip = FALSE)+
  theme(panel.grid = element_blank())+ scale_x_discrete(limits=c('OA','OW', 'Interaction'))+
  scale_y_continuous(limits=c(-1,1))


temperate <- filter(dt, Region2=='Temperate')

V_0.9_temperate <- make_VCV_matrix(data = temperate, V = "variance",cluster = "share_organism", 
                                   obs = "Effect_size_ID", type = "vcv", rho = 0.9) 
temperate_model <- rma.mv(yi = Main_effect, V = V_0.9_temperate, mods = ~Stressor-1,random = list(~1 | StudyID, ~1 | Effect_size_ID), 
                          data = temperate, method = "REML")

orchard_plot(temperate_model, mod = "Stressor", group = "Study_ID", xlab = "LnRR (main effect size)",transfm = "none", flip = FALSE)+
  theme(panel.grid = element_blank())+ scale_x_discrete(limits=c('OA','OW', 'Interaction'))+
  scale_y_continuous(limits=c(-1,1))












