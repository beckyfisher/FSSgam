# A simple function for full subsets multiple regression in ecology with R
# 
# R. Fisher
# S.K. Wilson
# S.M. Sin
# A.C. Lee
# T.J. Langlois

# Reproducible example for:
# Case Study 2: The role of large reef-associated predators in structuring adjacent soft-sediment communities

# A re-analysis of data presented in:
#   Langlois, T. J., M. J. Anderson, and R. C. Babcock. 2005. Reef-associated predators influence adjacent soft-sediment communities. Ecology 86: 1508–1519.

# Script information----
# This script is designed to:
# 1. Import and format importance scores 
# 2. Make heat map of the importance scores - which indicates those variables in the most parsimonious model


# librarys----
detach("package:plyr", unload=TRUE)#will error - no worries
library(tidyr)
library(dplyr)
library(ggplot2)


rm(list=ls())
study<-"Clams"

# Where the data sits----
# Set your specific work.dir here-
work.dir=("~/Dropbox/FSS_paper/case_study_TL clams") # for Tim's mac

model.out=paste(work.dir,"ModelOut",sep="/")
plots="~/Google Drive/Manuscripts/Manuscript_Fisher_FSS/Tim-Figures"


# Import importance scores ----
setwd(model.out)
dir()
dat.taxa<-read.csv("clams_all.var.imp.csv")%>%
  rename(resp.var=X)%>%
  gather(key=predictor,value=importance,2:ncol(.))
head(dat.taxa,5)


# Plotting defaults----

# Theme-
Theme1 <-
  theme( # use theme_get() to see available options
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    legend.background = element_rect(fill="white"),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    legend.text = element_text(size=8),
    legend.title = element_text(size=8, face="bold"),
    legend.position = "top",
    legend.direction="horizontal",
    text=element_text(size=10),
    strip.text.y = element_text(size = 10,angle = 0),
    axis.title.x=element_text(vjust=0.3, size=10),
    axis.title.y=element_text(vjust=0.6, angle=90, size=10),
    axis.text.x=element_text(size=10,angle = 90, hjust=1,vjust=0.5),
    axis.text.y=element_text(size=10,face="italic"),
    axis.line.x=element_line(colour="black", size=0.5,linetype='solid'),
    axis.line.y=element_line(colour="black", size=0.5,linetype='solid'),
    strip.background = element_blank())

# colour ramps-
re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)

# Labels-
legend_title<-"Importance"

# Annotations-
dat.taxa.label<-dat.taxa%>%
  mutate(label=NA)%>%
  mutate(label=ifelse(predictor=="Distance"&resp.var=="BDS","X",ifelse(predictor=="Status"&resp.var=="BDS","X",ifelse(predictor=="sqrt.X500um"&resp.var=="BDS","X",label))))%>%
  mutate(label=ifelse(predictor=="lobster"&resp.var=="BMS","X",label))%>%
  mutate(label=ifelse(predictor=="sqrt.X4mm"&resp.var=="CPN","X",ifelse(predictor=="lobster"&resp.var=="CPN","X",label)))
head(dat.taxa.label,2)


# Plot gg.importance.scores ----
setwd(plots)
dev.off()
gg.importance.scores <- ggplot(dat.taxa.label, aes(x=predictor,y=resp.var,fill=importance))+ 
  geom_tile(show.legend=T) + 
  scale_fill_gradientn(legend_title,colours=c("white", re), na.value = "grey98",
                       limits = c(0, max(dat.taxa.label$importance)))+
  scale_x_discrete(limits=c("Distance",
                            "Status",
                            "lobster",
                            "snapper",
                            "fetch",
                            "org",
                            "sqrt.X4mm",
                            "sqrt.X2mm",
                            "sqrt.X1mm",
                            "sqrt.X500um"),
               labels=c(
                 "Distance",
                 "Status",
                 "Lobster",
                 "Snapper",
                 "Fetch (km)",
                 "Organic content",
                 "Grain size: 4mm",
                 "            2mm",
                 "            1mm",
                 "            500um"
               ))+
scale_y_discrete(limits = c("CPN",
                            "BMS",
                            "BDS"),
                 labels=c("P. novizelandiae",
                          "M. striata",
                          "D. subrosea"))+
  xlab(NULL)+
  ylab(NULL)+
  theme_classic()+
  Theme1+
  geom_text(aes(label=label))
gg.importance.scores

ggsave("Langlois.importance.scores.png",width = 15, height = 7,units = "cm")
