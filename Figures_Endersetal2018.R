# This code will recreate figures from Enders et al 2018
# Necessary files are listed on lines 34-36
# Working directory is set on line 31. Files read and written here.
# Run to line 73 first. Then, each following sections uses the "long_data" variable to create a figure, up to line 545.
# Rest of script starting at line 546 makes Figure 3 from a different data frame ("data_opt").

##############################
# load in necessary packages #
##############################

library(tidyr)
library(reshape2)
library(ggplot2)
library(plyr)
library(pheatmap)
library(dplyr)
library(Hmisc)
options(scipen=999, digits=4)

#########################
# read in and prep data #
#########################

# read in data, setwd or provide path to directory where files are located on line 31 !!
#
# necessary files for this script are:
# Output_stress_duration_optimization.csv
# Output_genotype_survey.csv
# Supplemental Table 1 (I downloaded to desired directory and renamed as SupplementalTable1.csv)

setwd("S:/Groups/LAB-springer/5. Papers in progress/Enders et al Phingerprint/Scripts")

# read in combined csv data and information for annotations in heatmap
data <- read.delim(file="Output_genotype_survey.csv", head=T, sep=",")
data_opt <- read.delim(file="Output_stress_duration_optimization.csv", head=T, sep=",")
geno_info <- read.delim(file="SupplementalTable1.csv", head=T, sep=",")

# remove "Treatment_", make factor
data$treatment <- gsub("Treatment_", "", data$treatment)
data$treatment <- factor(data$treatment, levels=c("Control", "Cold"))

# make genotypes factors
data$genotype <- factor(data$genotype, levels=unique(data$genotype))

# correct switched pots
data$plant <- as.character(data$plant)
data[(data$plot==2325 & data$day==15 & data$plant=="B" & data$area==37250),][, "plant"] <- "C"
data[(data$plot==2325 & data$day==15 & data$plant=="C" & data$area==34544),][, "plant"] <- "B"
data[(data$plot==2376 & data$day==9 & data$plant=="A" & data$area==25150),][, "plant"] <- "B"
data[(data$plot==2376 & data$day==9 & data$plant=="B" & data$area==24466),][, "plant"] <- "A"
data[(data$plot==2446 & data$day==9 & data$plant=="A" & data$area==25315),][, "plant"] <- "B"
data[(data$plot==2446 & data$day==9 & data$plant=="B" & data$area==14974),][, "plant"] <- "A"
data[(data$plot==2525 & data$day==12 & data$plant=="A" & data$area==27599),][, "plant"] <- "B"
data[(data$plot==2525 & data$day==12 & data$plant=="B" & data$area==26520),][, "plant"] <- "C"
data[(data$plot==2525 & data$day==12 & data$plant=="C" & data$area==25549),][, "plant"] <- "A"
data[(data$plot==2778 & data$day==16 & data$plant=="B" & data$area==75522),][, "plant"] <- "C"
data[(data$plot==2778 & data$day==16 & data$plant=="C" & data$area==76765),][, "plant"] <- "B"

# split data into different categories
morph_data <- data[, c(1:17, 19:24)]

# remove data
rm(data)

# make long format
long_data <- gather(morph_data, key=trait, value=value, 7:23)
long_data$experiment <- factor(long_data$experiment, levels=c("78", "79", "80", "81", "83", "84", "85", "95", "99", "100", "102", "104", "106", "111"))

rm(morph_data)

long_data$trait <- gsub("_", ".", long_data$trait)
long_data <- long_data[order(long_data$genotype), ]

#######################################
# Supplemental Table 2 - sample sizes #
#######################################

n1 <- ddply(long_data[long_data$trait=="area",c(2,3,4,7,8)],
            .(genotype, day, treatment, trait),
            summarise, n=length(value))

n1 <- spread(n1, key=day, value=n)

write.csv(n1[,-3], "SupplementalTable2.csv", row.names=F)

rm(n1)

#######################
# Figure 5 - necrosis #
#######################

necrosis_data <- long_data[long_data$trait=="percent.necrosis",]
necrosis_13 <- necrosis_data[necrosis_data$day==13,]
necrosis_anova <- aov(value ~ genotype + treatment + genotype:treatment, data=necrosis_13)
necrosis_Tukey <- TukeyHSD(necrosis_anova)
tuk_geno <- as.data.frame(necrosis_Tukey$genotype)
tuk_treat <- as.data.frame(necrosis_Tukey$treatment)
tuk_genotreat <- as.data.frame(necrosis_Tukey$`genotype:treatment`)

temp <- as.data.frame(row.names(tuk_genotreat))
temp[,2:3] <- colsplit(temp$`row.names(tuk_genotreat)`, "-", names=c("one", "two"))
temp[,4:5] <- colsplit(temp$one, ":", names=c("genotype","four"))
temp[,6:7] <- colsplit(temp$two, ":", names=c("five","six"))

tuk_genotreat[, 5:8] <- temp[,c(4:7)]
# want rows where genotype is the same
tuk_genotreat <- tuk_genotreat[tuk_genotreat$genotype==tuk_genotreat$five,]


necrosis_avg <- ddply(necrosis_data, .(genotype, day, treatment), summarise, mean=mean(value), sd=sd(value), 
                      se=sd(value)/sqrt(length(value)), n=length(value))

combo <- necrosis_avg[necrosis_avg$day==13,]
combo$pvalue <- NA

for (i in unique(combo$genotype)){
  combo[which(combo$genotype==i), 8] <- tuk_genotreat[which(tuk_genotreat$genotype==i), 4] 
}

# order based on signifance
yes <- combo[combo$pvalue<=0.05,]
yes_cold <- yes[yes$treatment=="Cold",]
yes_cold <- yes_cold[order(yes_cold$mean),]

yes_order <- yes_cold$genotype

no <- combo[combo$pvalue>0.05,]
no_cold <- no[no$treatment=="Cold",]
no_cold <- no_cold[order(no_cold$mean),]

no_order <- no_cold$genotype

new_order <- c(as.character(no_order), as.character(yes_order))
combo$genotype <- factor(combo$genotype, levels=new_order)

pdf('Figure5B.pdf', width=6, height=2, onefile=F)
ggplot(combo, aes(x=genotype, y=mean, color=treatment, shape=treatment))+
  geom_rect(aes(xmin="Hp301", 
                xmax="Ki11", 
                ymin=-Inf, ymax=Inf), show.legend=FALSE, fill="gray88", alpha=0.1, color=NA)+
  geom_point(size=1)+
  expand_limits(y=0)+
  scale_y_continuous(labels=function(x)x*100)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.3, size=0.3)+
  theme_bw()+
  theme(panel.grid=element_blank(),
        plot.title=element_text(size=6, colour="black", hjust=0.5), 
        text=element_text(size=6, colour="black"),
        strip.text = element_text(size=6, colour="black", margin = margin(t=0.1, r=0.1, b=0.1, l=0.1, unit = "pt")),
        axis.title = element_text(size=6, colour="black"),
        axis.text.y = element_text(size=6, colour="black"),
        axis.text.x = element_text(size=6, colour="black", angle=45, hjust=1),
        legend.text = element_text(size=6, colour="black"),
        legend.title = element_blank(),
        legend.position=c(0.1, 0.75),
        legend.margin=margin(t=1, r=1, b=1, l=1, unit="mm"),
        legend.box.background = element_rect(colour = "black"))+
  scale_color_manual(values = c("black","steelblue3"))+
  ylab("Mean % Necrosis at 13 DAS")+
  xlab("Genotype")
dev.off()

# example genotypes over time

necrosis_avg <- necrosis_avg[(necrosis_avg$genotype %in% c("B73", "MoG", "Ki11", "M162W", "NC350", "F2")),]
necrosis_avg$genotype <- factor(necrosis_avg$genotype, levels=c("Ki11", "MoG", "F2", "B73", "M162W", "NC350")) 

pdf('Figure5A.pdf', width=6, height=1.5, onefile=F)
ggplot(necrosis_avg, 
       aes(x=day, y=mean, color=treatment, shape=treatment))+
  geom_rect(aes(xmin=9, xmax=11, ymin=-Inf, ymax=Inf), show.legend=FALSE, fill="#c6dbef", alpha=0.1, color=NA)+
  geom_point(size=1)+
  geom_line(size=0.3)+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se),width=0.3, size=0.3)+
  expand_limits(y=0)+
  theme_bw()+
  theme(panel.grid=element_blank(),
        plot.title=element_text(size=6, colour="black", hjust=0.5), 
        text=element_text(size=6, colour="black"),
        strip.text = element_text(size=6, colour="black", margin = margin(t=0.1, r=0.1, b=0.1, l=0.1, unit = "pt")),
        axis.title = element_text(size=6, colour="black"),
        axis.text = element_text(size=6, colour="black"),
        legend.text = element_text(size=6, colour="black"),
        legend.title = element_blank(),
        legend.position = c(0.92,0.65),
        legend.margin=margin(t=1, r=1, b=1, l=1, unit="mm"),
        legend.box.background = element_rect(colour = "black"))+
  scale_color_manual(values = c("black","steelblue3"))+
  scale_y_continuous(labels=function(x)x*100)+
  ylab("Mean % Necrosis")+
  xlab("Days After Sowing (DAS)")+
  scale_x_continuous(breaks=c(8,10,12,14,16))+
  facet_wrap(~genotype, nrow=1)
dev.off()

rm(combo, i, necrosis_13, necrosis_anova, necrosis_avg, necrosis_data, necrosis_Tukey, new_order, no, no_cold, temp, tuk_geno,
   tuk_genotreat, tuk_treat, yes, yes_cold, yes_order, no_order)

########################################
# line plots of means for supplemental #
########################################

mean_data <- ddply(long_data, .(genotype, trait, treatment, day), summarise,
                   mean_value=mean(value), sd_value=sd(value),n=length(value),
                   se_value=sd(value)/sqrt(length(value)))

mean_data$genotype <- as.character(mean_data$genotype)
mean_data <- mean_data[order(mean_data$genotype),]

for(i in unique(mean_data$trait)){
  temp <- mean_data[mean_data$trait==i, ]
  pdf(paste("SuppLineplot_", i, ".pdf", sep=""), width=6, height=8, onefile =F)
  print(ggplot(temp, aes(x=day, y=mean_value, color=treatment, shape=treatment))+
          geom_rect(aes(xmin=9, xmax=11, ymin=-Inf, ymax=Inf), show.legend=FALSE, fill="#c6dbef", alpha=0.1, color=NA)+
          geom_point(size=0.7)+
          geom_line(size=0.3)+
          geom_errorbar(aes(ymin=mean_value-se_value, ymax=mean_value+se_value), width=0.25, size=0.3)+
          facet_wrap(~genotype, ncol=5)+
          ggtitle(i)+
          expand_limits(y=0)+
          theme_bw()+
          theme(panel.grid=element_blank(),
                plot.title=element_text(size=6, colour="black", hjust=0.5), 
                text=element_text(size=6, colour="black"),
                strip.text = element_text(size=6, colour="black", margin = margin(t=0.1, r=0.1, b=0.1, l=0.1, unit = "pt")),
                axis.title = element_text(size=6, colour="black"),
                axis.text = element_text(size=6, colour="black"),
                legend.text = element_text(size=6, colour="black"),
                legend.title = element_text(size=6, colour="black"))+
          scale_color_manual(values = c("black", "steelblue3"))+
          xlab("Days After Sowing (DAS)")+
          ylab("Mean Value (pixels)")+
          scale_x_continuous(breaks=c(8,10,12,14,16))
  )
  dev.off()
}

rm(i, temp, mean_data)

####################################
# Figure 6 - plots with time shift #
####################################

mean_data <- ddply(long_data, .(genotype, trait, treatment, day), summarise,
                   mean_value=mean(value), sd_value=sd(value),n=length(value),
                   se_value=sd(value)/sqrt(length(value)))

ex.plot <- mean_data
c2 <- ex.plot[ex.plot$treatment=="Cold",]
c2$day <- c2$day-2 
c2$treatment <- "Cold-2d"
ex.plot <- rbind(ex.plot, c2)

ex.plot$label <- NA
ex.plot[ex.plot$treatment=="Control",9] <- ex.plot[ex.plot$treatment=="Control",4]
ex.plot[ex.plot$treatment=="Cold",9] <- ex.plot[ex.plot$treatment=="Cold",4]
ex.plot[ex.plot$treatment=="Cold-2d",9] <- ex.plot[ex.plot$treatment=="Cold-2d",4]+2
ex.plot <- ex.plot[!ex.plot$day %in% c(6,7),]

contcold <- ex.plot[!ex.plot$treatment=="Cold-2d",]

pdf(paste("Figure6A.pdf", sep=""), width=6.5, height=1.3, onefile = F)
ggplot(contcold[(contcold$trait=="area" & contcold$genotype %in% c("B73", "Mo17", "Ki11", "Tzi8", 
                                                                   "PH207", "MoG", "P39")),], 
       aes(x=day, y=mean_value, color=treatment, fill=treatment,label=label))+
  geom_rect(aes(xmin=9, xmax=11, ymin=-Inf, ymax=Inf), show.legend=FALSE, fill="#c6dbef", alpha=0.1, color=NA)+
  geom_errorbar(aes(ymin=mean_value-se_value, ymax=mean_value+se_value), width=0.3, size=0.3)+
  geom_line(size=0.3)+
  geom_point(size=1.5, shape=21, fill="white", color="white")+
  geom_text(size=1.5, aes(color=treatment))+
  facet_wrap(~genotype, nrow=1)+
  expand_limits(y=0)+
  theme_bw()+
  theme(panel.grid=element_blank(),
        plot.title=element_text(size=6, colour="black", hjust=0.5), 
        text=element_text(size=6, colour="black"),
        strip.text = element_text(size=6, colour="black", margin = margin(t=0.1, r=0.1, b=0.1, l=0.1, unit = "pt")),
        axis.title = element_text(size=6, colour="black"),
        axis.text = element_text(size=6, colour="black"),
        legend.text = element_text(size=6, colour="black"),
        legend.title = element_blank())+
  scale_color_manual(values = c("black","steelblue3"))+
  xlab("Days After Sowing (DAS)")+
  ylab("Mean Area (pixels)")+
  xlim(7.9,16.1)
dev.off()

contcold2 <- ex.plot[!ex.plot$treatment=="Cold",]

pdf(paste("Figure6B.pdf", sep=""), width=6.5, height=1, onefile = F)
ggplot(contcold2[(contcold2$trait=="area" & contcold2$genotype %in% c("B73", "Mo17", "Ki11", "Tzi8", 
                                                                      "PH207", "MoG", "P39")),], 
       aes(x=day, y=mean_value, color=treatment, fill=treatment,label=label))+
  geom_errorbar(aes(ymin=mean_value-se_value, ymax=mean_value+se_value), width=0.3, size=0.3)+
  geom_line(size=0.3)+
  geom_point(size=1.5, shape=21, fill="white", color="white")+
  geom_text(size=1.5, aes(color=treatment))+
  facet_wrap(~genotype, nrow=1)+
  expand_limits(y=0)+
  theme_bw()+
  theme(panel.grid=element_blank(),
        plot.title=element_text(size=6, colour="black", hjust=0.5), 
        text=element_text(size=6, colour="black"),
        strip.text = element_text(size=6, colour="black", margin = margin(t=0.1, r=0.1, b=0.1, l=0.1, unit = "pt")),
        axis.title = element_text(size=6, colour="black"),
        axis.text = element_text(size=6, colour="black"),
        legend.text = element_text(size=6, colour="black"),
        legend.title = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_color_manual(values = c("black","red"))+
  xlab("Days After Sowing (DAS)")+
  ylab("Mean Area (pixels)")+
  xlim(7.9,16.1)
dev.off()

pdf(paste("Figure6C.pdf", sep=""), width=2.5, height=1.3, onefile=F)
ggplot(contcold2[(contcold2$trait=="area" & contcold2$genotype %in% c("Ki11")),], 
       aes(x=day, y=mean_value, color=treatment, fill=treatment,label=label))+
  geom_errorbar(aes(ymin=mean_value-se_value, ymax=mean_value+se_value), width=0.3, size=0.3)+
  geom_line(size=0.3)+
  geom_point(size=1.5, shape=21, fill="white", color="white")+
  geom_text(size=1.5, aes(color=treatment))+
  facet_wrap(~genotype, nrow=1)+
  expand_limits(y=0)+
  theme_bw()+
  theme(panel.grid=element_blank(),
        plot.title=element_text(size=6, colour="black", hjust=0.5), 
        text=element_text(size=6, colour="black"),
        strip.text = element_text(size=6, colour="black", margin = margin(t=0.1, r=0.1, b=0.1, l=0.1, unit = "pt")),
        axis.title = element_text(size=6, colour="black"),
        axis.text = element_text(size=6, colour="black"),
        legend.text = element_text(size=6, colour="black"),
        legend.title = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  scale_color_manual(values = c("black","red"))+
  xlab("Days After Sowing (DAS)")+
  ylab("Mean Area (pixels)")+
  xlim(7.9,16.1)
dev.off()

rm(c2, contcold, contcold2, ex.plot, mean_data)

######################
# Figure 7 - heatmap #
######################

long_data2 <- long_data

long_data2$interval <- rep(NA, nrow(long_data2))
long_data2[(long_data2$day==8 & long_data2$treatment=="Control"), ][, "interval"] <- "before.8-9"
long_data2[(long_data2$day==9 & long_data2$treatment=="Control"), ][, "interval"] <- "before.8-9"
long_data2[(long_data2$day==10 & long_data2$treatment=="Control"), ][, "interval"] <- "cont.9-11_cold.11-13"
long_data2[(long_data2$day==11 & long_data2$treatment=="Control"), ][, "interval"] <- "cont.9-11_cold.11-13"
long_data2[(long_data2$day==12 & long_data2$treatment=="Control"), ][, "interval"] <- "cont.11-14_cold.13-16"
long_data2[(long_data2$day==13 & long_data2$treatment=="Control"), ][, "interval"] <- "cont.11-14_cold.13-16"
long_data2[(long_data2$day==14 & long_data2$treatment=="Control"), ][, "interval"] <- "cont.11-14_cold.13-16"

long_data2[(long_data2$day==8 & long_data2$treatment=="Cold"), ][, "interval"] <- "before.8-9"
long_data2[(long_data2$day==9 & long_data2$treatment=="Cold"), ][, "interval"] <- "before.8-9"
long_data2[(long_data2$day==10 & long_data2$treatment=="Cold"), ][, "interval"] <- "during.9-11"
long_data2[(long_data2$day==11 & long_data2$treatment=="Cold"), ][, "interval"] <- "cont.9-11_cold.11-13"
long_data2[(long_data2$day==12 & long_data2$treatment=="Cold"), ][, "interval"] <- "cont.9-11_cold.11-13"
long_data2[(long_data2$day==13 & long_data2$treatment=="Cold"), ][, "interval"] <- "cont.9-11_cold.11-13"
long_data2[(long_data2$day==14 & long_data2$treatment=="Cold"), ][, "interval"] <- "cont.11-14_cold.13-16"
long_data2[(long_data2$day==15 & long_data2$treatment=="Cold"), ][, "interval"] <- "cont.11-14_cold.13-16"
long_data2[(long_data2$day==16 & long_data2$treatment=="Cold"), ][, "interval"] <- "cont.11-14_cold.13-16"

temp <- rbind(long_data2[(long_data2$day==9 & long_data2$treatment=="Control"), ],
              long_data2[(long_data2$day==11 & long_data2$treatment=="Control"), ],
              long_data2[(long_data2$day==9 & long_data2$treatment=="Cold"), ],
              long_data2[(long_data2$day==11 & long_data2$treatment=="Cold"), ],
              long_data2[(long_data2$day==13 & long_data2$treatment=="Cold"), ])
temp[(temp$day=="9" & temp$treatment=="Control" & temp$interval=="before.8-9"), 9] <- "cont.9-11_cold.11-13"
temp[(temp$day=="11" & temp$treatment=="Control" & temp$interval=="cont.9-11_cold.11-13"), 9] <- "cont.11-14_cold.13-16"
temp[(temp$day=="9" & temp$treatment=="Cold" & temp$interval=="before.8-9"), 9] <- "during.9-11"
temp[(temp$day=="11" & temp$treatment=="Cold" & temp$interval=="cont.9-11_cold.11-13"), 9] <- "during.9-11"
temp[(temp$day=="13" & temp$treatment=="Cold" & temp$interval=="cont.9-11_cold.11-13"), 9] <- "cont.11-14_cold.13-16"

long_data2 <- rbind(long_data2, temp)

long_data2$day <- as.numeric(as.character(long_data2$day))
long_data2 <- long_data2[order(long_data2$day),]
long_data2$plot_plant <- paste(long_data2$plot, long_data2$plant, sep="_")

long_data2$grouping <- paste(long_data2$genotype, long_data2$treatment, long_data2$trait, long_data2$plot_plant,
                             long_data2$interval, sep="_")

fit <- dlply(long_data2, .(grouping), lm, formula=value~day)
slopes <- ldply(fit, coef)
colnames(slopes) <- c("grouping", "intercept", "slope")
slopes[, 4:9] <- colsplit(slopes$grouping,"_", names=c("genotype", "treatment","trait", 
                                                       "plot", "plant","interval"))

slopes$slope <- round(slopes$slope, 2)
slopes <- slopes[, -2]
slopes <- slopes[, -1]

slopes <- slopes[slopes$interval %in% c("cont.9-11_cold.11-13", "cont.11-14_cold.13-16"),]

slopes <- na.omit(slopes) #83718 to 81648

slopes <- ddply(slopes, .(genotype, trait, treatment, interval), summarise,
                slope=mean(slope))

wide_avg <- spread(slopes, key=treatment, value=slope)
wide_avg$ratio <- log2(wide_avg$Cold/wide_avg$Control)

wide_avg <- wide_avg[!wide_avg$trait=="percent.necrosis",]
wide_avg$interval <- factor(wide_avg$interval, levels=c("cont.9-11_cold.11-13", "cont.11-14_cold.13-16"))

stuff <- wide_avg
brown_data <- long_data[long_data$trait=="percent.necrosis",]
brown_avg <- ddply(brown_data, .(genotype, day, treatment), summarise, mean=mean(value), sd=sd(value))
brown_avg <- brown_avg[brown_avg$day=="13",]
brown_avg <- spread(brown_avg[,-5], key=treatment, value=mean)
brown_avg$ratio <- log2(brown_avg$Cold/brown_avg$Control)
brown_avg$trait <- "percent_necrosis_d13"
colnames(brown_avg) <- c("genotype", "interval", "Cold", "Control", "ratio", "trait")
brown_avg$interval <- factor(brown_avg$interval, levels=unique(brown_avg$interval))

stuff <- rbind(wide_avg[,c(1,2,3,6)], brown_avg[,c(1,2,5,6)])

stuff <- spread(stuff, key=genotype, value=ratio)
row_info <- stuff[,1:2]
rownames(row_info) <- paste(stuff$trait, stuff$interval, sep="_")
rownames(stuff) <- rownames(row_info)

col_info <- as.data.frame(colnames(stuff[,3:42]))
colnames(col_info) <- c("genotype")
rownames(col_info) <- colnames(stuff[,3:42])

col_info <- left_join(col_info, geno_info)
rownames(col_info) <- colnames(stuff[,3:42])

inbred <- stuff[(stuff$trait %in% c("percent_necrosis_d13", "width","height","area","perimeter","hull.area")), ]

col_info <- col_info[,c(4,5,15)]
col_info$latitude <- as.numeric(as.character(col_info$latitude))
col_info$kernel.texture <- factor(col_info$kernel.texture, levels=c("Dent", "Flint", "Flint-Dent Mixed",
                                                                    "Popcorn/Sweetcorn"))

kernel.texture <- c("#ef3b2c","#fed976", "#fd8d3c","#800026")
names(kernel.texture) <- levels(col_info$kernel.texture)

col_info$population.group <- factor(col_info$population.group, levels=c("NSS","SS","Iodent","NSS/Iodent","SS/Iodent",
                                                                        "TS","Popcorn/Sweetcorn", "Mixed"))

population.group <- c("#80cdc1","#f768a1","#fc8d59","#7fbc41","#1d91c0","#c994c7","#800026","#737373")
names(population.group) <- levels(col_info$population.group)

# latitude
latitude <- colorRampPalette(c('#d9d9d9','#252525'))(20)
latitude <- colorRampPalette(c('#d9f0a3','#004529'))(20)
names(latitude) <- sort(unique(col_info$latitude))

anno_colors <- list(kernel.texture=kernel.texture,
                    population.group=population.group, latitude=latitude)

rownames(inbred) <- c("area_early_recovery","area_late_recovery","height_early_recovery",
                      "height_late_recovery","hull.area_early_recovery","hull.area_late_recovery",
                      "percent_necrosis_d13","perimeter_early_recovery" ,"perimeter_late_recovery",
                      "width_early_recovery","width_late_recovery")

breaks <- c(seq(-2, 2, length.out = 50))
my_palette <- colorRampPalette(c('#2c7bb6','#ffffbf','#d7191c'))(49)

breaks <- c(seq(-2, 2, length.out = 50))
breaks2 <- breaks
breaks2[length(breaks)] <- max(apply(inbred[,3:42],2,function(x) max(as.numeric(na.omit(x)))), max(breaks))
breaks2[1] <- min(apply(inbred[,3:42],2,function(x) min(as.numeric(na.omit(x)))),min(breaks))

pdf('Figure7.pdf', width=8, height=8, onefile=F)
pheatmap(inbred[, c(3:42)],
         annotation_col = col_info,
         cluster_rows = T,
         cluster_cols = T,
         color=my_palette,
         breaks = breaks2,
         legend_breaks = c(-2,0,2),
         legend_labels = c(-2,0,2),
         annotation_colors = anno_colors,
         cellwidth = 6,
         cellheight= 6,
         fontsize=6)
dev.off()

### revisions---------------------------------------------------------------------------------

# correlations of traits (log2(FC) values from heatmap)

temp <- inbred
temp$combined <- rownames(inbred)
temp$combined <- gsub("_recovery", "", temp$combined)

temp <- gather(temp, key=genotype, value=value, 3:42)
temp <- temp[,3:5]
temp <- spread(temp, key=combined, value=value)

# pearson
temp_corr <- rcorr(as.matrix(temp[,2:12]), type="pearson")

pdf('SupplementalFigure8.pdf', width=4.5, height=4, onefile = F)
pheatmap(temp_corr$r,
         color=colorRampPalette(c('#d7191c','#ffffbf','#2c7bb6'))(49),
         breaks = c(seq(-1, 1, length=50)),
         legend_breaks = c(-1,-0.5,0,0.5,1),
         legend_labels = c(-1,-0.5,0,0.5,1),
         display_numbers = matrix(ifelse(temp_corr$P<=0.05 & temp_corr$P>0.01, "*",
                                         ifelse(temp_corr$P<=0.01 & temp_corr$P>0.001, "**",
                                                ifelse(temp_corr$P<=0.001, "***",""))), nrow(temp_corr$P)),
         cellwidth = 9,
         cellheight= 9,
         fontsize=6,
         fontsize_number = 6)
dev.off()

# Normalize, then cluster.
temp2 <- as.data.frame(apply(inbred[, 3:42], 1, 
                             function(x)((x-min(x))/(max(x)-min(x)))))

temp2 <- t(temp2)

png("heatmap_normalized_01.png", width=8, height=8, units="in", res=600)
pheatmap(temp2,
         annotation_col = col_info,
         cluster_rows = T,
         cluster_cols = T,
         color=colorRampPalette(c('#ffffbf','#2c7bb6'))(49),
         legend_breaks = c(0,1),
         legend_labels = c(0,1),
         annotation_colors = anno_colors,
         cellwidth = 6,
         cellheight= 6,
         fontsize=6)
dev.off()

# end of new figures for revisions ----------------------------------------------------------

rm(anno_colors, breaks, breaks2, brown_avg, brown_data, col_info, fit, geno_info, inbred, kernel.texture, latitude, long_data, long_data2, 
   my_palette, population.group, row_info, slopes, stuff, temp, wide_avg)

###########################################
# Figure 3 - stress duration optimization #
###########################################

# prep data from Output_stress_duration_optimization.csv

data_opt <- data_opt[data_opt$experiment==108, ]

data_opt$treatment <- gsub("Treatment_", "", data_opt$treatment)
data_opt$treatment <- factor(data_opt$treatment, levels=c("Control", "1day_cold", "2day_cold", "3day_cold", "4day_cold"))

data_opt$genotype <- factor(data_opt$genotype, levels=unique(data_opt$genotype))

# split data into different categories
morph_data <- data_opt[, c(1:17, 19:24)]

# make long format
long_data <- gather(morph_data, key=trait, value=value, 7:23)

long_data$trait <- gsub("_", ".", long_data$trait)

# average and plot
mean_data <- ddply(long_data, .(genotype, trait, treatment, day), summarise,
                   mean_value=mean(value), sd_value=sd(value),n=length(value),
                   se_value=sd(value)/sqrt(length(value)))

temp <- mean_data[(mean_data$trait %in% c("area", "height", "width", "percent.necrosis")), ]
temp$trait <- factor(temp$trait, levels=c("area", "height", "width", "percent.necrosis"))

pdf(paste("Figure3.pdf", sep=""), width=3, height=5, onefile=F)
ggplot(temp, 
       aes(x=day, y=mean_value, color=treatment, shape=treatment))+
  geom_errorbar(aes(ymin=mean_value-se_value, ymax=mean_value+se_value), width=0.3, size=0.3)+
  geom_line(size=0.3)+
  geom_point(size=1)+
  facet_grid(trait~genotype, scales = "free_y")+
  #  expand_limits(y=0)+
  theme_bw()+
  theme(panel.grid=element_blank(),
        plot.title=element_text(size=6, colour="black", hjust=0.5), 
        text=element_text(size=6, colour="black"),
        strip.text = element_text(size=6, colour="black", margin = margin(t=0.1, r=0.1, b=0.1, l=0.1, unit = "pt")),
        axis.title = element_text(size=6, colour="black"),
        axis.text = element_text(size=6, colour="black"),
        legend.text = element_text(size=6, colour="black"),
        legend.title = element_blank(),
        legend.position = "top")+
  guides(color = guide_legend(label.position = "top"))+
  scale_color_manual(values = c("black", "steelblue1","steelblue3", "steelblue4", "midnightblue"), 
                     breaks = c("Control", "1day_cold", "2day_cold", "3day_cold", "4day_cold"),
                     labels = c("Control", "1 d", "2 d", "3 d", "4 d"))+
  scale_shape_manual(values = c(15, 16, 17, 18, 4), 
                     breaks = c("Control", "1day_cold", "2day_cold", "3day_cold", "4day_cold"),
                     labels = c("Control", "1 d", "2 d", "3 d", "4 d"))+
  xlab("Days After Sowing (DAS)")+
  ylab("Mean Value (pixels)")+
  scale_x_continuous(breaks=c(8,10,12,14,16))
dev.off()
