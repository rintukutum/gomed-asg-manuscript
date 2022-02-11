rm(list=ls())
#---------
# load packages
library("RColorBrewer")
library("readr")
library("ggplot2")
library("ggpubr")
library("reshape2")
library("plyr")
#---------
# Figure 1A
fig01.data <- read_csv(
  'data/phenotype-sca-count.csv'
)
fig01.data$perc <- round(
  fig01.data$count/sum(fig01.data$count),
  3)*100
count.log10 <- log10(fig01.data$count)
names(count.log10) <- fig01.data$phenotype
png('figures/fig-1A.png',
    width = 800,height = 800,res = 150)
par(mar=c(5,10,4,2)+1)
br <- barplot(
  count.log10,las=2,horiz = TRUE,
  col=rev(c('grey50',
            brewer.pal(7,'Oranges')[4],
            brewer.pal(7,'Blues')
  )),
  xaxt = "n"
)
text(y=br[,1],x=0.05,
     labels = paste0(
       fig01.data$count,
       '(',
       fig01.data$perc, '%)'
     ),
     adj = 0,
     cex = 0.75,
     col = c('white','white',rep('black',6),
             'white'))
axis(1,at=c(0,1,2,3),labels=c(1,10,100,1000))
mtext(side = 2,line = 6,
      text = 'Ataxia subtypes',
      cex = 1.1
)
mtext(side = 1,line = 2,
      text = '# Samples',
      cex = 1.1
)
dev.off()
#---------------
rm(list=ls())
#-------------
# Figure 1B
fig01b.data <- read_csv('data/site-diagnosis.csv')
colnames(fig01b.data)[2:3] <- c('Undiagnosed','Diagnosed')
site.diag <- melt(fig01b.data[,1:3])
site.diag <- ddply(site.diag,'Site',function(x){
  x$perc <- round((x$value/sum(x$value))*100,2)
  x
})
sites <- paste0(
  fig01b.data$Site,"\n(N=",
  fig01b.data$`total samples`,")"
)
names(sites) <- fig01b.data$Site

png('figures/fig-1B.png',width = 700,height = 900,res = 150)
ggplot(site.diag,aes(Site,perc,fill=variable)) +
  geom_bar(stat='identity') +
  geom_text(position = position_stack(vjust = 0.5),
            aes(label=
                  paste0(perc,' %')),
            angle=90,
            size=4,
            col='white'
  ) +
  scale_fill_manual(name='',values =c('grey60',"#238B45")) +
  scale_x_discrete(
    limit = names(sites),
    labels = sites
  ) +
  theme_pubr() +
  ylab('Number of cases (in %)') +
  xlab('')
dev.off()
#---------------
rm(list=ls())
#---------------
# Figure 2
fig02.data <- read_csv(
  'data/sca-type-allele.csv'
)
fig02.data$var <- paste0(
  fig02.data$Status,'-',
  fig02.data$Allele)
fig02.data$type <- factor(
  fig02.data$type,
  levels = unique(fig02.data$type))
figNeg <- data.frame(
  fig02.data[fig02.data$var == "negative-A1/A2",]
)
figNeg.max <- plyr::ddply(
  figNeg,'type',function(x){
    x[which.max(x$value),]
  })
figNeg.max$value <- figNeg.max$value + 1
pdf('figures/fig-02.pdf',
    width = 8,height = 6)
ggplot(fig02.data,aes(x=var,y=value)) +
  geom_jitter(aes(fill=var),
              stat = "identity",
              size=1,
              width = 0.25,
              shape=21,alpha=0.75) +
  geom_hline(data=figNeg.max,
             aes(yintercept=value),
             color = 'red',linetype=2) +
  facet_wrap(facets = 'type',ncol = 7,scales = 'free_y') +
  scale_fill_manual(
    name='',
    labels = c('Negative','A1','A2'),
    values = c("#F46D43","#FFFFBF","#3288BD")) +
  ylab('CAG Length') + xlab('') +
  scale_x_discrete(labels=c('Negative (A1/A2)',
                            'PA1','PA2')) +
  theme_pubr() +
  theme(legend.position = 'none',
        axis.text.x = element_text(
          angle = 45,vjust = 1,hjust = 1)
  )
dev.off()
#----------
# Figure 3
fig03.D1 <- read_csv('data/sca-type-age-small.csv')
fig03.D2 <- read_csv('data/sca-type-age-big.csv')
fig03.data <- rbind(
  data.frame(fig03.D1),
  data.frame(fig03.D2)
)
# format data
types <- paste0('SCA',c(1:3,6,7,12,17))
sca.m <- list()
for(i in 1:length(types)){
  type <- paste0(types[i],'.',c('LA','UA'))
  sca <- melt(
    fig03.data[,c('SAMPLE.ID','AGE',type)],
    id=c('SAMPLE.ID','AGE')) 
  sca$variable <- types[i]
  sca.m[[i]] <- na.omit(sca)
}
sca.df <- ldply(sca.m)
sca.df$age.group <- cut(
  sca.df$AGE,
  breaks = c(0,17,29,39,49,59,69,79,89,99),
  labels = c('0-17','18-29','30-39',
             '40-49','50-59','60-69',
             '70-79','80-89','90+')
)
sca.stats <- ddply(
  sca.df,.variables = 'variable',
  function(x){
    data.frame(table(x[,c('age.group','value')]))
  }
)
idx.ok <- sca.stats$Freq != 0
sca.stats <- sca.stats[idx.ok,]
sca.stats$value <- factor(
  sca.stats$value,
  levels = sort(as.numeric(levels(sca.stats$value))
  ))
sca.stats$variable <- factor(
  sca.stats$variable,
  levels = paste0('SCA',c(1:3,6,7,12,17))
)
figNeg.max$value <- factor(figNeg.max$value)
colnames(figNeg.max)[5] <- 'variable'
png('figures/fig-03.png',
    width = 1400,height = 1400,res = 150)
ggplot(sca.stats,aes(x=value,y=age.group)) +
  geom_tile(aes(fill=log10(Freq)),color='black') +
  geom_vline(data = figNeg.max,
             aes(xintercept = value),
             color = 'red',linetype=2) +
  facet_grid(variable~.,) +
  scale_fill_gradient(
    low = "#FFFFFF",high = "#BD0026") +
  xlab('CAG Length') + ylab('Age') +
  theme_pubr() +
  theme(axis.text.x = element_text(
    angle = 90,vjust = 0.5,hjust = 1,size=8),
    axis.text.y = element_text(size=10))
dev.off()
png('figures/fig-03-without-log.png',
    width = 1400,height = 1400,res = 150)
ggplot(sca.stats,aes(x=value,y=age.group)) +
  geom_tile(aes(fill=Freq),color='black') +
  geom_vline(data = figNeg.max,
             aes(xintercept = value),
             color = 'red',linetype=2) +
  facet_grid(variable~.,) +
  scale_fill_gradient(
    low = "#FFFFCC",high = "#BD0026") +
  xlab('CAG Length') + ylab('Age') +
  theme_pubr() +
  theme(axis.text.x = element_text(
    angle = 90,vjust = 0.5,hjust = 1,size=8),
    axis.text.y = element_text(size=10))
dev.off()