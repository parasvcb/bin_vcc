library(ggplot2)
library(scales)
res1='../../../../plotting/difference_finaledges.csv'
res2='../../../../plotting/raw_distribution.csv'
df1=read.csv(res1)
df2=read.csv(res2)
print (head(df1))
gg <- ggplot(data=df1,aes(x=values)) 
#gg <- gg + stat_bin(aes(y=..density..), breaks = seq(min(df1$values), max(df1$values), by = .1), color="white")
#gg <- gg + geom_density(aes(y=..density..))
gg <- gg + stat_bin(aes(y=..count../sum(..count..)), breaks = seq(-1,1, by = .1), colour="black", fill="grey")
gg <- gg + geom_density(aes(x=values,..scaled..),size=0.5)
gg <- gg + scale_x_continuous('cij value', breaks = seq(-1,1, by = .1))
gg <- gg + scale_y_continuous('Frequency', breaks = seq(0,1, by = .2))
gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
ggsave(filename = '../../../../plotting/FigR2.pdf')

print (head(df2))

gg <- ggplot() 
gg <- gg + geom_density(data=df2,aes(x=cij,..scaled..,color=group),size=0.5) + geom_vline(xintercept=0.6, color='black',size=1)
gg <- gg + scale_color_manual(values=c("blue","cyan","skyblue","red","darkorange","maroon"))
gg <- gg + scale_x_continuous('cij value', breaks = seq(-0.2,1, by = .1))
gg <- gg + scale_y_continuous('Frequency', breaks = seq(0,1, by = .2))
gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
ggsave(filename = '../../../../plotting/FigR1.pdf')


#follwing is the new addition
library(ggplot2)

res3='../../../../plotting/raw_distribution_pythonpseudoHist_with0bin.csv'
df3=read.csv(res3)
gg <- ggplot() 
gg <- gg + geom_line(data=df3,aes(x=bins,y=frequency, color=reptype),size=0.5) + geom_vline(xintercept=0.6, color='black',size=1)
gg <- gg + geom_point(data=df3,aes(x=bins,y=frequency, color=reptype),size=0.5)
gg <- gg + scale_color_manual(values=c("blue","cyan","skyblue","red","darkorange","maroon"))
gg <- gg + scale_x_continuous('cij value', breaks = seq(-0.2,1, by = .1))
gg <- gg + scale_y_continuous('Frequency', breaks = seq(0,1, by = .2))
gg <- gg + xlab('Bins')
gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
ggsave(filename = '../../../../plotting/FigR1_with0bin.pdf')

res3='../../../../plotting/raw_distribution_pythonpseudoHist_without0bin.csv'
df3=read.csv(res3)
gg <- ggplot() 
gg <- gg + geom_line(data=df3,aes(x=bins,y=frequency, color=reptype),size=0.5) + geom_vline(xintercept=0.6, color='black',size=1)
gg <- gg + geom_point(data=df3,aes(x=bins,y=frequency, color=reptype),size=0.5)
gg <- gg + scale_color_manual(values=c("blue","cyan","skyblue","red","darkorange","maroon"))
gg <- gg + scale_x_continuous('cij value', breaks = seq(-0.2,1, by = .1))
gg <- gg + scale_y_continuous('Frequency', breaks = seq(0,1, by = .2))
gg <- gg + xlab('Bins')
gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 10, angle = 45),axis.text.y = element_text( hjust = 1, size = 10), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
ggsave(filename = '../../../../plotting/FigR1_without0bin.pdf')


   

# gg <- ggplot(data=df1,aes(x=hbx,y=mean,fill=hbtype))
# gg <- gg + geom_bar(stat="identity",position=position_dodge())
# gg <- gg + geom_text(aes(label=mean),vjust=1, color="black",position = position_dodge(0.9), size=3)
# gg <- gg + geom_errorbar(data= dfgen, aes(x=hbx,y=mean,ymin=mean-std, ymax=mean+std),position = position_dodge(0.9),size=0.3,width=0.5)
# gg <- gg + theme (axis.text.x = element_text( hjust = 1, size = 7, angle = 45),axis.text.y = element_text( hjust = 1, size = 7), legend.position="top", panel.background = element_rect(fill = "white", colour = "grey50"),panel.grid.major = element_line(colour = "grey90"),panel.grid.minor = element_line(colour = "grey95",size = 0.25))
# gg <- gg + facet_wrap(~polymer,ncol=1,scale='free_x') 
# ggsave(filename = outputfile, height=15, width=5)
