
cm <- read.xlsx("./Input/countmatrix202412/counts.distri.Pf.all.xlsx") ###Double check if includes total columns

df <- as.data.frame(table(cm$Total))
df2 <- df[c(1:20),]
histo.freq.figure <- ggplot(df2, aes(x=Var1, y=Freq)) + 
  xlab('Number of reads mapped to UII') + 
  geom_bar(stat='identity') +
  ylab('Count')+
  ggtitle("")
print(histo.freq.figure)


cm_sense 
cm_antisense 
