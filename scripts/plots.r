setwd("/Users/cora/git_repos/NetProject/results")

df <- read.csv("data.csv")
df$minutes <- df$time/60
require(ggplot2)

timeplot <- ggplot(data=df, aes(y = minutes, x = method)) + 
  geom_point() + ylab("Time (minutes)") +
  xlab("Method") + ggtitle("Comparing Methods by Time"); timeplot

timeplot_p <- ggplot(data=df[df$method=="parsimony",], aes(y = minutes, x = distance)) + 
  geom_point() + ylab("Time (minutes)") +
  xlab("Robinson-Foulds Distance") + ggtitle("Softwired Parsimony Method"); timeplot_p

timeplot_s <- ggplot(data=df[df$method=="snaq",], aes(y = minutes, x = distance)) + 
  geom_point() + ylab("Time (minutes)") +
  xlab("Robinson-Foulds Distance") + ggtitle("PseudoLikelihood SNaQ Method"); timeplot_s

distplot <- ggplot(data=df, aes(x = method, y = distance)) + 
  geom_point() + ylab("Robinson-Foulds Distance") +
  xlab("Method") + ggtitle("Distance by Method"); distplot
