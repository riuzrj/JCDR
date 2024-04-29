library(ggplot2)
data = lol.sims.toep(n=300, d=2000)
X = data$X
Y = data$Y
#scatter plot for data
ggplot(data.frame(X), aes(x=X[,1], y=X[,2], color=Y)) + geom_point(size=2) + theme_bw()+
  xlab("") + ylab("")
