library(ggplot2)
library(viridis)
# data
data <- data.frame(
  Dataset = c("yeoh", "yeoh", "yeoh", "yeoh","yeoh",
              "burczynski", "burczynski", "burczynski", "burczynski", "burczynski",
              "nakayama", "nakayama", "nakayama", "nakayama", "nakayama"),
  Method = rep(c("PCA", "LRLDA", "PLS", 'LOL', 'PCA+PLS'),3),
  Error = c(0.012, 0.024, 0.00801, 0.00801, 0.004,
            0.0397435897435897, 0.126923076923077, 0.0397435897435897,0.0705128205128205, 0.03205128,
            0.292727272727273, 0.333636363636364, 0.284545454545455,0.294545454545455, 0.2381818)
)

#horizontal plot
ggplot(data, aes(x=Error, y=Dataset, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_bw() +
  theme(legend.position="bottom") +
  labs(x="Error", y="Dataset", fill="Method")+
  scale_fill_viridis(discrete = TRUE, option = "D") +  # 使用Viridis的配色方案
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 12),
    legend.position = "bottom"
  )

# vertical plot
ggplot(data, aes(x=Dataset, y=Error, fill=Method)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_bw() +
  theme(legend.position="bottom") +
  labs(x="Dataset", y="Error", fill="Method")+
  scale_fill_viridis(discrete = TRUE, option = "D") +  # 使用Viridis的配色方案
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 12),
    legend.position = "bottom"
  )
