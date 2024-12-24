
results <- read.csv("/Users/ruijuanzhong/Project/gaus_noise_sun_inK.csv")
library(tidyverse)
head(results)
library(dplyr)
library(tidyr) # 为了使用 unnest 函数

# 假设 desired_r 是已经定义好的
desired_r <- seq(min(results$r), 50, by = 1)

# 对数据进行分组，并应用 spline 函数进行平滑处理
filterd_results <- results %>%
  group_by(sim, alg, r) %>%
  summarize(lhat = mean(lhat)) %>%
  #filter(!(alg %in% c('pca+pls') & !(r %in% rs))) %>%
  #把大于 30 的 r 值过滤掉
  filter(r <= 50)
#   mutate(smoothed = list(spline(rs, lhat, xout = desired_r, method = 'fmm', ties = mean))) %>%
# # 使用 unnest 获得平滑后的数据
#   unnest(c(smoothed))

#用 spline 函数进行平滑处理
# smoothed_results <- filterd_results %>%
#   group_by(sim, alg) %>%
#   summarize(x = list(desired_r), y = list(spline(r, lhat, xout = desired_r,
#                                                  method = 'fmm', ties = mean)$y)) %>%
#   unnest(c(x, y))

# dt_plot <- smoothed_results %>%
#   group_by(sim, alg) %>%
#   mutate(best.lhat = min(y),
#          best.r = list(x[y == best.lhat]),
#          lhat.thresh = best.lhat * 1.05,
#          lhat.beats = y <= lhat.thresh,
#          r.star = best.r)

#用 loess 函数进行平滑处理
smoothed_results <- filterd_results %>%
  group_by(sim, alg) %>%
  do({
    # 对每个组应用 loess 平滑
    smooth_span = 0.8  # 这个值可以调整来控制平滑程度，数值范围是 0 到 1
    data.frame(r = .$r, lhat = predict(loess(lhat ~ r, data = ., span = smooth_span)))
  }) %>%
  ungroup()

dt_plot <- smoothed_results %>%
  group_by(sim, alg) %>%
  mutate(best.lhat = min(lhat),
         best.r = r[lhat == best.lhat],
         lhat.thresh = best.lhat * 1.05,
         lhat.beats = lhat <= lhat.thresh,
         #r.star = min(r[lhat.beats]))
         r.star = best.r,
         sim = recode_factor(sim, "unif_noise" = "gaussian_noise")) #%>%
  #把 LRLDA 改为 rrLDA
#mutate(sim = recode_factor(sim, "uni_nois" = "gaussian_noise"))


#用 spline 函数进行平滑处理
# smoothed_results <- filterd_results %>%
#   group_by(sim, alg) %>%
#   summarize(x = list(desired_r), y = list(spline(r, lhat, xout = desired_r, method = 'fmm', ties = mean)$y)) %>%
#   unnest(c(x, y))
#
# dt_plot <- smoothed_results %>%
#   group_by(sim, alg) %>%
#   mutate(best.lhat = min(y),
#          best.r = list(x[y == best.lhat]),
#          lhat.thresh = best.lhat * 1.05,
#          lhat.beats = y <= lhat.thresh,
#          r.star = min(x[lhat.beats]))



# Continue with your plotting code
# ...

# dt_plot <- results %>%
#   group_by(sim,alg,r) %>%
#   summarise(lhat=mean(lhat)) %>%
#   #filter(!(alg %in% c('pca+pls') & !(r %in% rs))) %>%
#   mutate(best.lhat=min(lhat),
#          best.r=list(r[lhat == best.lhat]),
#          lhat.thresh=best.lhat*1.05,
#          lhat.beats=lhat <= best.lhat,
#          r.star=min(r[lhat.beats]))

#mutate(lhat_adjust = spline(rs, lhat, xout = desired_r, method = 'fmm', ties = mean))

#把 alg 中‘pca+pls’这个类中 r 大于 50 的数据去掉
#filter(!(alg %in% c('pca+pls') & r > 50))
#筛选出alg 中‘pca+pls’这个类中 r等于 rs 的数值



# mutate(alg=factor(recode_factor(alg, "QOL"="QOQ", "LRLDA"="rrLDA"),
#                   levels=algs.name,
#                   ordered=TRUE)) %>%
#filter(!((alg %in% c("RLOL")) & (sim == "Cross"))) %>%
# filter(!(alg %in% c("MVI"))) %>%
#filter(!((alg %in% c("QOQ")) & (sim %in% c("Trunk-3", "Robust")))) %>%
#mutate(
# sim=recode_factor(sim, "Trunk-3"="(A) Trunk-3", "Robust"="(B) Robust", "Cross"="(C) Cross",
# "class.diff-3"="(D) class.diff-3"),
# sim=factor(sim, levels=c("(A) Trunk-3", "(B) Robust", "(C) Cross","(D) class.diff"), ordered=TRUE))
# sim=factor(sim, levels=c("Trunk-3", "Robust", "Cross","class.diff-3"), ordered=TRUE))
#sim=factor(sim, levels=sims.name, ordered=TRUE))
head(dt_plot)
# nan.idx<- which(rowSums(is.na(dt_plot))==TRUE)


unique(dt_plot$alg)
color_lib <- c("PCA"='#1f77b4',"LOL"='#ff7f0e', "PLS"='#2ca02c',
               "rrLDA"= '#9467bd', 'pca+pls'='#d62728', 'Naive_pca+pls'= '#8c564b')


line_lib <- c("PCA"="solid", "LOL"="solid", "PLS"="solid",
              "rrLDA"="solid", 'pca+pls'='solid', 'Naive_pca+pls'='solid')

dt.star <- dt_plot %>% filter(r == r.star)
dt_plot %>%
  # slice(nan.idx) %>%
  ggplot(aes(x=r, y=lhat, color=alg)) +
  geom_line(aes(linetype=alg), linewidth=0.7) +
  scale_linetype_manual(values=line_lib, name="Algorithm") +
  scale_color_manual(values=color_lib, name="Algorithm") +
  geom_point(data=dt.star, size=3, shape=16) +
  # facet_grid(.~sim) +
  facet_wrap(.~sim,scales = 'free',nrow = 2) +
  xlab("r") +
  ylab("error") +
  guides(color=guide_legend(ncol=1), linetype=guide_legend(ncol=1))
# labs(
#   title = paste('Classifier:',classifier.name,', n=',n,', d=',dim(X[2]),', rep=')  ) +
# theme(text=element_text(size=18))
# strip.text.x=element_blank(),
# strip.background=element_blank())


