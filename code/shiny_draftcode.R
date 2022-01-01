library(dplyr)
library(tidyr)
dat <- data.frame(id = 1:1000)
set.seed(2022)
dat$x <- rbinom(1000, 1, 0.5)
table(dat$x)
set.seed(2021)
for(i in 1:length(dat$x)){
  dat$y[i] <- ifelse(dat$x[i] == 1, 
                     sample(c(1,2,3,4), 1000, T, c(0.1, 0.15, 0.25, 0.5)),
                     sample(c(1,2,3,4), 1000, T, c(0.3, 0.4, 0.2, 0.1)))
}
table(dat$y, dat$x)

fit1<-glm(y_lev1 ~ x, family = binomial, 
          data = dat %>% 
            mutate(y_lev1 = factor(ifelse(y == 1, 1, 0), levels=c(0,1))))
fit2<-glm(y_lev2 ~ x, family = binomial, 
          data = dat %>% 
            filter(y > 1) %>% 
            mutate(y_lev2 = factor(ifelse(y == 2, 1, 0), levels = c(0,1))))
fit3<-glm(y_lev3 ~ x, family = binomial, 
          data = dat %>% 
            filter(y > 2) %>% 
            mutate(y_lev3 = factor(ifelse(y == 3, 1, 0), levels=c(0,1))))

plot_x0 <- data.frame(id=1, x = 0)
p_y1_given_x0 <- predict(fit1, newdata=plot_x0, type = "response")
p_y2_given_x0ge1 <- predict(fit2, newdata=plot_x0, type = "response")
p_y3_given_x0ge2 <- predict(fit3, newdata=plot_x0, type = "response")

plot_x0$p1 <- p_y1_given_x0
plot_x0$p2 <- (1-p_y1_given_x0) * p_y2_given_x0ge1
plot_x0$p3 <- (1-p_y1_given_x0)*(1-p_y2_given_x0ge1) * p_y3_given_x0ge2
plot_x0$p4 <- 1-plot_x0$p1-plot_x0$p2-plot_x0$p3
plot_x0 %>% 
  gather("level", "prob", p1:p4) %>% 
  ggplot()+
  geom_bar(aes(x = level, y = prob),stat="identity")+theme_classic()+ylim(c(0,1))

plot_x1 <- data.frame(id=2, x = 1)
p_y1_given_x1 <- predict(fit1, newdata=plot_x1, type = "response")
p_y2_given_x1ge1 <- predict(fit2, newdata=plot_x1, type = "response")
p_y3_given_x1ge2 <- predict(fit3, newdata=plot_x1, type = "response")

plot_x1$p1 <- p_y1_given_x1
plot_x1$p2 <- (1-p_y1_given_x1) * p_y2_given_x1ge1
plot_x1$p3 <- (1-p_y1_given_x1) * (1-p_y2_given_x1ge1) * p_y3_given_x1ge2
plot_x1$p4 <-1-plot_x1$p1-plot_x1$p2-plot_x1$p3
plot_x1 %>% 
  gather("level", "prob", p1:p4) %>% 
  ggplot()+
  geom_bar(aes(x = level, y = prob),stat="identity")+theme_classic()+ylim(c(0,1))


m <- polr(y ~x, data = dat %>% mutate(y=  factor(y, levels= c(1,2,3,4))))
predict(m, newdata = plot_x0, type = "probs")
predict(m, newdata = plot_x1, type = "probs")
