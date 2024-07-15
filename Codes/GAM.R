
##################################### preparation work #####################################
### packages
packages <- c('tidyverse','dlnm','lubridate','tsModel','mgcv','scales','splines','ggsci')
lapply(packages, require, character.only=T) 

### data
load("sepsis0118.rda")
dfA$t <- 1:length(dfA$date)
data=dfA
print(data)

################################## Analysis for air pollutants ################################
HTadjust <- function(dis, plt, df_t){
  
  cb.plt = crossbasis(data[[plt]], lag=3, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.coplt2 = crossbasis(data[["rh"]], lag=3, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.coplt1 = crossbasis(data[["temp"]], lag=3, argvar=list(fun="thr",thr=c(18,28)), arglag=list(fun="integer"))
  
  plist <-  "cb.plt+cb.coplt1 +cb.coplt2"
  
  eval(parse(text=paste("model <- gam(", dis," ~ ", plist, "+
        s(t,k=",df_t,"*18)+
        as.factor(dow)+as.factor(holiday),
        family=quasipoisson(link = 'log'), data=data)", sep='')))
  
  iqr<-IQR(data[[plt]],na.rm=TRUE)
  eval(parse(text=paste("pred <- crosspred(cb.plt, model,
                        at=iqr, cumul=TRUE)",sep='')))
  matRRfit <- data.frame(fit = "matfit", dis = dis, plt = plt, iqr=iqr,
                         coplt = "None",df_t=df_t,
                         pred$matfit, stringsAsFactors = FALSE)
  matRRhigh <- data.frame(fit = "mathigh", dis = dis, plt = plt, iqr=iqr,
                          coplt = "None",df_t=df_t,
                          pred$matfit+1.96*pred$matse, stringsAsFactors = FALSE)
  matRRlow <- data.frame(fit = "matlow", dis = dis, plt = plt,iqr=iqr,
                         coplt = "None",df_t=df_t,
                         pred$matfit-1.96*pred$matse, stringsAsFactors = FALSE)
  out <- bind_rows(matRRfit, matRRhigh, matRRlow)
  out
}

plist <- list(dis = c('sep'),
              plt = c( "co","fsp","o3h8max","no2"),
              df_t=c(7))

bb <- plist %>% cross_df()  %>%  pmap_df(HTadjust)

bbA <- bb %>%
  gather(lagA, value, starts_with("lag")) %>%
  spread(fit, value) %>%
  mutate(lag = parse_number(lagA)) %>%
  mutate(lag=factor(lag))

dat1 <- bbA

################################## Analysis for temperature ################################
T_COadjust <- function(dis, plt, df_t){
  
  cb.coplt1 = crossbasis(data[[plt]], lag=3, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.coplt2 = crossbasis(data[["rh"]], lag=3, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.plt = crossbasis(data[["temp"]], lag=3, argvar=list(fun="thr",thr=c(18,28)), arglag=list(fun="integer"))
  
  plist <-  "cb.plt+cb.coplt1+cb.coplt2"
  
  eval(parse(text=paste("model <- gam(", dis," ~ ", plist, "+
        s(t,k=",df_t,"*18)+
        as.factor(dow)+as.factor(holiday),
        family=quasipoisson(link = 'log'), data=data)", sep='')))
  
  iqr<-IQR(data[["temp"]],na.rm=TRUE)
  eval(parse(text=paste("pred <- crosspred(cb.plt, model,
                        at=36, cumul=TRUE)",sep='')))
  matRRfit <- data.frame(fit = "matfit", dis = dis, plt = "temp", iqr=iqr,
                         coplt = "None",df_t=df_t,
                         pred$matfit, stringsAsFactors = FALSE)
  matRRhigh <- data.frame(fit = "mathigh", dis = dis, plt = "temp", iqr=iqr,
                          coplt = "None",df_t=df_t,
                          pred$matfit+1.96*pred$matse, stringsAsFactors = FALSE)
  matRRlow <- data.frame(fit = "matlow", dis = dis, plt = "temp",iqr=iqr,
                         coplt = "None",df_t=df_t,
                         pred$matfit-1.96*pred$matse, stringsAsFactors = FALSE)
  out <- bind_rows(matRRfit, matRRhigh, matRRlow)
  out
}

plist <- list(dis = c('sep'),
              plt = c( "co"),
              df_t=c(7))

bb <- plist %>% cross_df()  %>%  pmap_df(T_COadjust)

bbA <- bb %>%
  gather(lagA, value, starts_with("lag")) %>%
  spread(fit, value) %>%
  mutate(lag = parse_number(lagA)) %>%
  mutate(lag=factor(lag))

dat2 <- bbA

dat2 <- dat2 %>% mutate(matfit = matfit/10, mathigh=mathigh/10, matlow=matlow/10)

################################## Analysis for relative humidity ################################
RH_COadjust <- function(dis, plt,df_t){
  
  cb.coplt1 = crossbasis(data[[plt]], lag=3, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.plt = crossbasis(data[["rh"]], lag=3, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.coplt2 = crossbasis(data[["temp"]], lag=3, argvar=list(fun="thr",thr=c(18,28)), arglag=list(fun="integer"))
  
  plist <-  "cb.plt+cb.coplt1+cb.coplt2"
  
  eval(parse(text=paste("model <- gam(", dis," ~ ", plist, "+
        s(t,k=",df_t,"*18)+
        as.factor(dow)+as.factor(holiday),
        family=quasipoisson(link = 'log'), data=data)", sep='')))
  
  iqr<-IQR(data[["rh"]],na.rm=TRUE)
  eval(parse(text=paste("pred <- crosspred(cb.plt, model,
                        at=iqr, cumul=TRUE)",sep='')))
  matRRfit <- data.frame(fit = "matfit", dis = dis, plt ="rh", iqr=iqr,
                         coplt = "None",df_t=df_t,
                         pred$matfit, stringsAsFactors = FALSE)
  matRRhigh <- data.frame(fit = "mathigh", dis = dis, plt = "rh", iqr=iqr,
                          coplt = "None",df_t=df_t,
                          pred$matfit+1.96*pred$matse, stringsAsFactors = FALSE)
  matRRlow <- data.frame(fit = "matlow", dis = dis, plt = "rh",iqr=iqr,
                         coplt = "None",df_t=df_t,
                         pred$matfit-1.96*pred$matse, stringsAsFactors = FALSE)
  out <- bind_rows(matRRfit, matRRhigh, matRRlow)
  out
}

plist <- list(dis = c('sep'),
              plt = c( "co"),
              df_t=c(7))

bb <- plist %>% cross_df()  %>%  pmap_df(RH_COadjust)

bbA <- bb %>%
  gather(lagA, value, starts_with("lag")) %>%
  spread(fit, value) %>%
  mutate(lag = parse_number(lagA)) %>%
  mutate(lag=factor(lag))

dat3 <- bbA

################################## plotting ################################
dat <- rbind(dat1,dat2,dat3)

dat <- dat %>% mutate(Significance=as.factor(ifelse(mathigh<=0,1,ifelse(matlow>=0,1,0))))

dat$plt <- factor(dat$plt, levels=c("co","o3h8max","fsp","no2","temp","rh"),
                  labels = c(
                    expression(paste("CO",sep = " ")),
                    expression(paste("O"[3],sep = " ")),
                    expression(paste("PM"[2.5],sep = " ")),
                    expression(paste("NO"[2],sep = " ")),
                    expression(paste("Temp",sep = " ")),
                    expression(paste("Humid",sep = " "))
                  ))

dat$lag <- as.numeric(as.character(dat$lag))  

P2_2 <- ggplot(dat, aes(x = lag, y = matfit)) +
  geom_hline(yintercept = 0, linetype = 'dashed', col="gray") +
  geom_ribbon(aes(ymin = matlow, ymax = mathigh), fill = "#708090", alpha = 0.3, color = NA) +
  geom_line(aes(group = plt), color = "black", alpha=0.7) +
  geom_errorbar(data = dat[dat$Significance == 1, ], aes(ymin = matlow, ymax = mathigh), width = 0.2, color = "red") +
  geom_point(aes(shape = factor(Significance), color = factor(Significance == 1)), size = 2, alpha=0.7) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  scale_shape_manual(values = c(1, 8)) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.001)) +
  facet_wrap(~ plt, scales = 'fixed', labeller = "label_parsed", ncol = 6) +
  ylab(expression(paste("Sepsis risk (", beta, ")", sep = ""))) +
  scale_x_continuous(breaks = c(0, 1, 2, 3), labels = c("0", "1", "2", "3"), limits = c(-0.5, 3.5)) +
  xlab("") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8, color = "darkgray"),
    axis.text.y = element_text(size = 8, color = "darkgray"),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.ticks.length = unit(0.2, 'cm'),
    axis.ticks.x = element_line(color = "darkgray"),
    axis.ticks.y = element_line(color = "darkgray"),
    axis.line = element_line(colour = "darkgray"),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.box = "vertical",
    legend.margin = margin(-15, unit = "pt"),
    legend.box.spacing = margin(15.5),
    legend.background = element_rect(fill = "transparent"),
    strip.background = element_blank(),
    panel.border = element_blank(),
    legend.key = element_rect(fill = "transparent"),
    strip.text = element_text(size = 15),
    panel.grid = element_blank(),
    plot.title = element_text(size = 40, hjust = 0.5)
  )

x_lims <- c(-0.4, 3.4)
P2_2 <- P2_2 + scale_x_continuous(limits = x_lims, breaks = c(0, 1, 2, 3))
print(P2_2)

