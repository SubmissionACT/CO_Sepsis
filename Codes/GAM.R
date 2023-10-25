################################## Load packages and data ###################################
library(tidyverse); library(dlnm);
library(lubridate); library(tsModel);
library(mgcv); library(scales); library(splines)
library(ggsci)

load("sepsis0118.rda")

dfA$t <- 1:length(dfA$date)
data=dfA

################################## Analysis for air pollutants ################################
FTjust <- function(dis, plt, df_t){
  
  cb.plt = crossbasis(data[[plt]], lag=3, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.coplt2 = crossbasis(data[["rh"]], lag=3, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.coplt1 = crossbasis(data[["temp"]], lag=3, argvar=list(fun="thr",thr=c(18,28)),
                         arglag=list(fun="integer"))
  
  plist <-  "cb.plt+cb.coplt1 +cb.coplt2"
  
  eval(parse(text=paste("model <- gam(", dis," ~ ", plist, "+
        s(t,k=",df_t,"*18,fx=T,bs='cr')+

        as.factor(dow)+as.factor(holiday),
 
        family=quasipoisson(link = 'log'), scale=-1,
        control=gam.control(epsilon=0.0000001,maxit=200),
        na.action= na.omit, data=data)", sep='')))
  
  
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

bb <- plist %>% cross_df()  %>%  pmap_df(FTjust)

bbA <- bb %>%
  gather(lagA, value, starts_with("lag")) %>%
  spread(fit, value) %>%
  mutate(lag = parse_number(lagA)) %>%
  mutate(lag=factor(lag))

dat1 <- bbA
dat1

################################## Analysis for temperature ################################
FTjust <- function(dis, plt, df_t){
  
  cb.coplt1 = crossbasis(data[[plt]], lag=3, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.coplt2 = crossbasis(data[["rh"]], lag=3, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.plt = crossbasis(data[["temp"]], lag=3, argvar=list(fun="thr",thr=c(18,28)),
                      arglag=list(fun="integer"))
  plist <-  "cb.plt+cb.coplt1+cb.coplt2"
  
  eval(parse(text=paste("model <- gam(", dis," ~ ", plist, "+
        s(t,k=",df_t,"*18,fx=T,bs='cr')+

        as.factor(dow)+as.factor(holiday),
 
        family=quasipoisson(link = 'log'), scale=-1,
        control=gam.control(epsilon=0.0000001,maxit=200),
        na.action= na.omit, data=data)", sep='')))
  
  
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

bb <- plist %>% cross_df()  %>%  pmap_df(FTjust)

bbA <- bb %>%
  gather(lagA, value, starts_with("lag")) %>%
  spread(fit, value) %>%
  mutate(lag = parse_number(lagA)) %>%
  mutate(lag=factor(lag))

dat2 <- bbA
dat2

################################## Analysis for relative humidity ################################
FTjust <- function(dis, plt,df_t){
  
  cb.coplt1 = crossbasis(data[[plt]], lag=3, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.plt = crossbasis(data[["rh"]], lag=3, argvar=list(fun="lin"), arglag=list(fun="integer"))
  cb.coplt2 = crossbasis(data[["temp"]], lag=3, argvar=list(fun="thr",thr=c(18,28)),
                         arglag=list(fun="integer"))
  
  plist <-  "cb.plt+cb.coplt1+cb.coplt2"
  
  eval(parse(text=paste("model <- gam(", dis," ~ ", plist, "+
        s(t,k=",df_t,"*18,fx=T,bs='cr')+

        as.factor(dow)+as.factor(holiday),
 
        family=quasipoisson(link = 'log'), scale=-1,
        control=gam.control(epsilon=0.0000001,maxit=200),
        na.action= na.omit, data=data)", sep='')))
  
  
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

bb <- plist %>% cross_df()  %>%  pmap_df(FTjust)

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
                    expression(paste("Temp.",sep = " ")),
                    expression(paste("Humid.",sep = " "))
                  ))

P2_2 <- ggplot(dat, aes(lag,matfit,ymin = matlow, ymax = mathigh)) +
  geom_hline(yintercept = 0,linetype='dashed') +
  geom_errorbar(aes(group = plt, col = plt, width=0), 
                size=1, position = position_dodge((width=0.5))) +
  scale_color_jco()+
  geom_point(aes(group = plt, col = plt,shape=Significance),size = 2.3, 
             position = position_dodge(width=0.5)) +
  scale_shape_manual(values = c(16,8)) +
  scale_y_continuous(labels=label_number(accuracy =0.01))+
  facet_wrap(~plt,scales='free',labeller = "label_parsed",ncol =6 ) +
  ylab(expression(paste(beta," (Risk estimates for env. var. → Sepsis)", sep = "")))+
  xlab("Lag")+
  theme(axis.text.x=element_text(size = 12, color= "black"),
        axis.text.y=element_text(size = 12, color= "black"),
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.ticks.length=unit(0.2,'cm'),
        axis.line = element_line(colour = "black"),
        legend.position="none",
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.box="vertical",
        legend.margin=margin(-15,unit="pt"),
        legend.box.spacing = margin(15.5),
        legend.background = element_rect(fill="transparent"),
        strip.background = element_rect(
          color = "white", fill = "white"),
        panel.border = element_blank(),
        legend.key = element_rect(fill = "transparent"),
        strip.text = element_text(size = 14),
        panel.grid = element_blank(),
        plot.title = element_text(size = 18, hjust=0.5))
P2_2

