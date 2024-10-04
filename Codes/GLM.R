##################################### preparation work #####################################
### packages
packages <- c('tidyverse','dlnm','lubridate','tsModel','mgcv','scales','splines','ggsci')
lapply(packages, require, character.only=T) 

### data
load("sepsis0118.rda")
dfA <- df

dfA$sep <- log(dfA$sep)

### data standardization --- dtrend, dseason, and normalizaiton
nomz <- function(x, normalization=T, dseason=T, season_sd=F, sea=365, dtrend=T, dTtype="linear"){
  x <- as.numeric(x)
  xt <- x
  # Detrend
  if(dtrend==T & dTtype=="first"){xt <- diff(xt)} else if (dtrend==T & dTtype=="linear"){
    lm.t <- lm(xt~c(1:length(xt)))
    xt <- xt-(lm.t$coefficients[1]+lm.t$coefficients[2]*c(1:length(xt)))}
  # Deseason
  if(dseason==T){
    xs <- as.numeric(apply(matrix(xt[1:(sea*length(xt)%/%sea)],ncol=sea,byrow=T),2,mean,na.rm=T))
    xsd <- as.numeric(apply(matrix(xt[1:(sea*length(xt)%/%sea)],ncol=sea,byrow=T),2,sd,na.rm=T))
    xt <- xt-c(rep(xs,1+length(xt)%/%sea))[1:length(xt)]
    if(season_sd==T){xt <- xt/(c(rep(xsd,1+length(xt)%/%sea))[1:length(xt)])}}
  # Normalization (zero mean & unity variance)
  if(normalization==T){xt <- (xt-mean(xt,na.rm=T))/sd(xt,na.rm=T)}
  return(xt)
}

names(dfA)

dfB=dfA %>% select(date,dow,holiday,sep,temp,rh,co,fsp,o3=o3h8max) 

data <- dfB %>% mutate_at(4:ncol(.),nomz) %>% na.omit()

sumA <- function(fitX) {
  
  summ=summary(fitX)
  
  outA <- data.frame(
    beta = summ$coefficients[2,1],
    se   = summ$coefficients[2,2],
    t    = summ$coefficients[2,3],
    p    = summ$coefficients[2,4]
  )
  
}

### co 
fit_co_0=glm(formula = sep ~ lag(co,0) + lag(co,1) + lag(co,2) + lag(co,3) + temp+ factor(dow) + factor(holiday), data=data, family=gaussian(), na.action=na.omit)
out_co_0=sumA(fit_co_0) %>% mutate(plt='co',lag=0);  out_co_0

fit_co_1=glm(formula = sep ~ lag(co,1) + lag(co,0) + lag(co,2) + lag(co,3)+ factor(dow) + factor(holiday), data=data, family=gaussian(), na.action=na.omit)
out_co_1=sumA(fit_co_1) %>% mutate(plt='co',lag=1);  out_co_1

fit_co_2=glm(formula = sep ~ lag(co,2) + lag(co,0) + lag(co,1) + lag(co,3) + lag(temp,3)+ factor(dow) + factor(holiday), data=data, family=gaussian(), na.action=na.omit)
out_co_2=sumA(fit_co_2) %>% mutate(plt='co',lag=2);  out_co_2

fit_co_3=glm(formula = sep ~ lag(co,3) + lag(co,0) + lag(co,1) + lag(co,2) + lag(temp,3)+ factor(dow) + factor(holiday), data=data, family=gaussian(), na.action=na.omit)
out_co_3=sumA(fit_co_3) %>% mutate(plt='co',lag=3);  out_co_3

### o3
fit_o3_0=glm(formula = sep ~ lag(o3,0) + lag(o3,1) + lag(o3,2) + lag(o3,3)+ factor(dow) + factor(holiday), data=data, family=gaussian(), na.action=na.omit)
out_o3_0=sumA(fit_o3_0) %>% mutate(plt='o3',lag=0);  out_o3_0

fit_o3_1=glm(formula = sep ~ lag(o3,1) + lag(o3,0) + lag(o3,2) + lag(o3,3)+ factor(dow) + factor(holiday), data=data, family=gaussian(), na.action=na.omit)
out_o3_1=sumA(fit_o3_1) %>% mutate(plt='o3',lag=1);  out_o3_1

fit_o3_2=glm(formula = sep ~ lag(o3,2) + lag(o3,1) + lag(o3,0) + lag(o3,3)+ factor(dow) + factor(holiday), data=data, family=gaussian(), na.action=na.omit)
out_o3_2=sumA(fit_o3_2) %>% mutate(plt='o3',lag=2);  out_o3_2

fit_o3_3=glm(formula = sep ~ lag(o3,3) + lag(o3,1) + lag(o3,2) + lag(o3,0)+ factor(dow) + factor(holiday), data=data, family=gaussian(), na.action=na.omit)
out_o3_3=sumA(fit_o3_3) %>% mutate(plt='o3',lag=3);  out_o3_3

### fsp
fit_fsp_0=glm(formula = sep ~ lag(fsp,0) + lag(fsp,1) + lag(fsp,2) + lag(fsp,3)+ lag(temp,0)+ factor(dow) + factor(holiday), data=data, family=gaussian(), na.action=na.omit)
out_fsp_0=sumA(fit_fsp_0) %>% mutate(plt='fsp',lag=0);  out_fsp_0

fit_fsp_1=glm(formula = sep ~ lag(fsp,1) + lag(fsp,0) + lag(fsp,2) + lag(fsp,3) + lag(temp,0)+ lag(co,1)+ factor(dow) + factor(holiday), data=data, family=gaussian(), na.action=na.omit)
out_fsp_1=sumA(fit_fsp_1) %>% mutate(plt='fsp',lag=1);  out_fsp_1

fit_fsp_2=glm(formula = sep ~ lag(fsp,2) + lag(fsp,1) + lag(fsp,0) + lag(fsp,3)+ lag(temp,3)+ factor(dow) + factor(holiday), data=data, family=gaussian(), na.action=na.omit)
out_fsp_2=sumA(fit_fsp_2) %>% mutate(plt='fsp',lag=2);  out_fsp_2

fit_fsp_3=glm(formula = sep ~ lag(fsp,3) + lag(fsp,1) + lag(fsp,2) + lag(fsp,0)+ lag(temp,3)+ factor(dow) + factor(holiday), data=data, family=gaussian(), na.action=na.omit)
out_fsp_3=sumA(fit_fsp_3) %>% mutate(plt='fsp',lag=3);  out_fsp_3

### temp
fit_temp_0=glm(formula = sep ~ lag(temp,0) + lag(temp,1) + lag(temp,2) + lag(temp,3)+ factor(dow) + factor(holiday), data=data, family=gaussian(), na.action=na.omit)
out_temp_0=sumA(fit_temp_0) %>% mutate(plt='temp',lag=0);  out_temp_0

fit_temp_1=glm(formula = sep ~ lag(temp,1) + lag(temp,0) + lag(temp,2) + lag(temp,3)+ factor(dow) + factor(holiday), data=data, family=gaussian(), na.action=na.omit)
out_temp_1=sumA(fit_temp_1) %>% mutate(plt='temp',lag=1);  out_temp_1

fit_temp_2=glm(formula = sep ~ lag(temp,2) + lag(temp,1) + lag(temp,0) + lag(temp,3)+ factor(dow) + factor(holiday), data=data, family=gaussian(), na.action=na.omit)
out_temp_2=sumA(fit_temp_2) %>% mutate(plt='temp',lag=2);  out_temp_2

fit_temp_3=glm(formula = sep ~ lag(temp,3) + lag(temp,1) + lag(temp,2) + lag(temp,0)+ factor(dow) + factor(holiday), data=data, family=gaussian(), na.action=na.omit)
out_temp_3=sumA(fit_temp_3) %>% mutate(plt='temp',lag=3);  out_temp_3

### rh
fit_rh_0=glm(formula = sep ~ lag(rh,0) + lag(rh,1) + lag(rh,2) + lag(rh,3) + lag(temp,3)+ factor(dow) + factor(holiday), data=data, family=gaussian(), na.action=na.omit)
out_rh_0=sumA(fit_rh_0) %>% mutate(plt='rh',lag=0);  out_rh_0

fit_rh_1=glm(formula = sep ~ lag(rh,1) + lag(rh,0) + lag(rh,2) + lag(rh,3)+ factor(dow) + factor(holiday), data=data, family=gaussian(), na.action=na.omit)
out_rh_1=sumA(fit_rh_1) %>% mutate(plt='rh',lag=1);  out_rh_1

fit_rh_2=glm(formula = sep ~ lag(rh,2) + lag(rh,1) + lag(rh,0) + lag(rh,3)+ factor(dow) + factor(holiday), data=data, family=gaussian(), na.action=na.omit)
out_rh_2=sumA(fit_rh_2) %>% mutate(plt='rh',lag=2);  out_rh_2

fit_rh_3=glm(formula = sep ~ lag(rh,3) + lag(rh,1) + lag(rh,2) + lag(rh,0)+ factor(dow) + factor(holiday), data=data, family=gaussian(), na.action=na.omit)
out_rh_3=sumA(fit_rh_3) %>% mutate(plt='rh',lag=3);  out_rh_3

###################################################### plotting ##########################################
res <- rbind(out_co_0,
             out_co_1,
             out_co_2,
             out_co_3,
             out_temp_0,
             out_temp_1,
             out_temp_2,
             out_temp_3,
             out_rh_0,
             out_rh_1,
             out_rh_2,
             out_rh_3,
             out_o3_0,
             out_o3_1,
             out_o3_2,
             out_o3_3,
             out_fsp_0,
             out_fsp_1,
             out_fsp_2,
             out_fsp_3)

GLMresCI  =  res %>%
  mutate(
         betalow=beta-1.96*se, 
         betahigh=beta+1.96*se) %>% 
  select(lag,beta,betalow,betahigh,p,plt) %>% mutate(sig=ifelse(p<0.05,1,0))

iqr <- data %>%
  select(-c(date, dow, holiday, sep)) %>%
  gather(key = "plt", value = "value") %>%
  group_by(plt) %>%
  summarise(iqr1 = IQR(value, na.rm = TRUE))

GLMresCI <- GLMresCI %>% left_join(iqr, by="plt") %>% mutate(beta=beta*iqr1,
                                                             betalow=betalow*iqr1,
                                                             betahigh=betahigh*iqr1)


GLMresCI$sig <- factor(GLMresCI$sig)
GLMresCI$plt <- factor(GLMresCI$plt, levels = c("co","o3","fsp","temp","rh"), labels = c(expression(paste("CO",sep = " ")),
                                                                                         expression(paste("O"[3],sep = " ")),
                                                                                         expression(paste("PM"[2.5],sep = " ")),
                                                                                         expression(paste("Temp",sep = " ")),
                                                                                         expression(paste("Humid",sep = " "))))

GLMresCI$lag <- as.numeric(as.character(GLMresCI$lag)) 

P2_2 <- ggplot(GLMresCI, aes(x = lag, y = beta)) +
  geom_hline(yintercept = 0, linetype = 'dashed', col="gray") +
  geom_ribbon(aes(ymin = betalow, ymax = betahigh), fill = "#708090", alpha = 0.3, color = NA) +
  geom_line(aes(group = plt), color = "black", alpha=0.7) +
  geom_errorbar(data = GLMresCI[GLMresCI$sig == 1, ], aes(ymin = betalow, ymax = betahigh), width = 0.2, color = "red") +
  geom_point(aes(shape = factor(sig), color = factor(sig == 1)), size = 2, alpha=0.7) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  scale_shape_manual(values = c(1, 8)) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.001)) +
  facet_wrap(~ plt, scales = 'fixed', labeller = "label_parsed", ncol = 5) +
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


