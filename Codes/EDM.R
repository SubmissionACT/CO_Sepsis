##################################### preparation work #####################################
### packages
packages <- c('tidyverse','knitr','lubridate','dplyr','purrr','rEDM','doParallel',
              'foreach')
lapply(packages, require, character.only=T) 

### data and functions
load("sepsis0118.rda")
df <-  dfA %>% select(date, sep, temp,rh,co,fsp,o3h8max,no2)

### disease and pollutant variables
dis.name <- c("sep")
plt.name <- c("co","o3h8max","no2","fsp","temp","rh")

### parameters for multi_core running
cores_all <- detectCores()
cores <- ifelse(cores_all<9,4,18)
core_type <- 'PSOCK'

##################################### data preprocessing #####################################

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

df <- df %>% as.data.frame() %>% mutate_at(2:ncol(.),nomz) 

write.csv(df, file = "sepsis_TTTF.csv")

########################################## MFI ###########################################

### determine optimal E for S-map
fn_E_smapc <- function(data,y){
  df=data %>% select(date,y)
  M <- nrow(df)
  lib <- c(1,  M)
  pred <- c(1, M)
  E=EmbedDimension(dataFrame=df,lib=lib,pred=pred,
                   columns=y, target=y,maxE=10,showPlot=F)
  temp=data.frame(dis=y,E)
  temp
}

plist <- list(data=list(df),y=dis.name) %>% cross_df()
E_smapc_out <- plist %>% pmap_df(fn_E_smapc)

E_smapc=E_smapc_out %>%
  filter(E != 1) %>%
  group_by(dis) %>%
  mutate(row_number=1:n()) %>%
  filter(rho==max(rho))  %>%
  as.data.frame()

E_smapc

### determine optimal theta
fn_theta_justY <- function(data,y,theta){
  E <- E_smapc[E_smapc[,'dis']==y,'E']
  df <- data %>% select(date,y)
  m <- nrow(df)
  lib <- c(1,m)
  pred <- c(1,m)
  rho_theta <- PredictNonlinear(dataFrame = df,embedded = FALSE,
                                columns = y,target = y,Tp=1,
                                theta=theta,lib=lib,pred=pred,
                                showPlot = FALSE,E = E)
  best_theta_df <- rho_theta %>% mutate(dis=y,theta=theta)
}

plist <- list(data=list(df),y=dis.name,
              theta=c(0, 0.01, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3, 4,
                      5, 6, 7, 8, 9)) %>% cross_df()

cl <- makeCluster(cores[1],type = core_type)
registerDoParallel(cl)
theta_out <- foreach(i = 1:nrow(plist),
                     .packages = c("rEDM","tidyverse"),
                     .combine=rbind,
                     .inorder=FALSE)  %dopar% {
                       theta_out=plist[i,] %>% pmap_df(fn_theta_justY)
                     }
stopCluster(cl)

best_theta=theta_out %>%
  group_by(dis) %>%
  mutate(row_number=1:n()) %>%
  filter(rho==max(rho))  %>%
  select(-row_number) %>%
  as.data.frame()
best_theta

### forecast improvement
num_surr <- 100
num_sample <- 100

fn_mfi_1_smap_surr <- function(data,plt_1,dis,tp_value){
  E <- E_smapc[E_smapc$dis==dis,"E"]
  df <- data %>% mutate(plt=lag(.data[[plt_1]],-tp_value)) %>%
    filter(!(is.na(plt))) %>% select(date,dis,plt) %>%
    mutate(date=as.Date(date,'%Y-%m-%d'))
  set.seed(2019)
  surr_data <- SurrogateData(unlist(df[,'plt']), method = "ebisuzaki",
                             num_surr = num_surr, alpha = 1)
  all_data <- as.data.frame(cbind(df,as.data.frame(surr_data)))
  names(all_data) <- c("date", dis, 'T1',paste0("T", 2:(num_surr+1)))
  M <- nrow(all_data)
  out <- foreach(i = 1:(num_surr+1),
                 .packages = c("rEDM","tidyverse"),
                 .export=c('num_sample','num_surr','best_theta'),
                 .combine=rbind,
                 .inorder=FALSE)  %dopar% {
                   targetCol <- paste("T", i, sep = "")
                   embed_1 <- Embed(dataFrame = all_data, E = E, tau = -1, columns = dis )
                   dataFrame_1 <- cbind(all_data[E:M, 'date'],all_data[E:M, dis],embed_1[E:M, 1:(E)],
                                        all_data[E:M, targetCol]  ) %>% as.data.frame()
                   names(dataFrame_1) <- c('date',dis,letters[1:(E)],plt_1)
                   columns_1 <- paste(paste(letters[1:(E)],collapse =' '), plt_1 , sep=' ')
                   theta <- best_theta[best_theta[, "dis"] == dis, "theta"]
                   smap_1 <- SMap(dataFrame = dataFrame_1,embedded = TRUE,
                                  columns = columns_1,target = dis,
                                  lib = c(1, nrow(dataFrame_1)),
                                  pred = c(1, nrow(dataFrame_1)),
                                  theta = theta,Tp = 1)
                   out_1=compute_stats(smap_1$predictions$Observations,
                                       smap_1$predictions$Predictions) %>% mutate(mfi_012=1)
                   out=rbind(out_1) %>% mutate(i=i,plt_1=plt_1,dis=dis,tp_value=tp_value)
                   out
                 }
  out
}

plist <- list(data=list(df), dis=dis.name,plt_1=plt.name,
              tp_value=-3:0) %>% cross_df()

cl <- makeCluster(cores[1],type = core_type)
registerDoParallel(cl)
mfi_surr_out <- plist %>% pmap_df(fn_mfi_1_smap_surr)
stopCluster(cl)

p_function=function(x,pH,pL){
  dfg2_Q1=dfg2 %>% filter(grp=='surr') %>%
    group_by(dis,plt,tp_value,E) %>%
    summarise(ph=quantile(.data[[x]], pH,na.rm=T),
              pl=quantile(.data[[x]], pL,na.rm=T))

  dfg2_Q0=dfg2 %>% filter(grp=='raw') %>% select(-grp) %>%
    group_by(dis,plt,tp_value,E) %>%
    select(p0=x)

  dfg2_Q=left_join(dfg2_Q1,dfg2_Q0,by=c('E','tp_value',
                                        'dis','plt'))

  ggplot(dfg2_Q) +
    geom_hline(yintercept = 0,linetype='dashed') +
    geom_line(aes(tp_value,p0),size=0.5)+
    geom_point(aes(tp_value,p0),size=1)+
    geom_linerange(aes(x=tp_value,ymin =pl, ymax = ph),
                   size=1,color='green',alpha=0.5)+
    theme_bw() +
    ylab('MFI')+
    xlab("Lag")+
    theme(axis.text.x=element_text(size = 12, color= "black"),
          axis.text.y=element_text(size = 12, color= "black"),
          axis.title.x=element_text(size=15),
          axis.title.y=element_text(size=15),
          axis.ticks.length=unit(0.2,'cm'),
          legend.position="none",
          panel.background = element_blank(),
          plot.title = element_text(size = 18, hjust=0.5)) +
    facet_grid(dis~plt,scales='free')
}

dfg2 <- mfi_surr_out %>% filter(tp_value>=-3) %>% mutate(E=10) %>% mutate(i=i-1,grp=ifelse(i==0,'raw','surr'),
                                                                          rho=rho-best_theta$rho, plt=plt_1)

p_function(x='rho',pH=0.95,pL=0.05)
###################################### Effect size #######################################
best_theta = best_theta[, "theta"]

fn_smapc=function(data,plt,dis,tp_value){
  
  E = E_smapc[E_smapc[, "dis"] == dis, "E"]
  
  df=data %>%
    mutate(plt_tp=lag(.data[[plt]],-tp_value)) %>%
    filter(!(is.na(plt_tp))) %>%
    select(date,all_of(dis),plt_tp)
  
  M <- nrow(df)
  
  embed_1=Embed(dataFrame = df, E = E, tau = -1, columns = dis )
  
  dataFrame = cbind(df[E:M, 'date'],df[E:M, dis],
                    embed_1[E:M, 1:(E-1)],df[E:M, 'plt_tp']  ) %>%
    as.data.frame()
  
  names(dataFrame)=c('date',dis,letters[1:(E-1)],plt)
  
  columns = paste(paste(letters[1:(E-1)],collapse =' '), plt ,
                  sep=' ')
  
  m <- nrow(dataFrame)
  
  
  smap = SMap(dataFrame = dataFrame,
              embedded = TRUE,
              columns = columns,
              target = dis,
              lib = c(1, nrow(dataFrame)),
              pred=c(1, nrow(dataFrame)),
              theta=best_theta,
              Tp = 1,
              E = E)
  smapc_df=smap$coefficients[c(1,2+E)]
  
  names(smapc_df)=c('date','effect')
  smapc_df=smapc_df %>%
    mutate(date=lubridate::as_date(date,
                                   origin = lubridate::origin)) %>%
    mutate(dis=dis,plt=plt,E=E,tp_value=tp_value)
}

plist=list(data=list(df),
           dis=dis.name,
           plt=plt.name,
           tp_value=-3:0) %>% cross_df()

cl <- makeCluster(cores[1],type = core_type)
registerDoParallel(cl)

C_out=foreach(i = 1:nrow(plist),
              .packages = c("rEDM","tidyverse","lubridate"),
              .combine=rbind,
              # .export='num_sample',
              .inorder=FALSE)  %dopar% {
                C_out=plist[i,] %>% pmap_df(fn_smapc)
              }
stopCluster(cl)

tempA <- C_out %>%  filter(tp_value>=-3) %>% group_by(tp_value,E,dis,plt) %>%
  summarise(effect=mean(effect,na.rm=T))

ggplot(tempA,aes(tp_value,effect))+geom_line()+geom_hline(yintercept = 0) +
  ylab("Effect size")+
  xlab("lag")+
  facet_grid(~plt,scales='free')+
  theme_bw() +
  theme(axis.text.x=element_text(size = 12, color= "black"),
        axis.text.y=element_text(size = 12, color= "black"),
        axis.title.x=element_text(size=15),
        axis.title.y=element_text(size=15),
        axis.ticks.length=unit(0.2,'cm'),
        legend.position="none",
        panel.background = element_blank(),
        plot.title = element_text(size = 18, hjust=0.5))





