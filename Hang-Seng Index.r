



###### Leer DATOS - HANG SENG INDEX (HONG KONG) #####

data <- read.csv("^HSI.csv",header = TRUE,sep=";")

price <- data$HSI.Close
ret <- diff(price) #First Difference

par(mfrow=c(1,2))
plot.ts(price)
plot.ts(ret) #Series not stationary, transformation needed.

##### Dickey - Fuller test para ver si es estacionaria

adf.test(price)#accept null hypothesis, confirmed non stationarity.
adf.test(ret)#reject null hypothesis, series is now stationary


############################### MODELADO ###############################


#Correlograma
par(mfrow=c(1,2))
Acf(ret, lag.max = 50, main="ACF")
Pacf(ret, lag.max = 50, main="PACF")

fit.auto <- auto.arima(ret, ic = c("aic"), stationary = TRUE) #Choose model based on AIC
summary(fit.auto) #ARIMA(1,0,0) (e.g. AR(1) process)

fit.manual <- arima(ret,order=c(1,0,0))

#Correlograma residuales
par(mfrow=c(1,2))
acf(fit.manual$residuals,lag.max = 50,main="ACF") #Both plots behave like AR(1) process
pacf(fit.manual$residuals,lag.max = 50, main="PACF") 

#Test serial correlation
Box.test(fit.manual$residuals,  lag = 50, type = c("Ljung-Box")) #Accept null hypothesis, white noise among residuals

##### Revisar normalidad de los errores ####
checkresiduals(fit.manual)

resid.norm <- fit.manual$residuals/sqrt(fit.manual$sigma2)

par(mfrow = c(1,2))
hist(resid.norm, prob = TRUE, breaks = 100, xlab = "Normalized residuals", main="Histogram")
curve(dnorm(x, mean = 0, sd = 1), col="blue", lwd = 1, add = TRUE)
qqnorm(resid.norm)
qqline(resid.norm) #residuals are far from normal

#### Es Gaussiano? #####
jarque.bera.test(resid.norm) #null hyphotesis (NORMALITY) rejected
kurtosis(resid.norm) #fat tails: high chance of extreme events
skewness(resid.norm) #negative skewness


########################## MODELADO DE VOLATILIDAD #############################


resid.sqrd <- fit.manual$residuals^2

dev.off()
plot.ts(resid.sqrd, ylab = "Squared residuals") #Clusters de volatilidad encontrados

#Correlograma residuales al cuadrado
par(mfrow=c(1,2))
Acf(resid.sqrd, lag.max = 20, main="ACF")
Pacf(resid.sqrd, lag.max = 20, main="PACF")


##### ¿Hay efecto ARCH? #####
Box.test(resid.sqrd,  lag = 8, type = c("Ljung-Box")) #Serial correlation, reject null hypothesis
arch.test(fit.manual, output = T) #ARCH effects, reject null hypothesis as well

###################### Encontrar el mejor modelo de acuerdo a AIC ######################

p.max <- 5
q.max <- 5
aic.min <- Inf
best.p <- 0
best.q <- 0
tab <- matrix(rep(0,25),nrow=5)
for (i1 in 1:p.max) {
  for (i2 in 1:q.max) {
    ourSpec=ugarchspec(mean.model = list(armaOrder = c(1, 0), include.mean = FALSE), 
                       variance.model = list(garchOrder = c(i1, i2)))
    fit = ugarchfit(spec = ourSpec, data = fit.manual$residuals)
    tab[i1,i2]=infocriteria(fit)[1]
    inf.crit=infocriteria(fit)[1]
    aic.min=ifelse(inf.crit < aic.min, inf.crit, aic.min)
    best.p=ifelse(inf.crit == aic.min, i1, best.p)
    best.q=ifelse(inf.crit == aic.min, i2, best.q)
  }
}
rownames(tab) <- c('1','2','3','4','5')
colnames(tab) <- c('1','2','3','4','5')
tab <- as.table(tab) #AIC values are vaguely different
tab
c(best.p, best.q) #Not GARCH (1,1) but still as parsimonious


############################ ARIMA(1,0,0)- GARCH(2,1) ############################

ourSpec <- ugarchspec(mean.model = list(armaOrder = c(1, 0)), 
                      variance.model = list(garchOrder = c(best.p, best.q)))

fit_1 <- ugarchfit(spec = ourSpec, data = fit.manual$residuals)

fit_2 <- ugarchfit(spec = ugarchspec(mean.model = list(armaOrder = c(1, 0)), 
                                     variance.model = list(garchOrder = c(1, 1))),data = fit.manual$residuals)

fit_3 <- ugarchfit(spec = ugarchspec(mean.model = list(armaOrder = c(1, 0), include.mean = FALSE), 
                                     variance.model = list(model = "gjrGARCH", garchOrder = c(best.p, best.q)),
                                     distribution.model = "std"), data = fit.manual$residuals)

fit_4 <- ugarchfit(spec = ugarchspec(mean.model = list(armaOrder = c(1, 0), include.mean = FALSE), 
                                     variance.model = list(model = "tGARCH", garchOrder = c(best.p, best.q)),
                                     distribution.model = "std"), data = fit.manual$residuals)



par(mfrow=c(2,2))
plot(fit)
10
11
9
12
0


#HACER gjrGARCH tGARCH Y GARCH(1,1) PARA COMPARAR

