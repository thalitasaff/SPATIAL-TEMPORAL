# data : bismillah uass

# =========================================================
## 1. PACKAGE 
# =========================================================
library(car)
library(foreign)
library(shapefiles)
library(nlme)
library(maptools)
library(grid)
library(gstat)
library(latticeExtra)
library(splines)
library(spdep)
library(classInt)
library(RColorBrewer)
library(CARBayes)
library(CARBayesST)
library(INLA)
library(DCluster)
library(spdep)
library(CAMAN)
library(sp)
library(ctv)
library(rgdal)
library(Matrix)

# =========================================================
## 2. MASUKKAN DATA
# =========================================================
data <- read.csv(file.choose(), header=T, sep=";", dec=".")
attach(data)
head(data)

summary(data[1:30,])
summary(data[31:60,])
summary(data[61:90,])

# =========================================================
## 3. SET DATA MAPPING 
# =========================================================
Indo <- readRDS("/Users/ASUS/Downloads/gadm36_IDN_3_sp.rds")
head(Indo,10)
Bandung <- Indo[Indo$NAME_2 == "Kota Bandung",]
head(Bandung,10)
Kec <- Bandung$NAME_3 
row.names(Bandung) <- Kec
COOR <- coordinates(Bandung)

plot(Bandung, axes=T, col="gray90")
text(COOR[,1],COOR[,2], row.names(Bandung), col="black", cex=0.5, pos=1)
points(COOR[,1],COOR[,2], pch=19, cex=0.5,col="blue")

# =========================================================
## 4. MENGHITUNG STANDARDIZED MORBIDITY RASIO (SMR)
# =========================================================
# Estimasi SMR
SMR14 <- (data$JK[1:30]/data$Ei[1:30]) ;SMR14
SMR15 <- (data$JK[31:60]/data$Ei[31:60]) ;SMR15
SMR16 <- (data$JK[61:90]/data$Ei[61:90]) ;SMR16
SMR <- c(SMR14,SMR15,SMR16)
data.frame(Kecamatan,SMR14,SMR15,SMR16)

# Menghitung Standard Error SMR
var.SMR14 <- SMR14/Ei[1:30]
sd.SMR14 <- sqrt(var.SMR14) ;sd.SMR14
var.SMR15 <- SMR15/Ei[31:60]
sd.SMR15 <- sqrt(var.SMR15) ;sd.SMR15
var.SMR16 <- SMR16/Ei[61:90]
sd.SMR16 <- sqrt(var.SMR16) ;sd.SMR16

# Membuat boxplot dan Plot Taksiran SMR
boxplot(SMR14,SMR15,SMR16, main = "Boxplot SMR",
        names = c("SMR14","SMR15","SMR16"))

# =========================================================
## 5. MEMBUAT MATRIKS BOBOT SPASIAL 
# =========================================================
#option 1
W <- poly2nb(Bandung, queen=TRUE) ;W #Get W matrix 
WB <- nb2mat(W, style='B', zero.policy = TRUE) ;WB
Wls <- nb2listw(W, zero.policy = TRUE) ;Wls
Ws <- as(as_dgRMatrix_listw(Wls), "CsparseMatrix") ;Ws
Wm <- as.matrix(Ws)
plot(Bandung, axes=T, col="lightyellow")
plot(Wls, COOR, add=T, pch = 19, col="blue", lwd=1) #Uji autokorelasi spasial

# =========================================================
## 6. UJI AUTOKORELASI SPASIAL 
# =========================================================
# Option 1
Moran14<-moran.test(data$JK[1:30],listw=Wls) ; Moran14
Moran15<-moran.test(data$JK[31:60],listw=Wls) ; Moran15
Moran16<-moran.test(data$JK[61:90],listw=Wls) ; Moran16

# =========================================================
## 7. PENAKSIRAN MENGGUNAKAN SPATIO-TEMPORAL CAR MODEL
# =========================================================
data$Tahun <- as.factor(data$Tahun)
data$id <- as.factor(data$Id)

# Estimasi Model 
# Model1
Model1<-JK~ASI+GIZI
ModelRun1<-
  inla(Model1,family="poisson",data=data,E=Ei,control.predictor=list(compute=TRUE,link=1),
       control.compute=list(dic=TRUE,cpo=TRUE))
summary(ModelRun1)

# Model2
Model2<-JK~ASI+GIZI+f(Id,model="besag",graph=Wm)
ModelRun2<-
  inla(Model2,family="poisson",data=data,E=Ei,control.predictor=list(compute=TRUE,link=1),
       control.compute=list(dic=TRUE,cpo=TRUE))
summary(ModelRun2)

# Model3
Model3<-JK~ASI+GIZI+f(Tahun,model="rw2")
ModelRun3<-
  inla(Model3,family="poisson",data=data,E=Ei,control.predictor=list(compute=TRUE,link=1),
       control.compute=list(dic=TRUE,cpo=TRUE))
summary(ModelRun3)

# Model4
Model4<-JK~ASI+GIZI+f(Id,model="besag",graph=Wm)+f(Tahun,model="rw2")
ModelRun4<-
  inla(Model4,family="poisson",data=data,E=Ei,control.predictor=list(compute=TRUE,link=1),
       control.compute=list(dic=TRUE,cpo=TRUE))
summary(ModelRun4)

# Model5
Tahun1 <- Tahun
Model5<-JK~ASI+GIZI+f(Id,model="besag",graph=Wm)+f(Tahun,model="rw2"
  )+f(Tahun1,model="rw2",group=Id,control.group=list(model="iid"))
ModelRun5<-
  inla(Model5,family="poisson",data=data,E=Ei,control.predictor=list(compute=TRUE,link=1),
       control.compute=list(dic=TRUE,cpo=TRUE))
summary(ModelRun5)

# Menghitung Nilai Risiko Relatif 
ST14<-ModelRun5$summary.fitted.values$mean[1:30]
ST15<-ModelRun5$summary.fitted.values$mean[31:60]
ST16<-ModelRun5$summary.fitted.values$mean[61:90]
ST=c(ST14,ST15,ST16)
data.frame(Kecamatan[1:30],ST14,ST15,ST16)

# Menghitung Standard Error 
SD.ST14=ModelRun5$summary.fitted.values[,2][1:30]
SD.ST15=ModelRun5$summary.fitted.values[,2][31:60]
SD.ST16=ModelRun5$summary.fitted.values[,2][61:90]

# Membuat Boxplot Risiko Relatif 
boxplot(ST14,ST15,ST16, main = "Boxplot ST",
        names = c("ST14","ST15","ST16"))

# =========================================================
# 8. MENGHITUNG KECOCOKAN MODEL
# =========================================================
Y1=ModelRun1$summary.fitted.values$mean*Ei
Rsquare1=sum((Y1-mean(JK))^2)/sum((JK-mean(JK))^2)
Y2=ModelRun2$summary.fitted.values$mean*Ei
Rsquare2=sum((Y2-mean(JK))^2)/sum((JK-mean(JK))^2)
Y3=ModelRun3$summary.fitted.values$mean*Ei
Rsquare3=sum((Y3-mean(JK))^2)/sum((JK-mean(JK))^2)
Y4=ModelRun4$summary.fitted.values$mean*Ei
Rsquare4=sum((Y4-mean(JK))^2)/sum((JK-mean(JK))^2)
Y5=ModelRun5$summary.fitted.values$mean*Ei
Rsquare5=sum((Y5-mean(JK))^2)/sum((JK-mean(JK))^2)

DIC=c(ModelRun1$dic$dic,ModelRun2$dic$dic,
      ModelRun3$dic$dic,ModelRun4$dic$dic, ModelRun5$dic$dic)
Rsquare=c(Rsquare1,Rsquare2,Rsquare3,Rsquare4,
          Rsquare5)
data.frame(DIC,Rsquare)


# =========================================================
# 9. PERBANDINGAN RR UNTUK SMR DAN ST
# =========================================================
par(mfrow=c(1,3))
boxplot(SMR14,ST14, main="Boxplot RR SMR dan ST 2014",names=c("RR SMR","RR ST"))
boxplot(SMR15,ST15, main="Boxplot RR SMR dan ST 2015",names=c("RR SMR","RR ST"))
boxplot(SMR16,ST16, main="Boxplot RR SMR dan ST 2016",names=c("RR SMR","RR ST"))


# =========================================================
# . PERBANDINGAN SE UNTUK SMR DAN ST
# =========================================================
par(mfrow=c(1,3))
boxplot(sd.SMR14,SD.ST14, main="Boxplot SE SMR dan ST 2014",names=c("RR SMR","RR ST"))
boxplot(sd.SMR15,SD.ST15, main="Boxplot SE SMR dan ST 2015",names=c("RR SMR","RR ST"))
boxplot(sd.SMR16,SD.ST16, main="Boxplot SE SMR dan ST 2016",names=c("RR SMR","RR ST"))

# =========================================================
# 10. MEMBUAT PETA PERSEBARAN
# =========================================================
library(maptools)
library(spdep)
library(INLA)
library(sp)
library(RColorBrewer)
library(lattice)
library(gstat)
library(raster)
require(splancs)
library(ggplot2)
library(dplyr)
library(tidyr)
library(brinla)
library(sarima)
library(extrafont) 
library(ggsn)
Indo <- readRDS("/Users/ASUS/Downloads/gadm36_IDN_3_sp.rds")
head(Indo,10)
Bandung <- Indo[Indo$NAME_2 == "Kota Bandung",]
head(Bandung,10)
Kec <- Bandung$NAME_3 
row.names(Bandung) <- Kec
COOR <- coordinates(Bandung)

plot(Bandung, axes=T, col="gray90")
text(COOR[,1],COOR[,2], row.names(Bandung), col="black", cex=0.5, pos=1)
points(COOR[,1],COOR[,2], pch=19, cex=0.5,col="blue")

#Input Data Resiko Relatif
dataST <- data.frame(Id[1:30],Kecamatan[1:30],ST14,ST15,ST16)
attach(dataST)
dataST$id <- dataST$Id.1.30.

dataST$x<-COOR[,1]
dataST$y<-COOR[,2]

#Map Jabar
Bandung1 <- Bandung
Bandung1$ID <- c(1:30) 
Bandung1 <- fortify(Bandung1,region="ID") ;Bandung1

BandungMap <- fortify(Bandung1, regions="ID") ;BandungMap
BandungMap1 <- merge(BandungMap, dataST, by="id", all.x=TRUE)
BandungMap1 <- BandungMap1[order(BandungMap1$order), ] 

#Plot Mapping
p2014 <- ggplot(BandungMap1) +
  geom_polygon( aes(fill = ST14 , x = long,  y = lat, group=group), color="gray47", size=0.1) +  
  scale_fill_gradientn(colours = c("darkolivegreen1", "darkolivegreen2", "darkolivegreen3", "darkolivegreen4","darkolivegreen"), breaks=c(0.7,1.4,2.1,2.8,3.5)) + 
  theme_bw()+ ylab("Northing")+xlab("Easting")+   theme(axis.text.x = element_text(angle = 90))+
  theme(legend.position = "bottom") +  theme(axis.title.x=element_blank(),  axis.text.x=element_blank(),  axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(),  axis.text.y=element_blank(),  axis.ticks.y=element_blank())+ 
  geom_text(data = BandungMap1,aes(x=x, y=y, label = Kecamatan.1.30.), size=2, color="black")+ 
  labs(fill = "Taksiran Risiko Relatif Kasus Pneumonia pada Balita di Kota Bandung Tahun 2014")+ theme(legend.position="bottom", text = element_text(size=14))+
  theme(text = element_text(size=14))  +   guides(fill = guide_colorbar(title.position = "left", title.vjust = 1,  
                                                                        frame.colour = "black",   barwidth = 20,  barheight = 1.5))
p2014

#Plot Mapping
p2015 <- ggplot(BandungMap1) +
  geom_polygon( aes(fill = ST15 , x = long,  y = lat, group=group), color="gray47", size=0.1) +  
  scale_fill_gradientn(colours = c("darkolivegreen1", "darkolivegreen2", "darkolivegreen3", "darkolivegreen4","darkolivegreen"), breaks=c(0.8,1.6,2.4,3.2,4)) + 
  theme_bw()+ ylab("Northing")+xlab("Easting")+   theme(axis.text.x = element_text(angle = 90))+
  theme(legend.position = "bottom") +  theme(axis.title.x=element_blank(),  axis.text.x=element_blank(),  axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(),  axis.text.y=element_blank(),  axis.ticks.y=element_blank())+ 
  geom_text(data = BandungMap1,aes(x=x, y=y, label = Kecamatan.1.30.), size=2, color="black")+ 
  labs(fill = "Taksiran Risiko Relatif Kasus Pneumonia pada Balita di Kota Bandung Tahun 2015")+ theme(legend.position="bottom", text = element_text(size=14))+
  theme(text = element_text(size=14))  +   guides(fill = guide_colorbar(title.position = "left", title.vjust = 1,  
                                                                        frame.colour = "black",   barwidth = 10,  barheight = 1.5))
p2015

#Plot Mapping
p2016 <- ggplot(BandungMap1) +
  geom_polygon( aes(fill = ST16 , x = long,  y = lat, group=group), color="gray47", size=0.1) +  
  scale_fill_gradientn(colours = c("darkolivegreen1", "darkolivegreen2", "darkolivegreen3", "darkolivegreen4","darkolivegreen"), breaks=c(0.7,1.4,2.1,2.8,3.5)) + 
  theme_bw()+ ylab("Northing")+xlab("Easting")+   theme(axis.text.x = element_text(angle = 90))+
  theme(legend.position = "bottom") +  theme(axis.title.x=element_blank(),  axis.text.x=element_blank(),  axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(),  axis.text.y=element_blank(),  axis.ticks.y=element_blank())+ 
  geom_text(data = BandungMap1,aes(x=x, y=y, label = Kecamatan.1.30.), size=2, color="black")+ 
  labs(fill = "Taksiran Risiko Relatif Kasus Pneumonia pada Balita di Kota Bandung Tahun 2016")+ theme(legend.position="bottom", text = element_text(size=14))+
  theme(text = element_text(size=14))  +   guides(fill = guide_colorbar(title.position = "left", title.vjust = 1,  
                                                                        frame.colour = "black",   barwidth = 20,  barheight = 1.5))
p2016
