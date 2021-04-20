library(raster)
library(mgcv)
library(akima)

setwd('/home/danielgp/Dropbox/IMARPE/Cr1909-11')

datos <- read.csv('Bitacora de celdas 120 kHz Olaya + Flores + Humboldt.csv')
datos = datos[datos$Lon_M < -60,]

indcero = which(datos$ANC == 0)
indnocero = which(datos$ANC > 0)

datos$PA = NaN*datos$Lon_M

datos$PA[indcero]   = 0
datos$PA[indnocero] = 1

osm = raster(paste0(getwd(),'/grd/OSM.grd'))
ssm = raster(paste0(getwd(),'/grd/SSM.grd'))
tsm = raster(paste0(getwd(),'/grd/TSM.grd'))

datos$osm = extract(x = osm, cbind(datos$Lon_M,datos$Lat_M))
datos$ssm = extract(x = ssm, cbind(datos$Lon_M,datos$Lat_M))
datos$tsm = extract(x = tsm, cbind(datos$Lon_M,datos$Lat_M))



XY = coordinates(tsm)
XY = data.frame(Lon_M = XY[,1], Lat_M = XY[,2])


lonU = unique(XY$Lon)
latU = unique(XY$Lat)

XY$osm = matrix(osm, ncol = 1)
XY$ssm = matrix(ssm, ncol = 1)
XY$tsm = matrix(tsm, ncol = 1)

indnonan = which(is.na(XY$osm) == FALSE)
XY   = XY[indnonan,]

## Interpola
X = seq(-83, -70, by = 1/60)
Y = seq(-18.5, -3.5, by = 1/60)

XYT = expand.grid(X,Y)

indnonan = which(is.na(XY$osm) == FALSE)
XY   = XY[indnonan,]
zosm = interpp(x=XY$Lon_M, y = XY$Lat_M, z = XY$osm, xo = XYT$Var1, yo = XYT$Var2)

indnonan = which(is.na(XY$ssm) == FALSE)
XY   = XY[indnonan,]
zssm = interpp(x=XY$Lon_M, y = XY$Lat_M, z = XY$ssm, xo = XYT$Var1, yo = XYT$Var2)

indnonan = which(is.na(XY$tsm) == FALSE)
XY   = XY[indnonan,]
ztsm = interpp(x=XY$Lon_M, y = XY$Lat_M, z = XY$tsm, xo = XYT$Var1, yo = XYT$Var2)

XYT$osm = zosm$z
XYT$ssm = zssm$z
XYT$tsm = ztsm$z

names(XYT) = c('Lon_M','Lat_M','osm','ssm','tsm')
# a = rasterFromXYZ(XYT)

##

mod1 = gam(PA ~ s(osm) + s(ssm) + s(tsm) + s(Lon_M) + s(Lon_M,Lat_M), data = datos, family = binomial())
summary(mod1)

val = predict(mod1, XYT, type = 'response')

adf = data.frame(x = XYT$Lon, y = XYT$Lat, z = matrix(val))

# adf$z[which(adf$z < 0)] = 0
# adf$z[which(adf$z > 1)] = 1

write.csv(adf,paste0(getwd(),'/Resultado/Prob_V1.csv'),row.names = F)

a = rasterFromXYZ(adf)
adfM = matrix(adf$z, ncol = length(unique(adf$y)) )
adfM = adfM[,400:1]

image.plot(unique(adf$x), (unique(adf$y)), adfM)


YlOrBr <- c('white',"skyblue", 'yellow','red',"black")
 atp = pretty(adfM)
 att = atp
 filled.contour(unique(adf$x), rev(unique(adf$y)),adfM,
                color.palette = colorRampPalette(YlOrBr, space = "Lab"),
                key.axes=axis(4,at=att,labels = atp) , nlevels = 5)


