
library(geoR)
library(maps)
library(mapdata)
library(readxl)
library(rgdal)

diredatos = '/home/danielgp/Dropbox/IMARPE/DATA (copy)/pelagicos/1909-11/'
cruce     = '1909-11'
nombredata= 'Bitacora de celdas 120 kHz Olaya + Flores + Humboldt.csv'

carpeta = '/home/danielgp/Dropbox/IMARPE/Biomasatest/'

## Llamando librerias
source('/home/danielgp/Dropbox/codigos/CodigosVariosR/CodigoInterpolaKrig.R')
source(paste0(carpeta,'codigos/AsignaAreaISO.R'))
source(paste0(carpeta,'codigos/ComputeBiomass_Geo.R'))
source(paste0(carpeta,'codigos/ComputeBiomass.R'))

source('/home/danielgp/Documents/BajaDatos/Codigos/privado/FigImageVar.R')
## matriz para resultados

datos <- read.csv(paste0(diredatos,nombredata))

datos = datos[datos$Lon_M < -60,]

indcero = which(datos$ANC == 0)
indnocero = which(datos$ANC > 0)

x11()
plot(datos$Lon_M,datos$Lat_M, cex= 0.1, col = 'gray')
points(datos$Lon_M[indnocero],datos$Lat_M[indnocero], 
       cex= 1.5*datos$ANC[indnocero]/30000, col = 'red')
map('worldHires', add = T, col = 'khaki1', fill = T)
box()

## Borde de crucero 1902-03
bdr = locator(60, type='l')
bdr = cbind(bdr$x, bdr$y)
bdr = rbind(bdr, bdr[1,])
# # 
write.csv(bdr,paste0(diredatos,'Borde/Borde_190911_V1.csv'),row.names = F)
bdr = read.csv(paste0(diredatos,'Borde/Borde_190911_V1.csv'))

## 
datosposi = datos[indnocero,]
# datosposi$lanc = (datosposi$ANC)

# rm(datos)

#####################################
datos.geo = as.geodata(datosposi, coords.col = c('Lon_M','Lat_M'), data.col = 95)

x11()
plot(datos.geo, lambda=0)

lagvect = seq(from = 0,to = 3,by = 0.1)
varemp = variog(datos.geo, lambda = 0, uvec = lagvect)

x11()
plot(varemp)

var.fit = variofit(varemp, ini.cov.pars = c(6, 1), cov.model = 'exp', nugget = 3, fix.nugget = F)
lines(var.fit)
# var.fit$cov.pars[2] = 0.285

## Llamando a la grilla
gr = read.csv('/home/danielgp/Dropbox/IMARPE/Cr1909-11/Resultado/Prob_V1.csv')
## grilla
# xlim = c(-84, -70)
# ylim = c(-18, -3)
# 
# x = seq(from = xlim[1],to = xlim[2], by = 2/60)
# y = seq(from = ylim[1],to = ylim[2], by = 2/60)

XY = cbind(gr$x,gr$y)

latlimit = c(-3, -19)
parlat   = 10
overlat  = 0.15
lambda   = 0

abase = CodigoInterpolaKrig(datos.geo, datos, var.fit, XY, lambda, latlimit, parlat, overlat)
a = apply(abase, 1, mean, na.rm = T)

gr$pred = a

# gr$FN = gr$z*gr$pred

gr$PA = gr$z
gr$PA[which(gr$z < 0.7)] = 0
# gr$PA[which(gr$z < 0.15)] = 0



inddentro = point.in.polygon(gr$x, gr$y, bdr[,1], bdr[,2])
inddentro = which(inddentro == 0)

gr$FN[inddentro] = NaN
gr$z[inddentro] = NaN
gr$pred[inddentro] = NaN

gr$FN = gr$PA*gr$pred

gr$FN[which(gr$FN < 150 )] = 0 # 130
am = matrix(gr$z, ncol = length(unique(XY[,2])), nrow = length(unique(XY[,1])))

# an = rasterFromXYZ(gr)

xlim = c(-83,-70)
ylim = c(-18.5,-3.5)
xby = yby = 2

xlab = seq(xlim[1],xlim[2], by = xby)
xlabnew = as.character(abs(xlab))
xlabnew[which(xlab/round(xlab) == 1)] = paste0(abs(xlab[which(xlab/round(xlab) ==1)]),'.0')

ylab = seq(ylim[1],ylim[2], by = yby)
ylabnew = as.character(abs(ylab))
ylabnew[which(ylab/round(ylab) == 1)] = paste0(abs(ylab[which(ylab/round(ylab) ==1)]),'.0')



######
adfM = matrix(gr$FN, ncol = length(unique(gr$y)) )
# adfM = adfM[, length(unique(gr$y)) :1]
library(fields)
# image.plot(unique(gr$x), (unique(gr$y)), log10(adfM+1) )
# map('worldHires', add = T, col = 'khaki1', fill = T)
# box()

png('/home/danielgp/Dropbox/IMARPE/Cr1909-11/MapaGeo.png', width = 1050, height = 1200, pointsize = 30)
image.plot(unique(gr$x),unique(gr$y), log10(adfM+1),axes = F, xlab = '', ylab = '', 
           xlim = xlim, ylim = ylim, 
           legend.lab = 'log(NASC anch)')
map('worldHires', fill = T, add = T, col = 'khaki1')

axis(side = 1, at = xlab, cex.axis = 0.9, labels = paste0(xlabnew, '°W'), font = 2, cex.axis = 1)
axis(side = 2, at = ylab, las = 2, cex.axis = 0.9, labels = paste0(ylabnew, '°S'), 
     font = 2, cex.axis = 1)
box()
dev.off()

# , levels = c(0,10,100,50000)

YlOrBr <- c('cyan3','darkolivegreen4','green','green4',
            'yellow','darkorange','red','firebrick4','black')

YlOrBr <- c('white','cyan3','darkolivegreen4','green','green4',
            'yellow','darkorange','red','firebrick4','black')

Cpallete = colorRampPalette(YlOrBr, space = "Lab")

atp = c(0,10,100,250,500,1000,10000,100000)#c(0,10,50,100,500,1000, 20000,100000)#(pretty(adfM))
att = log10(atp+1)

filled.contour(unique(gr$x), (unique(gr$y)),log10(adfM+1),
               color.palette = Cpallete,
               key.axes=axis(4,at=att,labels = atp) , nlevels = length(atp)-1, levels = att ,
               axes = T, add = T, xlim = c(-84,-70), ylim = c(-20,-3), 
               plot.axes = {axis(1);axis(2); 
                 map('worldHires', xlim = c(-84,-70), ylim = c(-20,-3), add = T, fill = T,
                     col = 'khaki1')})


## para blablabla

gr = gr[,c('x','y','FN','z')]

names(gr) = c('Lon_M','Lat_M','NASC_ANC','Zvar')

gr$NASC_ANC[which(is.na(gr$NASC_ANC) == TRUE)] = 0
# gr = gr[which(is.na(gr$) == FALSE),]

## calculando la biomasa
carpetaguarda = paste0(carpeta,'Geo/')
basePaso = AsignaAreaISO(carpeta,carpetaguarda,gr)

datos = read.csv(paste0(carpeta,'Geo/datosparabio.csv'))
datos = datos[which(is.na(datos$ISOcode) == FALSE),]

datalm = read.csv(paste0(carpeta,'lm/','Crucero ',cruce,'_lm.csv'))#'lm/crucero 1802-04.csv')

nombresave = 'aTot.csv'
recurso    = 'Anchoveta'
a = ComputeBiomass(datos,datalm,carpeta,carpetaguarda,nombresave,recurso)

sum(a$Biom)

#  9591961

png('MapaFinal.png')
image(z,col=jet.colors(50))
image(z1,add=TRUE)

map('world',add=T, fill = T, col = 'khaki1')  
legend.krige(c(-90,-88),c(-16,-8),
             z$predict,vertical=T,col=jet.colors(50))

dev.off()
# 

# 6278







# YlOrBr <- c("#FFFFD4", "#FED98E", "#FE9929", "#D95F0E", "#993404")
YlOrBr <- c('white',"skyblue", 'yellow','red',"black")
atp = pretty(am)
att = log10(atp+1)
filled.contour(log10(am),
               color.palette = colorRampPalette(YlOrBr, space = "Lab"),
               key.axes=axis(4,at=att,labels = atp) )




