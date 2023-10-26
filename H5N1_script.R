## https://becarioprecario.bitbucket.io/spde-gitbook/ch-spacetime.html

library(tidyverse)
library(INLA)
library(ggpubr)

## DATA IMPORT
H5N1_JUN23 <- read.csv("https://raw.githubusercontent.com/ProfNascimento/H5N1/main/H5N1_JUN23.csv")
names(H5N1_JUN23)
H5N1_JUN23$TIME=as.Date(H5N1_JUN23$TIME, tryFormats = c("%Y-%m-%d"))

## WEEKLY GROWTH
H5N1_JUN23 %>% group_by(MONTH = lubridate::floor_date(TIME, "week") ) %>% 
  summarise(CUSUM=sum(RECORD)) %>% ggplot(aes(x=MONTH,y=CUSUM)) + geom_line() + 
  geom_point() + ylab("TOTAL OF ANIMALS (frequency)") + xlab("TIME PERIOD (weekly)")

## MONTHLY DYNAMIC
H5N1_JUN23 %>% group_by(MONTH = lubridate::floor_date(TIME, "month") ) %>% 
  ggplot(aes(x=as.factor(MONTH),y=RECORD)) + geom_boxplot(outlier.alpha = 0.3) +
  ylab("NUMBER OF ANIMALS (frequency per field' visits)") + xlab("TIME PERIOD (monthly)")

## TOTAL PER REGION
H5N1_JUN23 %>% group_by(STATE) %>% 
  summarise(TOTAL=sum(RECORD), MEAN=mean(RECORD)) %>% arrange(desc(TOTAL))

## PER SPECIES TYPE
library(stringr)
H5N1_week = H5N1_JUN23 %>% 
  mutate(WEEK = lubridate::week(TIME),
         MONTH = lubridate::month(TIME),
         CUSUM=sum(RECORD),
         type2 = str_extract(SPECIES, "\\b[A-Z]{2}"),
         type3 = str_sub(SPECIES, 6,))

H5N1_week$type2=as.factor(H5N1_week$type2)
levels(H5N1_week$type2)=c("Birds","Cetaceans","Undef","MU","PI","QU")

tapply(H5N1_week$RECORD, H5N1_week$type2, summary) # SUMMARY
tapply(H5N1_week$RECORD, H5N1_week$type2, sum) # TOTAL

H5N1_week %>% group_by(STATE, type3) %>% summarise(TOTAL=sum(RECORD)) %>% arrange(desc(TOTAL))

sort(tapply(H5N1_week$RECORD, as.factor(H5N1_week$type3),sum),decreasing = TRUE) # TOTAL

###########################################################################
# Model TMS data using SPDE model
# Import data
coords = data.frame(y=H5N1_week$LAT,x=H5N1_week$LON)
k=max(H5N1_week$MONTH)

#
bnd = inla.nonconvex.hull(as.matrix(coords))
mesh = inla.mesh.2d(boundary = bnd, max.edge = c(2,2)) #, boundary=bnd,
plot(mesh)
# SPDE model
spde <- inla.spde2.pcmatern(mesh = mesh,
                            prior.range = c(0.5, 0.01), # P(range < 0.05) = 0.01
                            prior.sigma = c(1, 0.01))   # P(sigma > 1) = 0.01

iset <- inla.spde.make.index(name="spatial.field", n.spde = spde$n.spde, n.group = k)

A <- inla.spde.make.A(mesh = mesh,
                      loc = as.matrix(coords),
                      group = H5N1_week$MONTH)

H5N1_week$STATE=as.factor(H5N1_week$STATE)
H5N1_week$type2=as.factor(H5N1_week$type2)

sdat2 <- inla.stack(
  data = list(y = H5N1_week$RECORD),
  A = list(A,1),
  effects = list(c(iset,list(Intercept=0)),
                 list(Type=H5N1_week$type2,
                      ID.Week = H5N1_week$WEEK,
                      State = H5N1_week$STATE,
                      Month = H5N1_week$MONTH)),
  tag = 'stk')

h.spec <- list(rho = list(prior = 'pc.cor1', param=c(0, 0.9)))

## MODEL 1
formula = y ~ -1 +
  f(spatial.field, model = spde, 
    control.group = list(model = 'ar1', hyper = h.spec)) +
  f(ID.Week, model = "ar1", 
    hyper = list(prec=list(param=c(10,100)))) +
  f(State, Month, model = "rw1")

output <- inla(formula = formula,
               family = "tweedie",
               data = inla.stack.data(sdat2),
               quantiles = c(0.1, 0.5, 0.9),
               control.predictor = list(compute = TRUE,A=inla.stack.A(sdat2)),
               control.compute=list(config = TRUE,cpo=TRUE, dic=TRUE, waic=TRUE,
                                    return.marginals.predictor=TRUE)
)

summary(output)
INLAutils::plot_inla_residuals(output,H5N1_week$RECORD)

## MODEL 2
formula2 = y ~ -1 + Type*State +
  f(spatial.field, model = spde, group = spatial.field.group, 
    control.group = list(model = 'ar1', hyper = h.spec)) +
  f(ID.Week, model = "ar1", 
    hyper = list(prec=list(param=c(10,100)))) +
  f(State, Month, model = "rw1")

output3 <- inla(formula = formula2,
                family = "tweedie",
                data = inla.stack.data(sdat2),
                quantiles = c(0.1, 0.5, 0.9),
                control.predictor = list(compute = TRUE,A=inla.stack.A(sdat2)),
                control.compute=list(config = TRUE,cpo=TRUE, dic=TRUE, waic=TRUE,
                                     return.marginals.predictor=TRUE)
)

summary(output3)
INLAutils::plot_inla_residuals(output3,H5N1_week$RECORD)


#
idat <- inla.stack.index(sdat2, 'stk')$data
cor(H5N1_week$RECORD, output3$summary.linear.predictor$mean[idat])

### Spatial plot
stepsize <- 0.4
nxy <- round(
  c(diff(range( coords[, 1])),
    diff(range( coords[, 2]))) / stepsize)

projgrid <- inla.mesh.projector(
  mesh, xlim = range(coords[, 1]),
  ylim = range(coords[, 2]), dims = nxy)

xmean <- list()
for(j in 1:k){
  xmean[[j]] = inla.mesh.project(
    projgrid, output3$summary.random$spatial.field$mean[iset$spatial.field.group==j])
}

library(splancs)
xy.in = inout(projgrid$lattice$loc, cbind(PRborder[,1],PRborder[,2]))


###################
library(sf)
library(rnaturalearthhires)
library(dplyr)

chile <- rnaturalearthhires::countries10 %>%
  st_as_sf() %>%
  filter(SOVEREIGNT %in% c("Chile")) %>% st_crop(
    xmin=-90, xmax=-30, ymin=-60, ymax=-10)

library(rgdal)
cord1 = data.frame(x=ResultDatframe$x,y=ResultDatframe$y)

coordinates(cord1) <- ~x+y
class(cord1)
proj4string(cord1) <- CRS("+proj=longlat +zone=19 +south +datum=WGS84")
utm2 <- spTransform(cord1,CRS("+proj=longlat +datum=WGS84"))
ss=as.data.frame(coordinates(utm2))

## DRAWING THE OCEAN IN THE PLOT
ocean <- rnaturalearth::ne_download(
  type = "ocean",
  category = "physical",
  returnclass = "sf")

#### UNLIST SPATIO ESTIMATION (PER MONTH --GROUP--)
ResultDatframe <- data.frame(df,
                             mean_s1=as.vector(xmean[[1]]),
                             mean_s2=as.vector(xmean[[2]]),
                             mean_s3=as.vector(xmean[[3]]),
                             mean_s4=as.vector(xmean[[4]]),
                             mean_s5=as.vector(xmean[[5]]),
                             mean_s6=as.vector(xmean[[6]]))

#Chile_map <- map_data("world", region="Chile")
G1=ggplot() +
  geom_raster(data = ss*(-1), aes(x=y,y=x,fill=ResultDatframe$mean_s1))+
  scale_fill_continuous(type = "viridis",limits=c(-3,3)) + labs(fill = "log(rate)") +
  xlab("") + ylab("") + theme(axis.text.x = element_text(angle = 45, hjust=1),
                              legend.position = "none") +
  # OCEAN DRAW  
  geom_sf(data = ocean,
          color = NA,
          fill = "white",
          size = 0.2) +
  # INLAND DRAW
  geom_sf(data = chile,
          color = "white",
          fill = NA,
          size = 0.2) +
  # CROP the whole thing to size
  coord_sf(xlim = c(-77.5, -67.5),
           ylim = c(-16, -55))

G2=ggplot() +
  geom_raster(data = ss*(-1), aes(x=y,y=x,fill=ResultDatframe$mean_s2))+
  scale_fill_continuous(type = "viridis",limits=c(-3,3)) + labs(fill = "log(rate)") +
  xlab("") + ylab("") + theme(axis.text.x = element_text(angle = 45, hjust=1),
                              legend.position = "none") +
  # OCEAN DRAW  
  geom_sf(data = ocean,
          color = NA,
          fill = "white",
          size = 0.2) +
  # INLAND DRAW
  geom_sf(data = chile,
          color = "white",
          fill = NA,
          size = 0.2) +
  # CROP the whole thing to size
  coord_sf(xlim = c(-77.5, -67.5),
           ylim = c(-16, -55))

G3=ggplot() +
  geom_raster(data = ss*(-1), aes(x=y,y=x,fill=ResultDatframe$mean_s3))+
  scale_fill_continuous(type = "viridis",limits=c(-3,3)) + labs(fill = "log(rate)") +
  xlab("") + ylab("") + theme(axis.text.x = element_text(angle = 45, hjust=1),
                              legend.position = "none") +
  # OCEAN DRAW  
  geom_sf(data = ocean,
          color = NA,
          fill = "white",
          size = 0.2) +
  # INLAND DRAW
  geom_sf(data = chile,
          color = "white",
          fill = NA,
          size = 0.2) +
  # CROP the whole thing to size
  coord_sf(xlim = c(-77.5, -67.5),
           ylim = c(-16, -55))

G4=ggplot() +
  geom_raster(data = ss*(-1), aes(x=y,y=x,fill=ResultDatframe$mean_s4))+
  scale_fill_continuous(type = "viridis",limits=c(-3,3)) + labs(fill = "log(rate)") +
  xlab("") + ylab("") + theme(axis.text.x = element_text(angle = 45, hjust=1),
                              legend.position = "none") +
  # OCEAN DRAW  
  geom_sf(data = ocean,
          color = NA,
          fill = "white",
          size = 0.2) +
  # INLAND DRAW
  geom_sf(data = chile,
          color = "white",
          fill = NA,
          size = 0.2) +
  # CROP the whole thing to size
  coord_sf(xlim = c(-77.5, -67.5),
           ylim = c(-16, -55))

G5=ggplot() +
  geom_raster(data = ss*(-1), aes(x=y,y=x,fill=ResultDatframe$mean_s5))+
  scale_fill_continuous(type = "viridis",limits=c(-3,3)) + labs(fill = "log(rate)") +
  xlab("") + ylab("") + theme(axis.text.x = element_text(angle = 45, hjust=1),
                              legend.position = "none") +
  # OCEAN DRAW  
  geom_sf(data = ocean,
          color = NA,
          fill = "white",
          size = 0.2) +
  # INLAND DRAW
  geom_sf(data = chile,
          color = "white",
          fill = NA,
          size = 0.2) +
  # CROP the whole thing to size
  coord_sf(xlim = c(-77.5, -67.5),
           ylim = c(-16, -55))

G6=ggplot() +
  geom_raster(data = ss*(-1), aes(x=y,y=x,fill=ResultDatframe$mean_s6))+
  scale_fill_continuous(type = "viridis",limits=c(-3,3)) + labs(fill = "log(rate)") +
  xlab("") + ylab("") + theme(axis.text.x = element_text(angle = 45, hjust=1),
                              legend.position = "none") +
  # OCEAN DRAW  
  geom_sf(data = ocean,
          color = NA,
          fill = "white",
          size = 0.2) +
  # INLAND DRAW
  geom_sf(data = chile,
          color = "white",
          fill = NA,
          size = 0.2) +
  # CROP the whole thing to size
  coord_sf(xlim = c(-77.5, -67.5),
           ylim = c(-16, -55))

cowplot::plot_grid(G1,G2,G3,G4,G5,G6,nrow=1)
