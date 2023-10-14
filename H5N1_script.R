
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
         type2 = str_extract(SPECIES, "\\b[A-Z]{2}"),
         type3 = str_sub(SPECIES, 6,)) 

H5N1_week$type2=as.factor(H5N1_week$type2)
levels(H5N1_week$type2)=c("Birds","Cetaceans","Undef","MU","PI","QU")

tapply(H5N1_week$RECORD, H5N1_week$type2, summary) # SUMMARY
tapply(H5N1_week$RECORD, H5N1_week$type2, sum) # TOTAL

H5N1_week %>% group_by(STATE, type3) %>% summarise(TOTAL=sum(RECORD)) %>% arrange(desc(TOTAL))

sort(tapply(H5N1_week$RECORD, as.factor(H5N1_week$type3),sum),decreasing = TRUE) # TOTAL

###########################################################################
W = matrix(c(0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,
             1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
             0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,
             0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,
             0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,
             0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,
             0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,
             0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,
             0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,
             0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,
             0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,
             0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,
             0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,
             1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,0,0,1,1,0,0,0,0,0,0,0),ncol=15,byrow=TRUE)

test1=H5N1_week %>% group_by(Región) %>% summarise(log(mean(RECORD)))
E=unlist(test1[,2])

scaled_x <- c(H5N1_week$week)
X <- model.matrix(~scaled_x)
O <- H5N1_week$RECORD

full_d <- list(n = dim(X)[1],
               p = dim(X)[2],
               X = X,               # design matrix
               y = O,               # observed number of cases
               W = W)               # adjacency matrix

full_fit <- stan(file="/home/diego/Documents/ARTICLEs/H5N1/car_prec.stan",
                 data = full_d,
                 iter = 500,
                 chains = 1)

fitChain=as.matrix(full_fit)
colMeans(fitChain)
