library(readr)
library(scales)
library(ggplot2)
library(ggcorrplot)
library(dplyr)
library("depmixS4")
library("lubridate")
library(ggbiplot)
library(robustbase)
library("devtools")
library("factoextra")
library(dplyr)
library(magrittr)  # to use piping %>%
library(ggplot2)   # for ploting
library(MASS)      # to calculate the pseudo-inverse of a matrix
library(caret)     # to center our data by subtracting its mean
library(reshape)  # for data manipulation



#part 1 PCA:
dataset <- read_csv("/Users/ericw/Downloads/Project_Dataset.txt",
                    col_types = cols(
                      Date = col_date("%d/%m/%Y"),
                      Time = col_character(),
                      Global_active_power = col_double(),
                      Global_reactive_power = col_double(),
                      Voltage = col_double(),
                      Global_intensity = col_double(),
                      Sub_metering_1 = col_double(),
                      Sub_metering_2 = col_double(),
                      Sub_metering_3 = col_double()
                    )
)


dataset <- na.omit(dataset)
colSums(is.na(dataset))
num_data <- dataset[,3:9]
head(num_data)

num_data <- num_data %>% mutate_all(~(scale(.) %>% as.vector))

corr_matrix <- cor(num_data)
ggcorrplot(corr_matrix)



data.pca <- princomp(corr_matrix)
summary(data.pca)

data.pca$loadings[, 1:2]

fviz_eig(data.pca, addlabels = TRUE)
ggbiplot(data.pca, addlabels = TRUE)

#none ggbiplot version
fviz_pca_var(data.pca, col.var = "black")

fviz_cos2(data.pca, choice = "var", axes = 1:2)

fviz_pca_var(data.pca, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)


# chosen response variables are: Voltage, global_intensity
# Can also include subm3 and global active power.


###################PART 2####################

#choosing a time window:
dataset$Date <-as.Date(dataset$Date,format = "%d/%m/%Y")
dataset$weekday <- wday(dataset$Date,label = TRUE)
dataset_day <- dataset[dataset$weekday == "Mon",]
dataset_filtered <- dataset_day[dataset_day$Time >= "20:00:00" & dataset_day$Time <= "23:30:00" ,]

colNames <- c('Voltage','Date',"Time")
row <- dataset_filtered[,colNames]
row$Global_active_power = scale(row$Voltage)
row[is.na(row)] <- 0

MonGroup = 0
graph_example <- cbind(row, MonGroup)
graph_example$MonGroup[which(graph_example$Date == "2006-12-18")] <- 1
graph_example$MonGroup[which(graph_example$Date == "2007-01-08")] <- 2
graph_example$MonGroup[which(graph_example$Date == "2007-02-12")] <- 3
graph_example$MonGroup[which(graph_example$Date == "2007-03-12")] <- 4
graph_example$MonGroup[which(graph_example$Date == "2007-04-16")] <- 5
graph_example$MonGroup[which(graph_example$Date == "2007-05-14")] <- 6
graph_example$MonGroup[which(graph_example$Date == "2007-06-18")] <- 7
graph_example$MonGroup[which(graph_example$Date == "2007-07-16")] <- 8
graph_example$MonGroup[which(graph_example$Date == "2007-08-20")] <- 9
graph_example$MonGroup[which(graph_example$Date == "2007-09-17")] <- 10
graph_example$MonGroup[which(graph_example$Date == "2007-10-15")] <- 11
graph_example$MonGroup[which(graph_example$Date == "2007-11-19")] <- 12
graph_example$MonGroup[which(graph_example$Date == "2007-12-17")] <- 13


graph_example <- filter(graph_example, MonGroup != 0)


graph_example %>% 
  group_by(MonGroup) %>% 
  ggplot()+geom_line(aes(Time, Voltage, group=MonGroup, color=MonGroup))


#Start traning and testing
dataset <- read_csv("/Users/ericw/Downloads/Project_Dataset.txt",
                    col_types = cols(
                      Date = col_date("%d/%m/%Y"),
                      Time = col_character(),
                      Global_active_power = col_double(),
                      Global_reactive_power = col_double(),
                      Voltage = col_double(),
                      Global_intensity = col_double(),
                      Sub_metering_1 = col_double(),
                      Sub_metering_2 = col_double(),
                      Sub_metering_3 = col_double()
                    )
)
dataset <- na.omit(dataset)

dataset$weekday <- wday(dataset$Date,label = TRUE)
dataset_day <- dataset[dataset$weekday == "Mon",]

dataset_filtered <- dataset_day[dataset_day$Time >= "20:00:00" & dataset_day$Time <= "23:30:00",]
colNames <- c('Voltage','Global_intensity','Date',"Time")
row <- dataset_filtered[,colNames]
row <-subset(row, select = -c(Date,Time) )
row <- row%>% mutate_all(~(scale(.) %>% as.vector))


train_set <- row[(1:26164),]
test_set <- row[(26165:32705),]


set.seed(1)

states<-c(5,7,9,11,13,15,17,19,21)
models <- vector()
BIC <- vector()
logs <-vector()

for(j in states){
  
  mod1<- depmix(list(train_set$Voltage ~1,train_set $Global_intensity~1), 
                data=train_set , nstates=j,
                family=list(gaussian(),gaussian()), 
                ntimes=length(train_set$Voltage))
  fm1 <- fit(mod1)
  print(j)
  print(fm1)
  models <- append(models, fm1);
  BIC <- append(BIC, BIC(fm1))
  logs <- append(logs, logLik(fm1))
}

# Apply state 19
mod1<- depmix(list(train_set$Voltage ~1,train_set $Global_intensity~1), 
              data=train_set , nstates=19,
              family=list(gaussian(),gaussian()), 
              ntimes=length(train_set$Voltage))
fm1 <- fit(mod1)
print(fm1)


ggplot()+
  geom_line(aes(x = states, y = BIC,color = "darkred"))+
  geom_line(aes(x = states, y = logs,color = "steelblue")) +
  scale_colour_manual (labels = c("BIC","Log likelihood"),values = c("darkred","steelblue"))+
  labs(title = "BIC and log likelyhood", y = "Unit values")

table <- data.frame(states,BIC,logs)

print(table)

mod1<- depmix(list(train_set$Voltage ~1,train_set $Global_intensity~1), 
              data=train_set , nstates=19,
              family=list(gaussian(),gaussian()), 
              ntimes=length(train_set$Voltage))
fm1 <- fit(mod1)


mod2<- depmix(list(test_set$Voltage ~1,test_set $Global_intensity~1), 
              data=test_set , nstates= 19,
              family=list(gaussian(), gaussian()), 
              ntimes=length(test_set$Voltage))

mod2<- setpars(mod2, getpars(fm1))

fb<-forwardbackward(mod2)
fb$logLike


##########################PART 3#####################3

a1 <- read_csv("/Users/ericw/Downloads/Dataset_with_Anomalies_1.txt",
                    col_types = cols(
                      Date = col_date("%d/%m/%Y"),
                      Time = col_character(),
                      Global_active_power = col_double(),
                      Global_reactive_power = col_double(),
                      Voltage = col_double(),
                      Global_intensity = col_double(),
                      Sub_metering_1 = col_double(),
                      Sub_metering_2 = col_double(),
                      Sub_metering_3 = col_double()
                    )
)

a2 <- read_csv("/Users/ericw/Downloads/Dataset_with_Anomalies_2.txt",
               col_types = cols(
                 Date = col_date("%d/%m/%Y"),
                 Time = col_character(),
                 Global_active_power = col_double(),
                 Global_reactive_power = col_double(),
                 Voltage = col_double(),
                 Global_intensity = col_double(),
                 Sub_metering_1 = col_double(),
                 Sub_metering_2 = col_double(),
                 Sub_metering_3 = col_double()
               )
)

a3 <- read_csv("/Users/ericw/Downloads/Dataset_with_Anomalies_3.txt",
               col_types = cols(
                 Date = col_date("%d/%m/%Y"),
                 Time = col_character(),
                 Global_active_power = col_double(),
                 Global_reactive_power = col_double(),
                 Voltage = col_double(),
                 Global_intensity = col_double(),
                 Sub_metering_1 = col_double(),
                 Sub_metering_2 = col_double(),
                 Sub_metering_3 = col_double()
               )
)



a1 <- na.omit(a1)
a2 <- na.omit(a2)
a3 <- na.omit(a3)



date <- a1$Date
time <- a1$Time

a1$weekday <- wday(a1$Date,label = TRUE)
a1_day <- a1[a1$weekday == "Mon",]

a1_filtered <- a1_day[a1_day$Time >= "20:00:00" & a1_day$Time <= "23:30:00",]
colNames <- c('Voltage','Global_intensity','Date',"Time")
a1 <- a1_filtered[,colNames]
a1 <-subset(a1, select = -c(Date,Time) )
a1 <- a1%>% mutate_all(~(scale(.) %>% as.vector))


date <- a2$Date
time <- a2$Time
a2$weekday <- wday(a2$Date,label = TRUE)
a2_day <- a2[a2$weekday == "Mon",]

a2_filtered <- a2_day[a2_day$Time >= "20:00:00" & a2_day$Time <= "23:30:00",]
colNames <- c('Voltage','Global_intensity','Date',"Time")
a2 <- a2_filtered[,colNames]
a2 <-subset(a2, select = -c(Date,Time) )
a2 <- a2%>% mutate_all(~(scale(.) %>% as.vector))

date <- a3$Date
time <- a3$Time
a3$weekday <- wday(a3$Date,label = TRUE)
a3_day <- a3[a3$weekday == "Mon",]

a3_filtered <- a3_day[a3_day$Time >= "20:00:00" & a3_day$Time <= "23:30:00",]
colNames <- c('Voltage','Global_intensity','Date',"Time")
a3 <- a3_filtered[,colNames]
a3 <-subset(a3, select = -c(Date,Time) )
a3 <- a3%>% mutate_all(~(scale(.) %>% as.vector))

set.seed(1)
moda1<- depmix(list(a1$Voltage~1,a1$Global_intensity~1), 
              data=a1 , nstates= 19,
              family=list(gaussian(), gaussian()), 
              ntimes=length(a1$Voltage))
moda1<- setpars(moda1, getpars(fm1))
fb<-forwardbackward(moda1)
fb$logLike

moda2<- depmix(list(a2$Voltage~1,a2$Global_intensity~1), 
              data=a2 , nstates= 19,
              family=list(gaussian(), gaussian()), 
              ntimes=length(a2$Voltage))
moda2<- setpars(moda2, getpars(fm1))
fb<-forwardbackward(moda2)
fb$logLike

moda3<- depmix(list(a3$Voltage~1,a3$Global_intensity~1), 
              data=a3 , nstates= 19,
              family=list(gaussian(), gaussian()), 
              ntimes=length(a3$Voltage))
moda3<- setpars(moda3, getpars(fm1))
fb<-forwardbackward(moda3)
fb$logLike



#########Plotting data

#change data.frame a1,a2,a3
X=as.data.frame(a3)
names(X)=c("Voltage", "Global_intensity")
XX=melt(X)
XX%>%ggplot(aes(x=value,fill=variable, color=variable))+
  geom_density(alpha = 0.3)+ggtitle('Distibution of X')

X%>%ggplot(aes(x=`Voltage`,y=`Global_intensity`))+
  geom_point(color='darkgreen')

X=as.data.frame(a1)
names(X)=c("Voltage", "Global_intensity")
XX=melt(X)
XX%>%ggplot(aes(x=value,fill=variable, color=variable))+
geom_density(alpha = 0.3)+ggtitle('Distibution of X')

X%>%ggplot(aes(x=`Voltage`,y=`Global_intensity`))+
  geom_point(color='blue')