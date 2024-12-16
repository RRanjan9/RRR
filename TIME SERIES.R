library(datasets) 
data=Nile 
plot(data) 
require(graphics) 
data1<-ts(data,start=c(1871,1), end=c(1970,12), frequency=12) 
decomp<-decompose(data1) # to decompose the data into different components 
plot(decomp) 

# (a) Autocorrelation up to order three
Autocorr<-acf(data1,lag.max = 3,plot=FALSE) 
Autocorr 

# (b) Partial Autocorrelation
Partial_Autocorr<-acf(data1,lag.max = 3,type=c("partial"),plot=FALSE) 
Partial_Autocorr

# (c) Plot the autocorrelation and partial autocorrelation functions
par(mfrow = c(2, 1)) # Set layout for two plots in one column
acf(data1, lag.max = 10, main = "Autocorrelation Function")
pacf(data1, main = "Partial Autocorrelation Function")