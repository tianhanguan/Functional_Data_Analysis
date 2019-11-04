library('fda')

##Step 1: Load the dataset
#Import the dataset
dat = read.csv('Altcoins.csv', header = TRUE)
dim(dat)
head(dat)

#time series plot for all 20 replicates
price = as.matrix(dat[,-1]) #365 rows(obs per replicate), 20 columns (20 replicates)
price = mapply(price, FUN = as.numeric)
price = matrix(data = price, ncol = 20, nrow = 365)
day = dat[, 1]  #time domain values
ts.plot(
  price,
  gpars = list(
    col = rainbow(20),
    xlab = 'Day',
    ylab = 'Price (log scale)',
    lty = c(1:20)
  ),
  main = 'Daily closing price of 20 Altcoins (2017-2018)'
)
########################################################################
#Step 2: Pre-smoothing using Penalized Spline models
#First, use GCV to select the number of basis functions
#Number of basis = (# interior knots) + (order 4 or cubic splines)
#nbasis can go from 4 (the smallest) up to (364+4=368)
rangeval0 = c(1, 365)
gcvlist = NULL
nbasis = seq(4, 100, 1)

for (i in 1:length(nbasis)) {
  basisobj = create.bspline.basis(rangeval0, nbasis[i]) #create the basis
  ys = smooth.basis(day, price, basisobj) #fit the smoothed data
  gcv = ys$gcv
  gcvlist = c(gcvlist, mean(gcv))
}

#Draw the plot - choose K = 50 (that is, 48 knots in total)
GCV = as.data.frame(cbind(nbasis, gcvlist))
library(ggplot2)
ggplot(data = GCV, aes(x = nbasis, y = gcvlist)) +
  geom_point(size = 2, shape = 16) +
  labs(x = "Number of basis", y = "GCV")
nbasis0 = 50

#Construct the basis system using nbasis = 50
nbasis0 = 50
basisobj = create.bspline.basis(rangeval0, nbasis = nbasis0) #create the basis object
plot(basisobj, main = 'Basis system of cubic splines (number of basis = 50)')

#Next, use GCV to choose smoothing parameter lambda
lambda0 = seq(0.1, 40, 0.5)
gcvlist = NULL

for (i in 1:length(lambda0)) {
  #roughness penalty with smoothing parameter lambda
  fdParobj = fdPar(fdobj = basisobj,
                   Lfdobj = 2,
                   lambda = lambda0[i])
  smoothlist = smooth.basis(day, price, fdParobj)
  gcv = smoothlist$gcv
  gcvlist = c(gcvlist, mean(gcv)) #average the gcv for all replicates
}

GCV = as.data.frame(cbind(lambda0, gcvlist))
ggplot(data = GCV, aes(x = lambda0, y = gcvlist)) +
  geom_point(size = 2, shape = 16) +
  labs(x = "Smoothing parameter lambda", y = "GCV")

lambda0 = lambda0[which.min(gcvlist)] #selected lambda value = 8.1

##Now ready to create the bsplie basis object
#Set up spline basis system
basisobj = create.bspline.basis(rangeval0, nbasis0)

#Set up roughness penalty with smoothing parameter
fdParobj = fdPar(fdobj = basisobj,
                 Lfdobj = 2,
                 lambda = lambda0)

#Smooth the data, outputting a list containing various quantities
smoothlist = smooth.basis(day, price, fdParobj)
xfd = smoothlist$fd   #finally create the fd object! (i.e. functional data!)

#Plot the functional data (all 20 replicates)
plot(
  xfd,
  xlab = 'Day',
  ylab = 'Price (log-scale)',
  main = 'Functional Data - Daily closing price of 20 Altcoins (2017-2018)')


##Add the fitted 'mean line' on functional data
plot(
  xfd,
  xlab = 'Day',
  ylab = 'Price (log-scale)',
  main = 'Functional Data (with the estimated mean function) 
  - Daily closing price of 20 Altcoins (2017-2018)')
lines(mean(xfd), lty=1, lwd=4, col=1)

## compare a curve to the data from which it was estimated by using
## the 'plotfit.fd' function (for the first 5 replicates)
plotfit.fd(price[,1], day, xfd[1],
           lty=1, lwd=1,col='deepskyblue', 
           main='Functional data - Altcoin #1')
plotfit.fd(price[,2], day, xfd[2],
           lty=1, lwd=2, col='palegreen3', main='Functional data - Altcoin #2')
plotfit.fd(price[,3], day, xfd[3],
           lty=1, lwd=2, col='peachpuff3', main='Functional data - Altcoin #3')
plotfit.fd(price[,4], day, xfd[4],
           lty=1, lwd=2, col='seagreen3', main='Functional data - Altcoin #4')
plotfit.fd(price[,5], day, xfd[5],
           lty=1, lwd=2, col='tomato3', main='Functional data - Altcoin #5')

#Covariance matrix (MDA)
cormat = cor(price,method='pearson')
library(reshape2)
melted_cormat <- melt(cormat)
head(melted_cormat)
library(ggplot2)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

#Covariance surface (FDA)
# define the variance-covariance bivariate fd object
pricevarbifd <- var.fd(xfd)
# evaluate the variance-covariance surface and plot
pricevarmat <- eval.bifd(day,day,pricevarbifd) #dim = 365*365

#First plot the contour plots
par(bg = "beige")
contour(day, day, pricevarmat, 
        xlab="Day",
        ylab="Day",
        main="Contour plot for the covariance function",
        cex.main=0.8, axes=TRUE,col='deepskyblue',lty = "solid",lwd=3)

#Next plot the 3D covariance surface
par(bg = "beige")
persp(day, day, pricevarmat, theta = 45, phi = 5,
      xlab="Day", ylab="Day", zlab="Covariance",
      border = 'chartreuse4',main="Covariance surface")

#############################################################
#Step 3: Functional PCA!
pricepca <- pca.fd(xfd, nharm=2, fdPar(xfd))
#pricepca2 = varmx.pca.fd(pricepca)
eigenfunc = pricepca$harmonics #coefs = 50 rows, 2 columns (50 basis, 2 PCs)
eigenval = pricepca$values #50 eigenvalues
scores = pricepca$scores #20 rows(replicates), 2 columns (2 PCs)
var = pricepca$varprop #first 2 PCs take 95% of variations

plot.pca.fd(pricepca, cex.main=0.9,xlab='Day',ylab=' ')
plot(pricepca$harmonics)  