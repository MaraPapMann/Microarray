X <- c(8.4395, 8.2963, 7.9803, 7.7198, 8.8918, 10.6271, 10.8426, 11.6703, 10.7754)
Xsum <- sum(X)
PX <- X/Xsum

Y <- c(5.3752, 5.6375, 5.5026, 5.9107, 5.5489, 5.6988, 5.782, 5.5134, 5.4937)
Ysum <- sum(Y)
PY <- Y/Ysum

HX <- - sum(PX* log2(PX))
HY <- - sum(PY* log2(PY))


MTX <- as.matrix(PX[1]*PY)

for(i in c(2:9)){
  Pxiyi = as.matrix(PX[i]*PY)
  MTX = cbind(MTX, Pxiyi)
}

colnames(MTX) <- c(1:9)

H = function(x){-(x*log2(x))}

HXYmatrix <- apply(MTX, 2, H)
HXY <- sum(HXYmatrix)

plot(c(0:8), X, type = "b", xlab = "Hour of Growth", ylab = "Expression Profile", main = "Profile Plot X")
plot(c(0:8), Y, type = "b", xlab = "Hour of Growth", ylab = "Expression Profile", main = "Profile Plot Y")
plot(c(0:8), PX, type = "b", xlab = "Hour of Growth", ylab = "Probability", main = "P_X")
plot(c(0:8), PY, type = "b", xlab = "Hour of Growth", ylab = "Probability", main = "P_Y")
