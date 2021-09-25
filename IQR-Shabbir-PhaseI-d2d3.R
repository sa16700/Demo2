library(mvtnorm);library(mnormt)
n = 5; m = 30; muy = 5; mux = 5; sy0 = sx0 = 1; sim = 1e5
mean2 = c(muy, mux); 

ryx_v = c(0.2,0.5,0.8,0.95) ; 
nryx_v = length(ryx_v)
d2mat = d3mat = matrix(, nrow=nryx_v, ncol = 9)

for(k in 1:nryx_v){
	ryx = ryx_v[k]
sigma2 = matrix(c(sy0^2,ryx*sy0*sx0,ryx*sy0*sx0,sx0^2),ncol=2)

Q1x = qnorm(0.25, 0, 1); Q3x = qnorm(0.75, 0, 1)

iqru = iqrr = iqrp = iqrapr = iqrg = iqrre = iqrpe = iqrrpow = iqrppow = matrix(,ncol=sim,nrow=m)
set.seed(1098)
for(i in 1:sim){
      for(j in 1:m){
	xn = rmvnorm(n = n, mean = mean2, sigma = sigma2)
	y = xn[,1]; x = xn[,2]
	Q1x = qnorm(0.25, mux, sx0); Q3x = qnorm(0.75,  mux, sx0)
	Q1y = qnorm(0.25, muy, sy0); Q3y = qnorm(0.75,  muy, sy0)

	q1y = quantile(y, 0.25, type = 6, names = F)
	q3y = quantile(y, 0.75, type = 6, names = F)
	q1x = quantile(x, 0.25, type = 6, names = F)
	q3x = quantile(x, 0.75, type = 6, names = F)
	P11xy1=pmnorm(c(Q1y,Q1x), mean2, sigma2, log = FALSE)
	fixy1=(P11xy1-0.25^2)/(0.25*(1-0.25))
	P11xy3=pmnorm(c(Q3y,Q3x), mean2, sigma2, log = FALSE)
	fixy3=(P11xy3-0.75^2)/(0.75*(1-0.75))

	fQy1=dnorm(Q1y, muy, sy0) ;fQx1=dnorm(Q1x, mux, sx0)
	fQy3=dnorm(Q3y, muy, sy0) ;fQx3=dnorm(Q3x, mux, sx0)

	byx1 = fixy1*fQx1/fQy1; byx3 = fixy3*fQx3/fQy3

	pow3 = Q3x*byx3/Q3y; pow1 = Q1x*byx1/Q1y
	pow3 = ryx^2; pow1 = ryx^2
	iqru[j,i] = q3y - q1y
	iqrr[j,i] = q3y * Q3x/q3x - q1y*Q1x/q1x
	iqrp[j,i] = q3y * q3x/Q3x - q1y*q1x/Q1x
	iqrapr[j,i] = (iqrr[j,i]+iqrp[j,i])/2
	iqrg[j,i] = q3y + byx3*(Q3x - q3x) - (q1y + byx1*(Q1x - q1x))
	iqrre[j,i] = q3y*exp((Q3x-q3x)/(Q3x+q3x)) - q1y*exp((Q1x-q1x)/(Q1x+q1x))
	iqrpe[j,i] = q3y*exp((q3x-Q3x)/(q3x+Q3x)) - q1y*exp((q1x-Q1x)/(q1x+Q1x))
	iqrrpow[j,i] = q3y * (Q3x/q3x)^pow3 - q1y*(Q1x/q1x)^pow1
	iqrppow[j,i] = q3y * (q3x/Q3x)^pow3 - q1y*(q1x/Q1x)^pow1

}
}
t1 = iqru; t2 = iqrr; t3 = iqrp; t4 = iqrapr; t5 = iqrg; t6 = iqrre; t7 = iqrpe
t8 = iqrrpow; t9 = iqrppow

d2e1 = mean(t1); d2e2 = mean(t2); d2e3 = mean(t3); d2e4 = mean(t4); 
d2e5 = mean(t5); d2e6 = mean(t6); d2e7 = mean(t7); d2e8 = mean(t8); d2e9 = mean(t9);
d2s0 = cbind(d2e1,d2e2,d2e3,d2e4,d2e5,d2e6,d2e7,d2e8,d2e9)

d3e1 = sd(t1); d3e2 = sd(t2);d3e3 = sd(t3); d3e4 = sd(t4);
d3e5 = sd(t5); d3e6 = sd(t6);d3e7 = sd(t7); d3e8 = sd(t8); d3e9 = sd(t9)
d3s0 = cbind(d3e1,d3e2,d3e3,d3e4,d3e5,d3e6,d3e7,d3e8,d3e9)

d2mat[k,] = d2s0; d3mat[k,] = d3s0
}
d2mat;d3mat
write.table(d2mat,"clipboard",sep="\t",col.names=F,row.names=F)
write.table(d3mat,"clipboard",sep="\t",col.names=F,row.names=F)

#library(beepr)
beep()

