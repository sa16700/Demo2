library(mvtnorm);library(mnormt)
data=read.table("D:\\2 - Backup Before QU\\Research + Backup\\Research\\New\\Shahid - CQD\\IQR-Aux-Phase I\\exp_data.txt")
y = data[,1]; x = data[,2]    # VapourPressure(Y)	InvTemp2 (Aux:X)
y1 = y1[1:250]; x1 = x1[1:250]

n = 5; m = 50; muy = mean(y); mux = mean(x); sy0 = sd(y); sx0 = sd(x); 
sim = 1e5
mean2 = c(muy, mux); 

ryx = cor(y,x)
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
t1 = iqru/sy0; t2 = iqrr; t3 = iqrp; t4 = iqrapr; t5 = iqrg/sy0; t6 = iqrre/sy0; t7 = iqrpe
t8 = iqrrpow; t9 = iqrppow

d2e1 = mean(t1); d2e2 = mean(t2); d2e3 = mean(t3); d2e4 = mean(t4); 
d2e5 = mean(t5); d2e6 = mean(t6); d2e7 = mean(t7); d2e8 = mean(t8); d2e9 = mean(t9);
d2s0 = cbind(d2e1,d2e2,d2e3,d2e4,d2e5,d2e6,d2e7,d2e8,d2e9)

d3e1 = sd(t1); d3e2 = sd(t2);d3e3 = sd(t3); d3e4 = sd(t4);
d3e5 = sd(t5); d3e6 = sd(t6);d3e7 = sd(t7); d3e8 = sd(t8); d3e9 = sd(t9)
d3s0 = cbind(d3e1,d3e2,d3e3,d3e4,d3e5,d3e6,d3e7,d3e8,d3e9)

mul = function(x) 
	x * shft1


m1v = 15#c(3,6,9,12)
shftv = 1#c(seq(1.25,3,by=0.25), 3.5,4,4.5,5,6,7)
nshft = length(shftv)
nm1v = length(m1v)
cnt = 0; 
prob1 = prob2 = prob3 = prob4 = double(nshft*nm1v)
prob5 = prob6 = prob7 = prob8 = prob9 = double(nshft*nm1v)
for(i1 in 1:nm1v){
	m1 = m1v[i1]
	m0 = m - m1
for(j1 in 1:nshft){
shft1 = c(rep(shftv[j1], m1), rep(1, m0))
cnt = cnt + 1

Lt1 = 4.273

t1s = apply(t1,2,mul)
pt1 = double(sim)
UCLt1 = Lt1 * apply(t1s,2,mean)/d2e1
for(k in 1:sim){
	pt1[k] = sum(t1s[,k] > UCLt1[k])/m
}
alps = mean(pt1); (alp1 = 1 - (1-alps)^m)
prob1[cnt] = alp1

Lt2 = 2.348;

t2s = apply(t2,2,mul)
pt2 = double(sim)
UCLt2 = Lt2*apply(t2s,2,mean)/d2e2
for(k in 1:sim){
	pt2[k] = sum(t2s[,k] > UCLt2[k])/m
}
alps = mean(pt2); (alp2 = 1 - (1-alps)^m)
prob2[cnt] = alp2

Lt3 = 5.302;

t3s = apply(t3,2,mul)
pt3 = double(sim)
UCLt3 = Lt3*apply(t3s,2,mean)/d2e3
for(k in 1:sim){
	pt3[k] = sum(t3s[,k] > UCLt3[k])/m
}
alps = mean(pt3); (alp3 = 1 - (1-alps)^m)
prob3[cnt] = alp3

Lt4 = 3.354;

t4s = apply(t4,2,mul)
pt4 = double(sim)
UCLt4 = Lt4*apply(t4s,2,mean)/d2e4 
for(k in 1:sim){
	pt4[k] = sum(t4s[,k] > UCLt4[k])/m
}
alps = mean(pt4); (alp4 = 1 - (1-alps)^m)
prob4[cnt] = alp4

Lt5 = 2.971
t5s = apply(t5,2,mul)
pt5 = double(sim)
UCLt5 = Lt5*apply(t5s,2,mean)/d2e5
for(k in 1:sim){
	pt5[k] = sum(t5s[,k] > UCLt5[k])/m
}
alps = mean(pt5); (alp5 = 1 - (1-alps)^m)

prob5[cnt] = alp5

Lt6 = 4.023;

t6s = apply(t6,2,mul)
pt6 = double(sim)
UCLt6 = Lt6 * apply(t6s,2,mean)/d2e6
for(k in 1:sim){
	pt6[k] = sum(t6s[,k] > UCLt6[k])/m
}
alps = mean(pt6); (alp6 = 1 - (1-alps)^m)
prob6[cnt] = alp6

Lt7 = 4.255;

t7s = apply(t7,2,mul)
pt7 = double(sim)
UCLt7 = Lt7 * apply(t7s,2,mean)/d2e7
for(k in 1:sim){
	pt7[k] = sum(t7s[,k] > UCLt7[k])/m
}
alps = mean(pt7); (alp7 = 1 - (1-alps)^m)
prob7[cnt] = alp7

Lt8 = 2.341;

t8s = apply(t8,2,mul)
pt8 = double(sim)
UCLt8 = Lt8 * apply(t8s,2,mean)/d2e8
for(k in 1:sim){
	pt8[k] = sum(t8s[,k] > UCLt8[k])/m
}
alps = mean(pt8); (alp8 = 1 - (1-alps)^m)
prob8[cnt] = alp8

Lt9 = 5.084;

t9s = apply(t9,2,mul)
pt9 = double(sim)
UCLt9 = Lt9 * apply(t9s,2,mean)/d2e9
for(k in 1:sim){
	pt9[k] = sum(t9s[,k] > UCLt9[k])/m
}
alps = mean(pt9); (alp9 = 1 - (1-alps)^m)
prob9[cnt] = alp9
}
}
L = c(Lt1,Lt2,Lt3,Lt4,Lt5,Lt6,Lt7,Lt8,Lt9)
rslt = cbind(prob1,prob2,prob3,prob4,prob5,prob6,prob7,prob8,prob9)
rslt
write.table(rslt,"clipboard",sep="\t",col.names=F,row.names=F)

write.table(rbind(L,rslt),"clipboard",sep="\t",col.names=F,row.names=F)
#write.table(rbind(L,d2s0,d3s0),"clipboard",sep="\t",col.names=F,row.names=F)

#write.table(rslt[,c(3,5)],"clipboard",sep="\t",col.names=F,row.names=F)



