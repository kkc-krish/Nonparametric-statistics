###Chi-squared Goodness-of-Fit Neyman-Pearson Test
### ch03 p48
chisq.test(c(315,101,108,32), p=c(9,3,3,1)/16)


###Kolmogorov-Smirnov test
### ch03 p65
healthy=c(0.4855,-0.005,-0.2762,1.2765,1.8634,-0.5226,
          0.1034,-0.8076,0.6804,-2.3646)
ks.test(healthy,pnorm,0,1)
