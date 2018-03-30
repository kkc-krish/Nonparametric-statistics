#Reference: Wang. J. L., Non-parametric statistical analysis,2006.4

###Brown-Mood Median Test
## page 55
BM.test=function(x,y,alternative)
{
    ##Accurate test, 
    # normal approximation, 
    # normal approximation after continuity correction
    ##coef alternative can be greater, less or both
    xy=c(x,y)
    md.xy=median(xy)
    t=sum(xy>md.xy)
    lx=length(x[x!=md.xy])
    ly=length(y[y!=md.xy])
    lxy=lx+ly
    A=sum(x>md.xy)  ##Test statistics A
    z1=(A-lx*t)/(lx+ly)/(lx*ly*t*(lx+ly-t)/(lx+ly)^3)^0.5  
    ##Normalized Statistics for Normal Approximation
    if(A>(min(lx,t)/2)){
        z2=(A+0.5-lx*t)/(lx+ly)/(lx*ly*t*(lx+ly-t)/(lx+ly)^3)^0.5
        ##Normalized statistics for normal approximation after continuity correction
    }else{z2=(A-0.5-lx*t)/(lx+ly)/(lx*ly*t*(lx+ly-t)/(lx+ly)^3)^0.5}
    if(alternative=="greater"){
        pv1=phyper(A,lx,ly,t) ##Accurate p-value
        pv2=pnorm(z1)         ##Normal approximation p-value
        pv3=pnorm(z2)         ##Normal approximation p-value after continuity correction
    }else if(alternative=="less"){
        pv1=1-phyper(A,lx,ly,t)
        pv2=1-pnorm(z1)
        pv3=1-pnorm(z2)
    }else if(alternative=="both"){
        pv1=2*min(phyper(A,lx,ly,t),1-phyper(A,lx,ly,t))
        pv2=2*min(pnorm(z1),1-pnorm(z1))
        pv3=2*min(pnorm(z2),1-pnorm(z2))
    }
    conting.table=matrix(c(A,lx-A,lx,t-A,ly-(t-A),ly,t,lxy-t,lxy),3,3)
    col.name=c("X","Y","X+Y")
    row.name=c(">MXY","<MXY","TOTAL")
    dimnames(conting.table)=list(row.name,col.name)
    list(contingency.table=conting.table,p.value=pv1,pvnorm=pv2,pvnr=pv3)
}

a=c(698,68,675,656,655,648,640,639,620)
b=c(780,754,740,712,693,680,621)
BM.test(a,b,alternative="greater")
BM.test(a,b,alternative="both")

women=c(28500,31000,22800,32350,30450,38200,34100,30150,33550,27350,
        25200,32050,26550,30650,35050,35600,26900,31350,28950,32900,
        31300,31350,35700,35900,35200,30450)
men=c(39700,33250,31800,38200,30800,32250,38050,34800,32750,38800,
      29900,37400,33700,36300,37250,33950,37750,36700,36100,26550,
      39200,41000,40400,35500)
BM.test(women,men,alternative="both")
BM.test(women,men,alternative="less")



###Mann-Whitney U Test & Wilcoxon Rank Sum Test
##  page 60 / page 75
## As wilcoxon rank sum statistic is equivalent to Man-Whitney statistic,
#  we only use one of them to test the hypothesis
## Also called Wilcoxon-Mann-Whitney U test
weight=c(70,83,85,94,97,101,104,107,112,113,118,119,123,124,129,132,134,146,161)
group=c(1,2,1,1,2,1,2,2,1,2,1,2,2,2,2,1,2,2,2)
by(weight, group, median)
wilcox.test(weight~group, exact=F,correct=F)

women=c(28500,31000,22800,32350,30450,38200,34100,30150,33550,27350,
        25200,32050,26550,30650,35050,35600,26900,31350,28950,32900,
        31300,31350,35700,35900,35200,30450)
men=c(39700,33250,31800,38200,30800,32250,38050,34800,32750,38800,
      29900,37400,33700,36300,37250,33950,37750,36700,36100,26550,
      39200,41000,40400,35500)
wilcox.test(women, men, exact = F, correct = F)
wilcox.test(women, men, exact = F, correct = F, alternative = "less")



###Ansari-Bradley Test
##  page 82
library(coin)
manmade=c(4.5,6.5,7,10,12)
machine=c(6,7.2,8,9,9.8)

ansari.test(manmade, machine, alternative="two.sided")


stock=c(1149,1152,1176,1149,1155,1169,1182,1160,1120,1171,
        1116,1130,1184,1194,1184,1147,1125,1125,1166,1151)
fmonth=factor(rep(0:1, c(10,10)), labels = c("nov", "dec"))
ansari_test(stock~fmonth, alternative="two.sided")


weight=c(4.95,5.45,5.50,5.75,5.48,5.26,5.33,5.20,5.12,
         5.16,5.20,5.22,5.38,5.90,5.36,5.25,5.28,5.35,
         5.20,5.22,4.98,5.45,5.70,5.34,5.18,
         5.22,5.30,5.34,5.28,5.29,5.25,5.30,5.27,5.16,
         5.38,5.34,5.35,5.19,5.35,5.05,5.36,5.28,5.33,
         5.30,5.28,5.30,5.20)
fmachine=factor(rep(1:2, c(25,22)))
by(weight, fmachine, sd)

ansari_test(weight~fmachine, alternative="two.sided")
ansari_test(weight~fmachine, alternative="g")



###Siegel-Tukey Test
##  page 83
## See [R-statistics blog](https://www.r-statistics.com/2010/02/siegel-tukey-a-
#  non-parametric-test-for-equality-in-variability-r-code/)
## Notice: cannot compute exact p-value with ties
source("https://raw.githubusercontent.com/talgalili/R-code-snippets/master/siegel.tukey.r")
weight=c(4.95,5.45,5.50,5.75,5.48,5.26,5.33,5.20,5.12,
         5.16,5.20,5.22,5.38,5.90,5.36,5.25,5.28,5.35,
         5.20,5.22,4.98,5.45,5.70,5.34,5.18,
         5.22,5.30,5.34,5.28,5.29,5.25,5.30,5.27,5.16,
         5.38,5.34,5.35,5.19,5.35,5.05,5.36,5.28,5.33,
         5.30,5.28,5.30,5.20)
machine=rep(0:1, c(25,22))
siegel.tukey(weight~machine)
