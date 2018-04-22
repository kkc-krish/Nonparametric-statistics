### Chi-square independence test
library(vcd)
mytable=matrix(c(83,70,45,91,86,15,41,38,10),
               byrow=T,ncol = 3)
chisq.test(mytable)


### Spearman Kendall tao 秩相关检验
x=c(65,79,67,66,89,85,84,73,88,80,86,75)
y=c(62,66,50,68,88,86,64,62,92,64,81,80)
cor.test(x,y,method = "spearman")
cor.test(x,y,method = "kendall")


### Kendall 协和系数
friedman.test.fixed=function (y, groups, blocks, ...) 
{
    # 解决区组有相同的秩问题
    # 计算Kendall 协和系数
    DNAME <- deparse(substitute(y))
    if (is.matrix(y)) {
        groups <- factor(c(col(y)))
        blocks <- factor(c(row(y)))
    }else {
        if (anyNA(groups) || anyNA(blocks)) 
            stop("NA's are not allowed in 'groups' or 'blocks'")
        if (any(diff(c(length(y), length(groups), length(blocks))) != 
                0L)) 
            stop("'y', 'groups' and 'blocks' must have the same length")
        DNAME <- paste(DNAME, ", ", deparse(substitute(groups)), 
                       " and ", deparse(substitute(blocks)), sep = "")
        if (any(table(groups, blocks) != 1)) 
            stop("not an unreplicated complete block design")
        groups <- factor(groups)
        blocks <- factor(blocks)
        o <- order(groups, blocks)
        y <- y[o]
        groups <- groups[o]
        blocks <- blocks[o]
    }
    k <- nlevels(groups)
    y <- matrix(unlist(split(c(y), blocks)), ncol = k, byrow = TRUE)
    y <- y[complete.cases(y), ]
    n <- nrow(y)
    r <- t(apply(y, 1L, rank))
    TIES <- tapply(c(r), row(r), table)
    STATISTIC0 <- ((12 * sum((colSums(r) - n * (k + 1)/2)^2))/(n * k * (k + 1) 
                 -(sum(unlist(lapply(TIES, function(u)u^3-u)))/(k - 1))))
    knot=as.numeric(apply(y, 1, function(x){as.numeric(table(rank(x)[duplicated(rank(x))])+1)}))
    STATISTIC=STATISTIC0/(1-sum(knot^3-knot, na.rm = T)/(k*n*(n^2-1)))
    coef=STATISTIC/k/(n-1)
    PARAMETER <- k - 1
    PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)
    names(STATISTIC) <- "Friedman chi-squared"
    names(PARAMETER) <- "df"
    names(coef)="Kendall coefficient W"
    structure(list(statistic = c(STATISTIC, coef), parameter = PARAMETER, 
                   p.value = PVAL, method = "Friedman rank sum test", data.name = DNAME), 
              class = "htest")
}

mydata=matrix(c(1.5,2,2,1,1,1,1,1, 2,
                1.5,1,1,2,2,3,4,2, 2,
                3  ,3,3,3,4,4,3,5, 2,
                4  ,4,4,4,3,2,2,3, 4,
                5  ,5,5,5,5,5,5,4, 5),
              ncol=5)
friedman.test.fixed(mydata)
