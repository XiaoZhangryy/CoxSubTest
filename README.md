# CoxSubTest
Subgroup testing in Cox proportional hazards model with a change-plane.

Maximum likelihood ratio test is proposed for Cox proportional hazard models with a change plane. Different testing methods are provided in this package. 

# Installation

    #install Rtools 3.5 (http://cran.r-project.org/bin/windows/Rtools)
    #install.packages("devtools")
    #install.packages("Rcpp")
    library(devtools)
    install_github("XiaoZhangryy/CoxSubTest")

# Usage

- [x] [CoxSubTest-manual](https://github.com/XiaoZhangryy/inst/CoxSubTest-manual.pdf) ------------ Details of the usage of the package.

# Example

    library(CoxSubTest)
    library(mvtnorm)

    n = 100
    p1 = 2
    p2 = 1
    p3 = 3
    alpha = rep(1, p1)
    beta  = rep(1, p2)/2
    gamma = c(1, seq(-1,1,length.out = p3-1)) 
    rho = 0.3
    cenRate = 0.2
    data = generate_cox_data(n, alpha, beta, gamma, rho, cenRate = cenRate)
    fit <- CoxSubTestLRT(data)



# References

Subgroup testing in Cox proportional hazards model with a change-plane. Manuscript.

# Development
The R-package is developed by Xiao Zhang (zhangxiao1994@cuhk.edu.cn), Panpan Ren, Xu Liu and Xingjie Shi.




