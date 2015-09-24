This package can be used to confirmatory test latent variable network models. See `?lvnet` for details. To install the package, use:

```r
library("devtools")
install_github("sachaepskamp/lvnet")
```

Requires OpenMx to be installed.

#### Example:
```r
# Load package:
library("lvnet")

# Load dataset:
library("lavaan")
data(HolzingerSwineford1939)
Data <- HolzingerSwineford1939[,7:15]

# Measurement model:
Lambda <- matrix(0, 9, 3)
Lambda[1:3,1] <- NA
Lambda[4:6,2] <- NA
Lambda[7:9,3] <- NA

# Fit CFA model:
CFA <- lvnet(Data, lambda = Lambda)

# Latent network:
Omega_psi <- matrix(c(
  0,NA,NA,
  NA,0,0,
  NA,0,0
),3,3,byrow=TRUE)

# Fit model:
LNM <- lvnet(Data, lambda = Lambda, omega_psi=Omega_psi)

# Compare fit:
lvnetCompare(cfa=CFA,lnm=LNM)
```
