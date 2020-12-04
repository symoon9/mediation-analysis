## Mediation Analysis

### Moderated Mediation Analysis using library(rosetta)
```
result <- gemm(dat = cpbExample,
                xvar = "procJustice",
                mvars = c("cynicism","trust"),
                yvar = "CPB",
                xmmod = "insecure",
                mymod = "gender" ,
                cmvars = c("age"),
                nboot = 500)
```
- `xvar`: predictor
- `yvar`: dependent variable
- `mvars`: a vector of mediators
- `xmmod`: moderated value in path x->m
- `mymod`: moderated value in path m->y
- `cmvars`: covariates for mediators
- `cyvars`: covariates for dependent variable
- `nboot`: bootstrapping