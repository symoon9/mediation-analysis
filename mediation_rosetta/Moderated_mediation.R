#############################
# Mediation Analysis
# rosetta
# DATA: GPS, ELS, BRAIN
#############################

library(rosetta)
data("cpbExample")

# Mediation 
result <- gemm(dat = cpbExample,
               xvar = "procJustice",
               mvars = c("cynicism","trust"),
               yvar = "CPB",
               nboot = 100)
print(result)

# Moderated Mediation
result2 <- gemm(dat = cpbExample,
               xvar = "procJustice",
               mvars = c("cynicism","trust"),
               yvar = "CPB",
               xmmod = "insecure",
               mymod = "gender" ,
               cmvars = c("age"),
               nboot = 500)
print(result2)
plotIMM(result2)
plotSS(result2)
