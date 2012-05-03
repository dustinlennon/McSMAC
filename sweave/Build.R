
library(cacheSweave)
setwd("..")
setCacheDir("sweave")
Sweave("./sweave/MCSMAC.Rnw",
       driver=cacheSweaveDriver
       )
