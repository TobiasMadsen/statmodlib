# statmodlib

Er en lille R pakke til faget Statistiske Modeller 2 på Aarhus Universitet. Den kan installeres med følgende kommandoer

```R
install.packages("devtools")
devtools::install_github('tobiasmadsen/statmodlib')
library(statmodlib)
```

Der er eksempler på løsninger til opgave 3.4, eksamen sommeren 2005 opgave 4 samt eksamen sommeren 2006 opgave 4. De kan ses ved

```R
vignette(topic = "opgave2_4", package = "statmodlib")
vignette(topic = "sommer05_opg4", package = "statmodlib")
vignette(topic = "sommer06_opg4", package = "statmodlib")
```