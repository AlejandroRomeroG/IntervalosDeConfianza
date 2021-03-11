library(greekLetters)

# INTERVALOS DE CONFIANZA PARA MEDIAS #

ConfianzaMediasVarianzaConocidaMediaDada <- function(n,certeza,media,varianza){
  alfa <- 1-certeza
  alfachida <- alfa/2
  z <- qnorm(alfachida)
  
  vecindadderecha <- media-z*(sqrt(varianza)/sqrt(n))
  vecindadizquierda <- media+z*(sqrt(varianza)/sqrt(n))
  
  paste("P(", vecindadizquierda, " < ", greeks("mu") ," < ",vecindadderecha, " ) = ", certeza)
}

ConfianzaMediasVarianzaConocidaMediaDesconocida <- function(datos,certeza,varianza){
  n <- length(datos)
  media <- mean(datos)
  alfa <- 1-certeza
  alfachida <- alfa/2
  z <- qnorm(alfachida)
  
  vecindadderecha <- media-z*(sqrt(varianza)/sqrt(n))
  vecindadizquierda <- media+z*(sqrt(varianza)/sqrt(n))
  
  paste("P(", vecindadizquierda, " < ", greeks("mu") ," < ",vecindadderecha, " ) = ", certeza)
}

ConfianzaMediasVarianzaDesconocidaMediaDada <- function(n,certeza,media,desvestandar){
  alfa <- 1-certeza
  alfachida <- alfa/2
  t <- qt(alfachida,n-1)
  
  vecindadderecha <- media-t*(desvestandar/sqrt(n))
  vecindadizquierda <- media+t*(desvestandar/sqrt(n))
  
  paste("P(", vecindadizquierda, " < ", greeks("mu") ," < ",vecindadderecha, " ) = ", certeza)
}

ConfianzaMediasVarianzaDesconocidaMediaDesconocida <- function(datos,certeza,desvestandar){
  n <- length(datos)
  media <- mean(datos)
  alfa <- 1-certeza
  alfachida <- alfa/2
  t <- qt(alfachida,n-1)
  
  vecindadderecha <- media-t*(desvestandar/sqrt(n))
  vecindadizquierda <- media+t*(desvestandar/sqrt(n))
  
  paste("P(", vecindadizquierda, " < ", greeks("mu") ," < ",vecindadderecha, " ) = ", certeza)
}



# INTERVALOS DE CONFIANZA PARA DIFERENCIA DE MEDIAS #

ConfianzaDiferenciaDeMediasConVarianzas <- function(media1,media2,varianza1,varianza2,n1,n2,certeza){
  alfa <- 1-certeza
  alfachida <- alfa/2
  z <- qnorm(alfachida)
  
  vecindadderecha <- (media1-media2)-z*(sqrt((varianza1/n1)+(varianza2/n2)))
  vecindadizquierda <- (media1-media2)+z*(sqrt((varianza1/n1)+(varianza2/n2)))
  
  paste("P(", vecindadizquierda, " < ", "µ_1 - µ_2"  ," < ",vecindadderecha, " ) = ", certeza)
}

ConfianzaDiferenciaDeMediasConDesviacionesEstandar <- function(media1,media2,desvest1,desvest2,n1,n2,certeza){
  alfa <- 1-certeza
  alfachida <- alfa/2
  z <- qnorm(alfachida)
  
  vecindadderecha <- (media1-media2)-z*(sqrt((desvest1^2/n1)+(desvest2^2/n2)))
  vecindadizquierda <- (media1-media2)+z*(sqrt((desvest1^2/n1)+(desvest2^2/n2)))
  
  paste("P(", vecindadizquierda, " < ", "µ_1 - µ_2"  ," < ",vecindadderecha, " ) = ", certeza)
}

ConfianzaDiferenciaDeMediasConLasMismasVarianzas <- function(media1,media2,desvest1,desvest2,n1,n2,certeza){
  alfa <- 1-certeza
  alfachida <- alfa/2
  t <- qt(alfachida,n1+n2-2)
  sigmap <- sqrt((((n1-1)*desvest1^2)+((n2-1)*desvest2^2))/(n1+n2-2))
  
  vecindadderecha <- (media1-media2)-t*sigmap*sqrt(1/n1 + 1/n2)
  vecindadizquierda <- (media1-media2)+t*sigmap*sqrt(1/n1 + 1/n2)
  
  paste("P(", vecindadizquierda, " < ", "µ_1 - µ_2"  ," < ",vecindadderecha, " ) = ", certeza)
}


# INTERVALO DE CONFIANZA PARA PROPORCIONES #

ConfianzaProporciones <- function(k,n,certeza){
  alfa <- 1-certeza
  alfachida <- alfa/2
  z <- qnorm(alfachida)
  estimador <- k/n
  
  vecindadderecha <- estimador-z*sqrt((estimador*(1-estimador))/n)
  vecindadizquierda <- estimador+z*sqrt((estimador*(1-estimador))/n)
  
  paste("P(", vecindadizquierda, " < ", "θ" ," < ",vecindadderecha, " ) = ", certeza)
}


# INTERVALO DE CONFIANZA PARA DIFERENCIA DE PROPORCIONES #

ConfianzaDiferenciaDeProporciones <- function(k1,k2,n1,n2,certeza){
  alfa <- 1-certeza
  alfachida <- alfa/2
  z <- qnorm(alfachida)
  estimador1 <- k1/n1
  estimador2 <- k2/n2
  
  vecindadderecha <- (estimador1-estimador2)-z*(sqrt((estimador1*(1-estimador1))/n1+(estimador2*(1-estimador2))/n2))
  vecindadizquierda <- (estimador1-estimador2)+z*(sqrt((estimador1*(1-estimador1))/n1+(estimador2*(1-estimador2))/n2))
  
  paste("P(", vecindadizquierda, " < ", "θ_1 - θ_2"  ," < ",vecindadderecha, " ) = ", certeza)
}


# INTERVALO DE CONFIANZA PARA VARIANZAS #

ConfianzaVarianzasDesviacionEstandar <- function(n,desvest,certeza){
  alfa <- 1-certeza
  alfachida <- alfa/2
  Xizquierda <- qchisq(alfachida,n-1)
  Xderecha <- qchisq(1-alfachida,n-1)
  varianza <- desvest^2
  
  vecindadderecha <- ((n-1)*varianza)/Xizquierda
  vecindadizquierda <- ((n-1)*varianza)/Xderecha
  
  paste("P(", vecindadizquierda, " < ", "σ^2" ," < ",vecindadderecha, " ) = ", certeza)
}

ConfianzaVarianza <- function(datos,certeza){
  n <- length(datos)
  media <- mean(datos)
  varianza <- var(datos)
  alfa <- 1-certeza
  alfachida <- alfa/2
  Xizquierda <- qchisq(alfachida,n-1)
  Xderecha <- qchisq(1-alfachida,n-1)
  
  vecindadderecha <- ((n-1)*varianza)/Xizquierda
  vecindadizquierda <- ((n-1)*varianza)/Xderecha
  
  paste("P(", vecindadizquierda, " < ", "σ^2" ," < ",vecindadderecha, " ) = ", certeza)
}


# INTERVALO DE CONFIANZA PARA RAZON DE VARIANZAS #

ConfianzaRazonDeVarianzas <- function(desvest1,desvest2,n1,n2,certeza){
  alfa <- 1-certeza
  alfachida <- alfa/2
  varianza1 <- desvest1^2
  varianza2 <- desvest2^2
  fizquierda <- 1/qf(alfachida,n2-1,n1-1)
  fderecha <- qf(alfachida,n2-1,n1-1)
  
  vecindadderecha <- fizquierda*(varianza2/varianza1)
  vecindadizquierda <- fderecha*(varianza2/varianza1)
  
  paste("P(", vecindadizquierda, " < ", "(σ_2)^2 / (σ_1)^2" ," < ",vecindadderecha, " ) = ", certeza)
}