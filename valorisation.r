# Rentabilite d'un portefeuille equipondere
func_1<-function(mu)
{
  if (is.null(ncol(mu)) != TRUE)
  stop("Le rendement est un vecteur")
  nb.actions <- length(mu)
  weights.vec <- rep(1,nb.actions)/nb.actions 
  mup <- crossprod(weights.vec,mu)
  cat("L'esperance de la rentabilité du portefeuille equipondere est : ", mup,"\n")
  mup
}


func_2 <- function(mu,Sigma,rentaCible){
  
  #Calcul des poids des actifs pour un niveau de rendement donne
  #Methode des cofacteurs --> R.C.Merton
  
  #Arguments d'entree :
  
  #mu : Considerons N actifs, on dispose d'un rendement de l'actif sur
  #la periode T, vecteur de N * 1 avec N le nombre de titres du portefeuille
  
  #Sigma : variance de dimension N*N, matrice symetrique, avec les variances
  #de chaque titre en diagonale
  
  #rentaCible : Vecteur 1*1 designant le pourcentage de rendement que
  #l'investisseur souhaite atteindre 
  
  #Arguments de sortie :
  
  #Determiner un vecteur de poids N*1 designant la part de chaque actif
  #dans le portefeuille en fonction de la rentaCible
  
  n <- rep(1, length(mu))
  
  if (is.null(ncol(mu)) != TRUE)
    stop("Le rendement doit etre un vecteur")
  
  if (isSymmetric.matrix(sigma) == FALSE)
    stop("La matrice de variance-covariance doit etre symetrique")
  
  if (length(mu) != dim(sigma)[1])
    stop("La dimension du vecteur de rentabilite est inapproprie")
  
  mat <- solve(Sigma) # inverse de la matrice
  A <- as.numeric(t(n) %*% mat %*% mu) 
  B <- as.numeric(t(mu) %*% mat %*% mu)
  C <- as.numeric(t(n) %*% mat %*% n)
  D <- B * C - A * A
  E <- mat %*% ( - A * n + C * mu) / D
  F <- mat %*% (B * n - A * mu) / D
  w <- E * rentaCible + F
  Risk <- sqrt(t(w) %*% Sigma %*% w) # calcul de l'ecart-type
  ExpReturn <- t(w) %*% mu
  Weights <- (round(w,5))
  list(Weights = t(Weights),
       ExpReturn = as.numeric(round(ExpReturn,5)),
       Risk = as.numeric(round(Risk,5)))
}


func_3 <- function(mu, Sigma, varianceCible){
  
  #Calcul des poids et de la rentabilite pour 
  #un niveau de risque donnee
  
  #Arguments d'entree :
  
  #mu : Considerons N actifs, on dispose d'un rendement de l'actif sur
  #la periode T, vecteur de N * 1 avec N le nombre de titres du portefeuille
  
  #Sigma : variance de dimension N*N, matrice symetrique, avec les variances
  #de chaque titre en diagonale
  
  #varianceCible : vecteur 1*1 designant le niveau de risque que l'investisseur 
  #est pret a assumer
  
  #Arguments de sortie :
  
  #Determiner un vecteur de poids N*1 designant la part de chaque actif
  #dans le portefeuille en fonction de la varianceCible
  
  n <- rep(1,length(mu))
  
  if (is.null(ncol(mu)) != TRUE)
    stop("Le rendement doit etre un vecteur")
  
  if (isSymmetric.matrix(sigma) == FALSE)
    stop ("La matrice de variance-covariance n'est pas symetrique")
  
  if (length(mu) != dim(sigma)[1])
    stop("La dimension du vecteur de rentabilite est inapproprie")
  
  if (varianceCible<0)
    stop("Le niveau de risque doit etre positif")
  
  mat <- solve(Sigma)
  A <- n %*% mat %*% mu
  B <- t(mu) %*% mat %*% mu
  C <- t(n) %*% mat %*% n
  D <- (B * C) - A^2
  F <- as.numeric(sqrt ((varianceCible*C-1)/D))
  G <- as.numeric((1/C) %*% (1 - A*F))
  W <- mat %*% (F * mu + G * n)
  Risk2.p <- t(W) %*% Sigma %*% W
  Risk <- sqrt(Risk2.p)
  Expreturn <- t(W) %*% mu
  Weights <- (round(W,5))
 list(Weights = t(Weights),
      Expreturn = as.numeric(round(Expreturn,5)),
      Risk = round(Risk,5))
}

#Portefeuille de variance minimum
func_4 <- function(mu, sigma){
  
  #Calcul les ponderations qui nous donne le niveau de risque le plus faible
  
  #Arguments d'entree :
  
  #mu : Considerons N actifs, on dispose d'un rendement de l'actif sur
  #la periode T, vecteur de N * 1 avec N le nombre de titres du portefeuille
  
  #Sigma : variance de dimension N*N, matrice symetrique, avec les variances
  #de chaque titre en diagonale
  
  #Arguments de sortie :
  #Poids : vecteur N*1
  
  n <- rep(1, length(mu)) # vecteur de 1
  
  if (is.null(ncol(mu)) != TRUE)
    stop("Le rendement doit etre un vecteur")
  
  if (isSymmetric.matrix(sigma) == FALSE)
    stop ("La matrice de variance-covariance n'est pas symetrique")
  
  if (length(mu) != dim(sigma)[1])
    stop("La dimension du vecteur de rentabilite est inapproprie")
  
  mat <- solve(sigma)     # inverse de la matrice de covariance
  w <- mat %*% n /as.numeric(t(n) %*% mat %*% n)
  Return <- t(w) * mu
  Weights <- (round(w, 5))
  list(Weights = t(Weights))
}


func_5 <- function(mu, sigma, riskFreeRate){
  
  #Calcul du portefeuille tangent 
  
  #Arguments d'entree :
  
  #mu : Considerons N actifs, on dispose d'un rendement de l'actif sur
  #la periode T, vecteur de N * 1 avec N le nombre de titres du portefeuille
  
  #Sigma : variance de dimension N*N, matrice symetrique, avec les variances
  #de chaque titre en diagonale
  
  #RiskFreeRate : vecteur 1*1 taux sans risque
  
  #Arguments de sortie :
  #Poids : vecteur N*1
  
  n <- rep(1, length(mu))
  
  if (is.null(ncol(mu)) != TRUE)
    stop("Le rendement doit etre un vecteur")
  
  if (isSymmetric.matrix(sigma) == FALSE)
    stop ("La matrice de variance-covariance n'est pas symetrique")
  
  if (length(mu) != dim(sigma)[1])
    stop("La dimension du vecteur de rentabilite est inapproprie")
  
  ExcessReturn <- (mu - riskFreeRate * n) # excedent de rendement
  mat <- solve(sigma) # inverse de la matrice
  w <- (mat %*% ExcessReturn)/ as.numeric(n%*%mat%*%ExcessReturn)
  list(weights = t(round(w,5)))
}

func_6 <- function(mu, sigma, rentaCible){
  
  #Calcul du portefeuille efficient en absence de vente a decouvert
  #Utilisation de la fonction solve.QP
  
  #Arguments d'entree :
  
  #mu : Considerons N actifs, on dispose d'un rendement de l'actif sur
  #la periode T, vecteur de N * 1 avec N le nombre de titres du portefeuille
  
  #Sigma : variance de dimension N*N, matrice symetrique, avec les variances
  #de chaque titre en diagonale
  
  #rentaCible : vecteur 1*1 de rendement espere par l'investisseur
  
  #Contrainte : le poids d'un actif ne doit pas etre inferieur a 0
                #la somme des poids doit etre egale a 1
  
  #Arguments de sortie :
  #Poids : vecteur N*1 des poids de chaque actif du portefeuille
  
  if (is.null(ncol(mu)) != TRUE)
    stop("Le rendement doit etre un vecteur")
  
  if (isSymmetric.matrix(sigma) == FALSE)
    stop ("La matrice de variance-covariance n'est pas symetrique")
  
  if (length(mu) != dim(sigma)[1])
    stop("La dimension du vecteur de rentabilite est inapproprie")
  
  D <- sigma
  target = rentaCible
  d <- rep(0,length(mu))
  meq <- 2
  A <- matrix(1, nrow=1, ncol= length(mu))
  AMAT <- t(rbind(mu,A,diag(length(mu))))
  b0 <- c(target,1,d)
  resultat <- solve.QP(Dmat=D, dvec=d, Amat= AMAT, bvec=b0, meq =2)
  w <- resultat$solution
  list(Weights = t(round(w,5)))
}


func_7 <- function(mu, sigma, RiskFreeRate, varianceCible){
  
  #Calcul du portefeuille efficient en presence d'un actif sans risque
  #et pour un niveau de risque donnee
  
  #Arguments d'entree :
  
  #mu : vecteur N*1 
  #sigma : vecteur N*N
  #RiskFreeRate vecteur 1*1 actif sans risque
  #varianceCible : vecteur 1*1
  
  #Arguments de sortie :
  
  #Poids : vecteur N * 1 
  
  if (is.null(ncol(mu)) != TRUE)
    stop("Le rendement doit etre un vecteur")
  
  if (isSymmetric.matrix(sigma) == FALSE)
    stop ("La matrice de variance-covariance n'est pas symetrique")
  
  if (length(mu) != dim(sigma)[1])
    stop("La dimension du vecteur de rentabilite est inapproprie")
  
  
  n <- rep(1,length(mu))
  ExcessReturn <- (mu - RiskFreeRate * n)
  mat <- solve(sigma)
  A <- n %*% mat %*% ExcessReturn
  B <- t(ExcessReturn) %*% mat %*% ExcessReturn
  C <- t(n) %*% mat %*% n
  D <- (B * C) - A^2
  F <- as.numeric(sqrt ((varianceCible*C-1)/D))
  G <- as.numeric((1/C) %*% (1 - A*F))
  W <- mat %*% (F * ExcessReturn + G * n)
  Risk2.p <- t(W) %*% sigma %*% W
  Risk <- sqrt(Risk2.p)
  Expreturn <- t(W) %*% ExcessReturn
  Weights <- (round(W,5))
  list(Weights = t(Weights),
       Expreturn = as.numeric(round(Expreturn,5)))
}



func_8 <- function(mu, sigma, delta){
  
  #portefeuille qui maximise w'R - (delta/2) * w'Ew avec w'e = 1 et delta=0.4
  
  if (is.null(ncol(mu)) != TRUE)
    stop("Le rendement doit etre un vecteur")
  
  if (isSymmetric.matrix(sigma) == FALSE)
    stop ("La matrice de variance-covariance n'est pas symetrique")
  
  if (length(mu) != dim(sigma)[1])
    stop("La dimension du vecteur de rentabilite est inapproprie")
  
  if (delta == 0)
    stop("Le delta doit etre different de 0")
  n <- rep(1, length(mu))
  mat <- solve(sigma) #Inverse de la matrice
  A <- as.numeric(t(mu) %*% mat %*% n)
  C <- as.numeric(t(n) %*% mat %*% n)
  W1 <- mat %*% (mu*C - (A-delta)*n)
  W2 <- as.numeric(delta*C)
  W <- W1/W2
  list(Weights = t(round(W,5)))
}