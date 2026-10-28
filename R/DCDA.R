#' Diagonal Canonical Discriminant Analysis (DCDA)
#'
#' Performs the Diagonal Canonical Discriminant Analysis of a term from the model defined by \strong{formula} as defined in Mahieu & Cariou (2025).
#'
#' @param formula A formula with no left term that specify the model from the elements of the \strong{design} argument.
#' @param design A data.frame that contains only factors specifying the design on which rely the specified model of \strong{formula} argument.
#' @param responses A matrix or data.frame that contains only numerics or integers being the responses variables to be explained by the model from \strong{formula}.
#' @param term A character specifying the term from \strong{formula} for which the MultLSD tests must be performed.
#'
#' @return Returns a list of the following elements:\cr
#' @return \item{scores}{A data.frame giving the scores of the individuals of the diagonal canonical variates.}
#' @return \item{eigen}{A matrix giving the eigenvalues, corresponding percentages of inertia and cumulative percentages of inertia of the DCDA.}
#' @return \item{level.coord}{A matrix giving the coordinates of the levels of the \strong{term} on the diagonal canonical variates.}
#' @return \item{var.coord}{A matrix giving the coordinates of the variables of the on the diagonal canonical variates.}
#' @return \item{raw.coef}{A matrix giving the raw diagonal canonical coefficients.}
#' @return \item{error}{A list containing several information relative to the error. Most of this information are used to compute the approximate ellipses when calling \code{\link[MultANOVA]{plot.DCDA}}.}
#' @return \item{svd}{Results of the svd of the weighted orthogonalized effect matrix of \strong{term}.}
#'
#' @import stats
#'
#' @references Mahieu, B., & Cariou, V. (2025). MultANOVA Followed by Post Hoc Analyses for Designed High‐Dimensional Data: A Comprehensive Framework That Outperforms ASCA, rMANOVA, and VASCA. Journal of Chemometrics, 39(7). \doi{https://doi.org/10.1002/cem.70039}
#' @references Dudoit, S., Fridlyand, J., & Speed, T. P. (2002). Comparison of Discrimination Methods for the Classification of Tumors Using Gene Expression Data. Journal of the American Statistical Association, 97(457), 77–87. \doi{https://doi.org/10.1198/016214502753479248}
#'
#' @export
#'
#' @examples
#' data(OTU)
#' acd=DCDA(~Lot+Atm+Time,OTU[,1:4],OTU[,-c(1:4)],"Time")
#' lsd=MultLSD(~Lot+Atm+Time,OTU[,1:4],OTU[,-c(1:4)],"Time")
#' fish=FisherS(~Lot+Atm+Time,OTU[,1:4],OTU[,-c(1:4)],"Time")
#' plot(acd,axes = c(1,2),pair.comp = lsd,expansion.var = 1.5,select.var = which(fish[2,]<=0.05))

DCDA=function(formula,design,responses,term){
  if (!inherits(formula,"formula")){
    stop("class(formula) must be formula")
  }
  if (is.data.frame(design)){
    for (j in 1:ncol(design)){
      if (!class(design[,j])%in%c("factor")){
        stop("design must be composed of only factors")
      }
    }
  }else{
    stop("class(design) must be data.frame")
  }
  if (is.data.frame(responses) | is.matrix(responses)){
    for (j in 1:ncol(responses)){
      if (!class(responses[,j])%in%c("numeric","integer")){
        stop("responses must be composed of only numerics or integers")
      }
    }
  }else{
    stop("class(responses) must be data.frame or matrix")
  }
  vari=apply(responses, 2, sd)
  if (any(vari<=1e-12)){
    ou.vari.nulle=which(vari<=1e-12)
    stop(paste("response(s) number ",paste(ou.vari.nulle,collapse = ", ")," have too low variance",sep=""))
  }
  old.contr = options()$contrasts
  on.exit(options(contrasts = old.contr))
  options(contrasts = c("contr.sum","contr.sum"))
  effect.names=attr(terms(formula),"term.labels")
  if (is.character(term)){
    if (length(term)==1){
      if (!term%in%c(effect.names)){
        stop("term must be a term in the formula")
      }
    }else{
      stop("length(term) must equal 1")
    }
  }else{
    stop("class(term) must be character")
  }
  pres.interact=any(regexpr(":",effect.names)>0)
  if (pres.interact){
    vec.f=NULL
    for (f in effect.names){
      vec.f=c(vec.f,strsplit(f,":")[[1]])
    }
    fact.names=unique(vec.f)
  }else{
    fact.names=effect.names
  }
  compo.term=strsplit(term,"[:]")[[1]]
  level.term=levels(interaction(design[,compo.term]))
  dans.term=NULL
  for (st in effect.names){
    if (all(strsplit(st,":")[[1]] %in% compo.term)){
      dans.term=c(dans.term,st)
    }
  }
  ordre=order(interaction(design[,compo.term]))
  design=design[ordre,,drop=FALSE]
  responses=responses[ordre,,drop=FALSE]
  responses=as.matrix(responses)
  myformula=as.formula(paste("responses",as.character(formula)[2],sep="~"))
  mod=lm(myformula,design)
  design.full=model.matrix(mod)
  error=residuals(mod)
  E=crossprod(error) ; vE=mod$df.residual ; SE=E/vE
  dimen=min(length(level.term)-1,ncol(responses))
  design.term.corresp=rep(attr(terms(formula),"term.labels"),as.numeric(table(attr(design.full,"assign")[-1])))
  select.fact=design.term.corresp%in%dans.term ; select.orth=!design.term.corresp%in%dans.term
  coef.fact=coef(mod)[c(FALSE,select.fact),,drop=FALSE]
  design.fact=design.full[,c(FALSE,select.fact),drop=FALSE]
  design.orth=design.full[,c(TRUE,select.orth),drop=FALSE]
  effet=design.fact%*%coef.fact ; effetO=effet-design.orth%*%solve(crossprod(design.orth))%*%crossprod(design.orth,effet)
  std.effetO=effetO%*%diag(1/sqrt(diag(E)))
  udv=svd(std.effetO,nu=dimen,nv=dimen)
  udv$d=udv$d[1:dimen]
  eig=matrix(0,dimen,3) ; rownames(eig)=paste("Dim",1:dimen,sep=".") ; colnames(eig)=c("eigenvalue","%inertia","cum%inertia")
  eig[,1]=udv$d^2 ; eig[,2]=(udv$d^2)/sum(udv$d^2)*100 ; eig[,3]=cumsum(eig[,2])
  if (dimen>1){
    scores=(effet+error)%*%diag(1/sqrt(diag(E)))%*%udv$v
    scores=cbind.data.frame(interaction(design[,compo.term]),scores)
    colnames(scores)=c("level",rownames(eig))
    error.coord=(error)%*%diag(1/sqrt(diag(E)))%*%udv$v ; colnames(error.coord)=rownames(eig)
  }else{
    scores=(effet+error)%*%diag(1/sqrt(diag(E)))%*%udv$v
    udv.error=svd(E,nu=1,nv=1)
    scores=cbind.data.frame(interaction(design[,compo.term]),scores,error%*%udv.error$v[,1,drop=FALSE])
    colnames(scores)=c("level",rownames(eig),"PC1 of Residuals")
    error.coord=(error)%*%diag(1/sqrt(diag(E)))%*%udv$v
    error.coord=cbind(error.coord,error%*%udv.error$v[,1,drop=FALSE]) ; colnames(error.coord)=c(rownames(eig),"PC1 of Residuals")
  }
  if (dimen>1){
    var.coord=udv$v
    rownames(var.coord)=colnames(responses) ; colnames(var.coord)=rownames(eig)
    rc=diag(1/sqrt(diag(E)))%*%udv$v
    rownames(rc)=colnames(responses) ; colnames(rc)=rownames(eig)
  }else{
    var.coord=udv$v
    rownames(var.coord)=colnames(responses) ; colnames(var.coord)=rownames(eig)
    udv.error=svd(E,nu=1,nv=1)
    var.coord=cbind(var.coord,udv.error$v[,1,drop=FALSE]) ; colnames(var.coord)[ncol(var.coord)]="PC1 of Residuals"
    rc=diag(1/sqrt(diag(E)))%*%udv$v
    rownames(rc)=colnames(responses) ; colnames(rc)=rownames(eig)
  }
  level.coord=aggregate(.~level,scores,mean) ; rownames(level.coord)=level.coord$level ; level.coord$level=NULL
  vec.leverage=rep(0,nrow(level.coord)) ; names(vec.leverage)=rownames(level.coord)
  for (l in rownames(level.coord)){
    oul=match(l,interaction(design[,compo.term]))
    Cl = design.full[oul,] ; Cl[c(FALSE,select.orth)]=0
    leverage.l=as.numeric(t(as.matrix(Cl))%*%solve(crossprod(design.full))%*%as.matrix(Cl))
    vec.leverage[l]=leverage.l
  }
  vec.leverage=round(vec.leverage,12)
  coeff.full=coef(mod)
  mat.T2.stat=matrix(0,length(level.term),length(level.term)) ; rownames(mat.T2.stat)=colnames(mat.T2.stat)=level.term
  for (i in 1:(length(level.term)-1)){
    for (j in (i+1):length(level.term)){
      lei=level.term[i] ; lej=level.term[j] ; oui=match(lei,interaction(design[,compo.term])) ; ouj = match(lej,interaction(design[,compo.term]))
      Ci = design.full[oui,] ; Ci[!c(TRUE,design.term.corresp%in%dans.term)]=0
      Cj = design.full[ouj,] ; Cj[!c(TRUE,design.term.corresp%in%dans.term)]=0
      Cdiff=Ci-Cj ; xbarre.diff=t(as.matrix(Cdiff))%*%coeff.full%*%rc
      leverage.diff=as.numeric(t(as.matrix(Cdiff))%*%solve(crossprod(design.full))%*%as.matrix(Cdiff))
      StdE.diff=(leverage.diff*t(rc)%*%SE%*%rc)
      T2.stat=as.numeric(xbarre.diff%*%solve(StdE.diff)%*%t(xbarre.diff))
      mat.T2.stat[lei,lej]=mat.T2.stat[lej,lei]=T2.stat
    }
  }
  retour.error=list(error.coord=as.matrix(error.coord),sigma.hat=SE,error.df=vE,leverage.levels=vec.leverage,projected.T2=mat.T2.stat)
  retour=list(eigen=as.matrix(eig),level.coord=as.matrix(level.coord),var.coord=as.matrix(var.coord),raw.coef=rc,scores=scores,error=retour.error,svd=udv)
  class(retour)="DCDA"
  return(retour)
}
