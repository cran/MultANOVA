#' Multiple ANOVAs (MultANOVA)
#'
#' Performs the MultANOVA omnibus test of each term from the model defined by \strong{formula} as defined in Mahieu & Cariou (2025).
#'
#' @param formula A formula with no left term that specify the model from the elements of the \strong{design} argument.
#' @param design A data.frame that contains only factors specifying the design on which rely the specified model of \strong{formula} argument.
#' @param responses A matrix or data.frame that contains only numerics or integers being the responses variables to be explained by the model from \strong{formula}.
#' @param graph A logical indicating whether or not multivariate type III r-squared by term must be plotted.
#' @param size.graph If \strong{graph=TRUE}, the overall size of labels and titles.
#'
#' @return An ANOVA-like table where each row corresponds to a term in the model or the residuals, the first column indicates the degrees of freedom, the second column corresponds to the multivariate type III r-squared and the third column corresponds to the p-value of the MultANOVA test.
#'
#'
#' @import stats
#'
#' @references Mahieu, B., & Cariou, V. (2025). MultANOVA Followed by Post Hoc Analyses for Designed High‐Dimensional Data: A Comprehensive Framework That Outperforms ASCA, rMANOVA, and VASCA. Journal of Chemometrics, 39(7). \doi{https://doi.org/10.1002/cem.70039}
#'
#' @export
#'
#' @examples
#'data(OTU)
#'an=MultANOVA(~Lot+Lactate+Atm+Time+Lactate:Atm+Lactate:Time,OTU[,1:4],OTU[,-c(1:4)])
#'print(an)


MultANOVA=function(formula,design,responses,graph=TRUE,size.graph=2.25){
  tested.terms="all"
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
  if (!is.logical(graph)){
    stop("class(graph) must be logical")
  }
  if (is.numeric(size.graph) | is.integer(size.graph)){
    if (size.graph<=0){
      stop("size.graph must be strictly positive")
    }
  }else{
    stop("class(size.graph) must be numeric or integer")
  }
  old.contr = options()$contrasts
  on.exit(options(contrasts = old.contr))
  options(contrasts = c("contr.sum","contr.sum"))
  effect.names=attr(terms(formula),"term.labels")
  if (any(!tested.terms%in%c("all",effect.names))){
    stop("tested.terms must be all or terms in the formula")
  }
  if (tested.terms=="all"){
    tested.terms=effect.names
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
  responses=as.matrix(responses)
  D=model.matrix(formula,design)
  myAnovaIII=function(vec){
    y=as.matrix(vec)
    b=solve(crossprod(D))%*%crossprod(D,y)
    E=y-D%*%b ; vE=length(y)-ncol(D)
    myan=rep(0,length(tested.terms)) ; names(myan)=tested.terms
    for (fa in tested.terms){
      ou.fact=which(effect.names==fa)
      design.fact=D[,attr(D,"assign")==ou.fact,drop=FALSE] ; vH=ncol(design.fact)
      design.orth=D[,attr(D,"assign")!=ou.fact,drop=FALSE]
      coef.fact=b[attr(D,"assign")==ou.fact,,drop=FALSE]
      coef.orth=b[attr(D,"assign")!=ou.fact,,drop=FALSE]
      effet=design.fact%*%coef.fact ; effet.O=effet-design.orth%*%solve(crossprod(design.orth))%*%crossprod(design.orth,effet)
      CMH=as.numeric(crossprod(effet.O)/vH) ; CMR=as.numeric(crossprod(E)/vE)
      myan[fa]=pf(CMH/CMR,vH,vE,lower.tail = FALSE)
    }
    return(myan)
  }
  if (length(tested.terms)>1){
    multan=apply(responses, 2, myAnovaIII) ; multan=t(apply(multan,1,p.adjust,method="fdr"))
  }else{
    multan=apply(responses, 2, myAnovaIII) ; multan=p.adjust(multan,method = "fdr")
  }
  retour.an=matrix(NA,length(tested.terms)+1,3)
  rownames(retour.an)=c(tested.terms,"Residuals") ; colnames(retour.an) = c("DF","R2","p.value")
  SST=sum(diag(cov(responses)*(nrow(responses)-1))) ; mod=lm(as.formula(paste("responses",as.character(formula)[2],sep="~")),design)
  retour.an[nrow(retour.an),2]=sum(diag(crossprod(residuals(mod))))/SST ; retour.an[nrow(retour.an),1]=mod$df.residual
  for (fa in tested.terms){
    ou.fact=which(effect.names==fa)
    design.fact=model.matrix(mod)[,attr(model.matrix(mod),"assign")==ou.fact,drop=FALSE] ; vH=ncol(design.fact)
    design.orth=model.matrix(mod)[,attr(model.matrix(mod),"assign")!=ou.fact,drop=FALSE]
    coef.fact=coef(mod)[attr(model.matrix(mod),"assign")==ou.fact,,drop=FALSE]
    coef.orth=coef(mod)[attr(model.matrix(mod),"assign")!=ou.fact,,drop=FALSE]
    effet=design.fact%*%coef.fact ; effet.O=effet-design.orth%*%solve(crossprod(design.orth))%*%crossprod(design.orth,effet)
    H=crossprod(effet.O)
    retour.an[fa,2]=sum(diag(H))/SST ; retour.an[fa,1]=vH
    if (length(tested.terms)>1){
      retour.an[fa,3]=min(multan[fa,],na.rm=TRUE)
    }else{
      retour.an[fa,3]=min(multan,na.rm=TRUE)
    }
  }
  retour.an=as.data.frame(retour.an)
  etoiles = rep("   ", length(retour.an$p.value))
  etoiles[retour.an$p.value <= 0.10] = ".  "
  etoiles[retour.an$p.value <= 0.05] = "*  "
  etoiles[retour.an$p.value <= 0.01] = "** "
  etoiles[retour.an$p.value <= 0.001] = "***"
  retour.an$Sig=etoiles
  names(retour.an)[4]="" ; retour.an["Residuals",4]=""
  retour.an$R2=round(retour.an$R2,3) ; retour.an$p.value=format.pval(retour.an$p.value,na.form = "")
  if (graph){
    df.plot=data.frame(Effect=rownames(retour.an),R2=retour.an[,2])
    df.plot$Effect=factor(df.plot$Effect,levels = rownames(retour.an))
    p=ggplot(df.plot,mapping = aes(x=df.plot[,1],y=df.plot[,2]))+theme_bw()
    p=p+xlab("Effect")+ylab("R2")+ggtitle("Variance Explained")
    p=p+theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 5*size.graph, face = "bold"),plot.title = element_text(hjust = 0.5, face = "bold",size = 7*size.graph),axis.text.y = element_text(size=3*size.graph),axis.text.x = element_text(face="bold",size=5*size.graph,angle = 30,hjust = 1))
    p=p+geom_col(color="black",fill="deepskyblue3",linewidth=1)
    p=p+ylim(0,1)
    print(p)
  }
  class(retour.an)=c("MultANOVA","data.frame")
  return(retour.an)
}
