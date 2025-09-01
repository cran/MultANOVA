#' Multiple Least-Squares Difference tests (MultLSD)
#'
#' Performs the Multiple Least-Squares Difference tests of a term from the model defined by \strong{formula} as defined in Mahieu & Cariou (2025).
#'
#' @param formula A formula with no left term that specify the model from the elements of the \strong{design} argument.
#' @param design A data.frame that contains only factors specifying the design on which rely the specified model of \strong{formula} argument.
#' @param responses A matrix or data.frame that contains only numerics or integers being the responses variables to be explained by the model from \strong{formula}.
#' @param tested.term A character specifying the term from \strong{formula} for which the MultLSD tests must be performed.
#'
#' @return A symmetric matrix having as rows and columns as there is levels of the \strong{tested.term} and that contains the pairwise fdr corrected p-values of the MultLSD tests.
#'
#' @import stats
#'
#'
#' @references Mahieu, B., & Cariou, V. (2025). MultANOVA Followed by Post Hoc Analyses for Designed High‐Dimensional Data: A Comprehensive Framework That Outperforms ASCA, rMANOVA, and VASCA. Journal of Chemometrics, 39(7). \doi{https://doi.org/10.1002/cem.70039}
#'
#' @export
#'
#' @examples
#' data(OTU)
#' lsd=MultLSD(~Lot+Atm+Time,OTU[,1:4],OTU[,-c(1:4)],"Time")
#' print(lsd)

MultLSD=function(formula,design,responses,tested.term){
  permute <- function(vecteur) {
    if (length(vecteur) == 1) {
      toutes_les_perms=as.data.frame(vecteur) ; colnames(toutes_les_perms)=as.character(1:ncol(toutes_les_perms))
      return(toutes_les_perms)
    }
    toutes_les_perms <- list()
    for (i in seq_along(vecteur)) {
      element_courant <- vecteur[i]
      reste_du_vecteur <- vecteur[-i]
      permutations_du_reste <- permute(reste_du_vecteur)
      for (perm_reste in permutations_du_reste) {
        toutes_les_perms <- c(toutes_les_perms, list(c(element_courant, perm_reste)))
      }
    }
    toutes_les_perms=as.data.frame(toutes_les_perms) ; colnames(toutes_les_perms)=as.character(1:ncol(toutes_les_perms))
    return(toutes_les_perms)
  }
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
  if (is.character(tested.term)){
    if (length(tested.term)==1){
      if (!tested.term%in%c(effect.names)){
        stop("tested.term must be a term in the formula")
      }
    }else{
      stop("length(tested.term) must equal 1")
    }
  }else{
    stop("class(tested.term) must be character")
  }
  if (regexpr(":",tested.term)>0){
    tested.term.interac=TRUE
    compo.interac=strsplit(tested.term,"[:]")[[1]]
  }else{
    tested.term.interac=FALSE
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
  design.full=model.matrix(formula,design)
  compo.term=strsplit(tested.term,"[:]")[[1]]
  level.term=levels(interaction(design[,compo.term]))
  coupe.term=strsplit(tested.term,"[:]")[[1]] ; coupe.term.perm=permute(coupe.term) ; term.version=apply(coupe.term.perm,2,paste,collapse=":")
  dans.term=NULL
  for (st in effect.names){
    if (any(regexpr(st,term.version)>0)){
      dans.term=c(dans.term,st)
    }
  }
  design.term.corresp=rep(attr(terms(formula),"term.labels"),as.numeric(table(attr(design.full,"assign")[-1])))
  myLSD=function(j){
    y=as.matrix(responses[,j])
    b=solve(crossprod(design.full))%*%crossprod(design.full,y)
    E=y-design.full%*%b ; vE=length(y)-ncol(design.full)
    SE=crossprod(E)/vE
    myposthoc=matrix(1,length(level.term),length(level.term)) ; rownames(myposthoc)=colnames(myposthoc)=level.term
    for (i in 1:(length(level.term)-1)){
      for (j in (i+1):length(level.term)){
        lei=level.term[i] ; lej=level.term[j] ; oui=match(lei,interaction(design[,compo.term])) ; ouj = match(lej,interaction(design[,compo.term]))
        Ci = design.full[oui,] ; Ci[!c(TRUE,design.term.corresp%in%dans.term)]=0
        Cj = design.full[ouj,] ; Cj[!c(TRUE,design.term.corresp%in%dans.term)]=0
        Cdiff=Ci-Cj ; xbarre.diff=t(as.matrix(Cdiff))%*%b
        leverage.diff=as.numeric(t(as.matrix(Cdiff))%*%solve(crossprod(design.full))%*%as.matrix(Cdiff))
        StdE.diff=(leverage.diff*SE)
        t.stat=as.numeric(xbarre.diff)/sqrt(as.numeric(StdE.diff))
        sur=pt(t.stat,vE,lower.tail = FALSE) ; sous=pt(t.stat,vE,lower.tail = TRUE)
        myposthoc[lei,lej]=myposthoc[lej,lei]=min(c(sur,sous))*2
      }
    }
    return(myposthoc)
  }
  multLSD=lapply(1:ncol(responses),myLSD) ; tensLSD=array(unlist(multLSD),dim=c(length(level.term),length(level.term),ncol(responses)),dimnames = list(levels=level.term,levels=level.term,y=colnames(responses)))
  retour.pval=matrix(1,length(level.term),length(level.term)) ; rownames(retour.pval)=colnames(retour.pval)=level.term
  for (i in 1:(length(level.term)-1)){
    for (j in (i+1):length(level.term)){
      lei=level.term[i] ; lej=level.term[j]
      pvalij=min(p.adjust(tensLSD[lei,lej,],method = "fdr"),na.rm=TRUE)
      retour.pval[lei,lej]=retour.pval[lej,lei]=pvalij
    }
  }
  retour.pval[upper.tri(retour.pval)]=p.adjust(retour.pval[upper.tri(retour.pval)],method = "fdr")
  retour.pval[lower.tri(retour.pval)]=p.adjust(retour.pval[lower.tri(retour.pval)],method = "fdr")
  class(retour.pval)="MultLSD"
  return(retour.pval)
}
