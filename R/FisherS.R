#' Variable-wise Type III Fisher tests
#'
#' Computes the variable-wise Type III Fisher tests of a term from the model defined by \strong{formula}. This enables investigating the most discriminant variables for that term and possibly leads to variable selection.
#'
#' @param formula A formula with no left term that specify the model from the elements of the \strong{design} argument.
#' @param design A data.frame that contains only factors specifying the design on which rely the specified model of \strong{formula} argument.
#' @param responses A matrix or data.frame that contains only numerics or integers being the responses variables to be explained by the model from \strong{formula}.
#' @param term A character specifying the term from \strong{formula} for which the MultLSD tests must be performed.
#' @param alpha The alpha risk to evaluate significance of the p-values of each Type III Fisher test. Variables will be colored differently according to their significance with \strong{graph=TRUE}.
#' @param graph A logical indicating whether or not variable-wise Type III Fisher statistics must be plotted together with their significance and average value.
#' @param size.graph If \strong{graph=TRUE}, the overall size of labels and titles.
#'
#' @return A matrix with two rows and as much columns as responses. The first row contains the Type III Fisher statistic and the second row contains the fdr correct p-value.
#'
#'
#' @import stats
#'
#'
#'
#' @export
#'
#' @examples
#' data(OTU)
#' fish=FisherS(~Lot+Atm+Time,OTU[,1:4],OTU[,-c(1:4)],"Time")
#' print(fish)


FisherS=function(formula,design,responses,term,alpha=0.05,graph=TRUE,size.graph=2.25){
  if (!is.logical(graph)){
    stop("class(graph) must be logical" )
  }
  if (!inherits(formula,"formula")){
    stop("class(formula) must be formula")
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
  if (is.numeric(alpha) | is.integer(alpha)){
    if (length(alpha)==1){
      if (alpha>1 | alpha<0){
        stop("alpha must be between 0 and 1")
      }
    }else{
      stop("length(alpha) must equal 1")
    }
  }else{
    stop("class(alpha) must be integer or numeric")
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
  myFIII=function(vec){
    y=as.matrix(vec)
    b=solve(crossprod(D))%*%crossprod(D,y)
    E=y-D%*%b ; vE=length(y)-ncol(D)
    ou.fact=which(effect.names==term)
    design.fact=D[,attr(D,"assign")==ou.fact,drop=FALSE] ; vH=ncol(design.fact)
    design.orth=D[,attr(D,"assign")!=ou.fact,drop=FALSE]
    coef.fact=b[attr(D,"assign")==ou.fact,,drop=FALSE]
    coef.orth=b[attr(D,"assign")!=ou.fact,,drop=FALSE]
    effet=design.fact%*%coef.fact ; effet.O=effet-design.orth%*%solve(crossprod(design.orth))%*%crossprod(design.orth,effet)
    CMH=as.numeric(crossprod(effet.O)/vH) ; CMR=as.numeric(crossprod(E)/vE)
    return(c(CMH/CMR,pf(CMH/CMR,vH,vE,lower.tail = FALSE)))
  }
  fisherS=apply(responses, 2, myFIII)
  fisherS[2,]=p.adjust(fisherS[2,],method = "fdr")
  if (graph){
    df.p=as.data.frame(fisherS[1,]) ; colnames(df.p)="Fisher"
    df.p$var=rownames(df.p) ; df.p$var=factor(df.p$var,levels = colnames(fisherS))
    df.p$Sig=ifelse(fisherS[2,]<=alpha,paste("p<",alpha*100,"%",sep=""),paste("p>",alpha*100,"%",sep="")) ; df.p$Sig=as.factor(df.p$Sig)
    p=ggplot(data=df.p,mapping = aes(x=df.p[,2],y=df.p[,1],fill=df.p[,3],color=df.p[,3]))+theme_bw()
    p=p+xlab("Variable")+ylab("Fisher Statistic")+ggtitle(paste(term,"term"),subtitle = "Blue line is the average Fisher statistic")
    p=p+theme(axis.title.x = element_blank(),axis.title.y = element_text(size = 5*size.graph, face = "bold"),plot.title = element_text(hjust = 0.5, face = "bold",size = 7*size.graph),plot.subtitle = element_text(hjust = 0.5, color="blue",face = "bold",size = 5*size.graph),axis.text.y = element_text(size=3*size.graph),axis.text.x = element_text(face="bold",size=3*size.graph,angle = 30,hjust = 1),legend.title = element_text(size = 5*size.graph, face = "bold"))
    p=p+scale_fill_manual(name="Significance",values=c("green3","orangered"))+scale_color_manual(name="Significance",values=c("green3","orangered"))
    p=p+geom_col()
    p=p+geom_hline(yintercept=mean(fisherS[1,]),color="blue",linewidth=1)
    print(p)
  }
  rownames(fisherS)=c("Fisher","Type III p.value (fdr corrected)")
  return(fisherS)
}
