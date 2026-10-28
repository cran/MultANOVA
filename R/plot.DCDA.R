#' Plots for \code{\link[MultANOVA]{DCDA}} objects
#'
#' Plots results from \code{\link[MultANOVA]{DCDA}} with possible customizations from the different arguments.
#'
#' @param x An object resulting from \code{\link[MultANOVA]{DCDA}}.
#' @param axes Which dimensions should be plotted?
#' @param pair.comp An optional objects resulting from \code{\link[MultANOVA]{MultLSD}}. When provided, no discriminated modalities of the factor at the \strong{alpha} risk will linked on the plot.
#' @param alpha The alpha risk to determine evaluate significance of the p-values from \strong{ pair.comp}.
#' @param select.var The indices of the variables to be plotted. By default all variables are plotted.
#' @param select.level The indices of the levels of the term to be plotted. By default all levels of the term are plotted.
#' @param title An optional character specifying a title for the plot.
#' @param size The overall size of labels, points, etc.
#' @param expansion.var The factor of expansion applied to variables coordinates to increase readability
#' @param ... further arguments passed to or from other methods.
#'
#' @import ggplot2
#' @import ggrepel
#' @import stats
#'
#' @return The required plot.
#'
#' @export
#'
#' @examples
#' data(OTU)
#' acd=DCDA(~Lot+Atm+Time,OTU[,1:4],OTU[,-c(1:4)],"Time")
#' lsd=MultLSD(~Lot+Atm+Time,OTU[,1:4],OTU[,-c(1:4)],"Time")
#' fish=FisherS(~Lot+Atm+Time,OTU[,1:4],OTU[,-c(1:4)],"Time")
#' plot(acd,axes = c(1,2),pair.comp = lsd,expansion.var = 1.5,select.var = which(fish[2,]<=0.05))

plot.DCDA=function(x,axes=c(1,2),pair.comp=NULL,alpha=0.05,select.var=1:nrow(x$var.coord),select.level=1:nrow(x$level.coord),title=NULL,size=2.25,expansion.var=1.25,...){
  if (is.numeric(axes) | is.integer(axes)){
    if (length(axes)==2){
      if (min(axes)<1){
        stop("min(axes) much be larger than or equal to 1")
      }
      if (max(axes)>ncol(x$level.coord)){
        stop("max(axes) much be lower than or equal to ncol(x$level.coord)")
      }
    }else{
      stop("length(axes) must equal 2")
    }
  }else{
    stop("class(axes) must be integer or numeric")
  }
  if (!inherits(pair.comp,"MultLSD") & !is.null(pair.comp)){
    stop("pair.comp must be NULL or a MultLSD object")
  }else{
    if (!is.null(pair.comp)){
      if (any(rownames(pair.comp)!=rownames(x$level.coord))){
        stop("levels of term do not match between x and pair.comp")
      }
    }
  }
  if (is.numeric(select.var) | is.integer(select.var) | is.null(select.var)){
    if (!is.null(select.var)){
      if (min(select.var)<1){
        stop("min(select.var) much be larger than or equal to 1")
      }
      if (max(select.var)>nrow(x$var.coord)){
        stop("max(select.var) much be lower than or equal to nrow(x$var.coord)")
      }
    }
  }else{
    stop("class(select.var) must be integer, numeric or NULL")
  }
  if (is.numeric(select.level) | is.integer(select.level) | is.null(select.level)){
    if (!is.null(select.level)){
      if (min(select.level)<1){
        stop("min(select.level) much be larger than or equal to 1")
      }
      if (max(select.level)>nrow(x$level.coord)){
        stop("max(select.level) much be lower than or equal to nrow(x$level.coord)")
      }
    }
  }else{
    stop("class(select.level) must be integer, numeric or NULL")
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
  if (!is.null(title) & !is.character(title)){
    stop("title must be NULL or of class character")
  }
  if (is.numeric(size)){
    if (length(size)>1){
      stop("length(size) must be 1")
    }
  }else{
    stop("size must be of class numeric")
  }
  if (is.numeric(expansion.var)){
    if (length(expansion.var)>1){
      stop("length(expansion.var) must be 1")
    }
  }else{
    stop("expansion.var must be of class numeric")
  }
  if (!is.null(pair.comp)){
    pval=pair.comp[upper.tri(pair.comp)]
    hot=x$error$projected.T2[upper.tri(x$error$projected.T2)]
    optimd=function(d){
      sortie=sum((pval-pf((x$error$error.df-d+1)/(x$error$error.df*d)*hot,d,x$error$error.df-d+1,lower.tail = FALSE))^2)
      return(sortie)
    }
    dd=optimize(optimd,interval = c(nrow(x$eigen),min(x$error$error.df+1,nrow(x$var.coord))))$minimum
    if (nrow(x$level.coord)>2){
      base.ell=t(x$raw.coef)%*%x$error$sigma.hat%*%x$raw.coef ; dilat.ell=sqrt(((x$error$error.df*dd)/(x$error$error.df-dd+1))*qf(1-alpha,dd,x$error$error.df-dd+1))
      ell=NULL
      for (l in rownames(x$level.coord)){
        cov.l=(base.ell*x$error$leverage.levels[l])[c(axes[1],axes[2]),c(axes[1],axes[2])]
        ell.l=ellipse::ellipse(cov.l,centre=x$level.coord[l,c(axes[1],axes[2])],t=dilat.ell)
        ell=rbind(ell,cbind(rep(l,nrow(ell.l)),ell.l))
      }
      ell=as.data.frame(ell) ; ell[,1] = as.factor(ell[,1]) ; ell[,2] = as.numeric(ell[,2]) ; ell[,3]=as.numeric(ell[,3])
    }else{
      poid=cbind(x$raw.coef,x$var.coord[,2,drop=FALSE])
      base.ell=t(poid)%*%x$error$sigma.hat%*%poid ; dilat.ell=sqrt(((x$error$error.df*dd)/(x$error$error.df-dd+1))*qf(1-alpha,dd,x$error$error.df-dd+1))
      ell=NULL
      for (l in rownames(x$level.coord)){
        cov.l=(base.ell*x$error$leverage.levels[l])[c(axes[1],axes[2]),c(axes[1],axes[2])]
        ell.l=ellipse::ellipse(cov.l,centre=x$level.coord[l,c(axes[1],axes[2])],t=dilat.ell)
        ell=rbind(ell,cbind(rep(l,nrow(ell.l)),ell.l))
      }
      ell=as.data.frame(ell) ; ell[,1] = as.factor(ell[,1]) ; ell[,2] = as.numeric(ell[,2]) ; ell[,3]=as.numeric(ell[,3])
    }
  }
  adjusted.var.coord=x$var.coord
  max.norm.level = max(abs(x$level.coord[,axes]))
  max.norm.var = max(abs(adjusted.var.coord[,axes]))
  adjusted.var.coord=adjusted.var.coord*max.norm.level/max.norm.var*expansion.var

  if (!is.null(pair.comp)){
    xmin=min(x$level.coord[,axes[1]],adjusted.var.coord[,axes[1]],ell[,2])
    xmax=max(x$level.coord[,axes[1]],adjusted.var.coord[,axes[1]],ell[,2])
    ymin=min(x$level.coord[,axes[2]],adjusted.var.coord[,axes[2]],ell[,3])
    ymax=max(x$level.coord[,axes[2]],adjusted.var.coord[,axes[2]],ell[,3])
    pmin=min(xmin,ymin)*1.05
    pmax=max(xmax,ymax)*1.05
  }else{
    xmin=min(x$level.coord[,axes[1]],adjusted.var.coord[,axes[1]])
    xmax=max(x$level.coord[,axes[1]],adjusted.var.coord[,axes[1]])
    ymin=min(x$level.coord[,axes[2]],adjusted.var.coord[,axes[2]])
    ymax=max(x$level.coord[,axes[2]],adjusted.var.coord[,axes[2]])
    pmin=min(xmin,ymin)*1.05
    pmax=max(xmax,ymax)*1.05
  }

  if (!is.null(select.var)){
    adjusted.var.coord=adjusted.var.coord[select.var,,drop=FALSE]
  }

  p=ggplot(as.data.frame(x$level.coord),aes(x=x$level.coord[,axes[1]],y=x$level.coord[,axes[2]]))+theme_bw()
  if (nrow(x$level.coord)>2){
    p=p+xlim(pmin,pmax)+ylim(pmin,pmax)+xlab(paste("Dim.",axes[1]," (",round(x$eigen[axes[1],2],2),"%)",sep=""))+ylab(paste("Dim.",axes[2]," (",round(x$eigen[axes[2],2],2),"%)",sep=""))+ggtitle(title)
  }else{
    p=p+xlim(pmin,pmax)+ylim(pmin,pmax)+xlab(paste("Dim.",axes[1]," (",round(x$eigen[axes[1],2],2),"%)",sep=""))+ylab("PC1 of Residuals")+ggtitle(title)
  }
  p=p+theme(axis.title.x = element_text(size = 5*size, face = "bold"),axis.title.y = element_text(size = 5*size, face = "bold"),plot.title = element_text(hjust = 0.5, face = "bold",size = 7*size),axis.text = element_text(size=3*size))
  p=p+geom_hline(yintercept=0,linetype="dashed",linewidth=1)+geom_vline(xintercept=0,linetype="dashed",linewidth=1)

  if (!is.null(pair.comp)){
    ell=ell[ell$V1%in%rownames(x$level.coord)[select.level],]
    ell=as.data.frame(ell) ; ell[,1] = as.factor(as.character(ell[,1])) ; ell[,2] = as.numeric(ell[,2]) ; ell[,3]=as.numeric(ell[,3])
    p=p+geom_path(data=as.data.frame(ell),aes(x=ell[,2],y=ell[,3],group=ell[,1]),colour="blue",linewidth=1)
    df.segment=NULL
    for (i in 1:(nrow(pair.comp)-1)){
      for (j in (i+1):ncol(pair.comp)){
        l.1=rownames(pair.comp)[i]
        l.2=colnames(pair.comp)[j]
        if (pair.comp[l.1,l.2]>alpha & i%in%select.level & j%in%select.level){
          l.1.coord=x$level.coord[l.1,axes,drop=FALSE]
          l.2.coord=x$level.coord[l.2,axes,drop=FALSE]
          sous.df.segment=cbind(l.1.coord,l.2.coord)
          df.segment=rbind(df.segment,sous.df.segment)
        }
      }
    }
    if (!is.null(df.segment)){
      colnames(df.segment)=as.character(1:ncol(df.segment))
      p=p+geom_segment(data=as.data.frame(df.segment),aes(x = df.segment[,1], y = df.segment[,2], xend = df.segment[,3], yend = df.segment[,4]),colour="blue",linewidth=1.3)
    }
  }

  if (!is.null(select.var)){
    df.fleche=as.matrix(adjusted.var.coord[,axes,drop=FALSE])
    col.fleche=rep("red",length(select.var))
    p=p+geom_segment(data = as.data.frame(df.fleche), aes(x=0, y=0,xend = df.fleche[,1], yend = df.fleche[,2]), arrow=arrow(length = unit(0.14*size, "cm"),type = "closed"), colour=col.fleche,linewidth=1)
  }

  if (!is.null(select.level)){
    df.point=x$level.coord[select.level,axes,drop=FALSE]
    col.point=rep("blue",length(select.level))
    p=p+geom_point(data=as.data.frame(df.point),aes(x=df.point[,1],y=df.point[,2]),colour=col.point,size=size)
  }

  lab=NULL
  col.lab=NULL
  if (!is.null(select.level)){
    lab=rbind(lab,df.point)
    col.lab=c(col.lab,rep("blue",nrow(df.point)))
  }
  if (!is.null(select.var)){
    lab=rbind(lab,df.fleche)
    col.lab=c(col.lab,rep("red",nrow(df.fleche)))
  }
  nudge=lab*0.01
  p=p+geom_label_repel(as.data.frame(lab),mapping=aes(x=lab[,1],y=lab[,2],label=rownames(lab)),linewidth=0,colour=col.lab,size=size*2,segment.size=1,label.padding = 0,
                       nudge_x = nudge[,1],nudge_y = nudge[,2],min.segment.length = 1)

  return(p)
}
