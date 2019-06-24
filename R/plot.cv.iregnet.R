##' Plot the variable weights and estimated train/validation loss as a
##' function of the model complexity.
##' @title Regularization path plot for cross-validated iregnet model
##' @param x object of class "cv.iregnet"
##' @param ... ignored
##' @export
##' @return object of class "ggplot"
##' @method plot cv.iregnet
##' @author Toby Dylan Hocking
##' @import ggplot2
plot.cv.iregnet <- function(x, ...){
  is.validation <- x$plot.data$stats$set.name == "validation"
  validation.stats <- x$plot.data$stats[is.validation, ]
  i.min <- x$selected[["min"]]
  min.stats <- validation.stats[i.min, ]
  min.upper <- with(min.stats, data.table(
    min.upper=mean+sd, facet="neg.loglik"))
  with(x$plot.data, {
    ggplot()+
      theme_bw()+
      theme(panel.spacing=grid::unit(0, "lines"))+
      facet_grid(facet ~ ., scales="free")+
      geom_line(aes(
        -log10(lambda),
        weight,
        group=variable),
        data=data.table(weights, facet="weight"))+
      geom_segment(aes(
        -log10(lambda),
        xend=-log10(lambda),
        y=mean-sd,
        yend=mean+sd,
        color=set.name),
        data=data.table(stats, facet="neg.loglik"))+
      geom_line(aes(
        -log10(lambda), neg.loglik,
        group=paste(set.name, validation.fold),
        color=set.name),
        data=data.table(likelihood, facet="neg.loglik"))+
      geom_point(aes(
        -log10(lambda),
        mean,
        color=set.name),
        data=data.table(validation.stats, facet="neg.loglik"))+
      ylab("-log10(likelihood)")+
      geom_point(aes(
        -log10(lambda),
        neg.loglik,
        color=set.name,
        fill=type),
        shape=21,
        data=data.table(selected, facet="neg.loglik"))+
      scale_fill_manual(values=c(min="black", "1sd"="white"))+
      geom_hline(aes(
        yintercept=min.upper),
        data=min.upper)+
      ylab("")+
      xlab("model complexity -log10(lambda)")
  })
}