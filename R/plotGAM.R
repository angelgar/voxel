#' GAM plotting using ggplot2
#'
#' @family Plotting
#' @param gamFit fitted gam model as produced by mgcv::gam()
#' @param smooth.cov (character) name of smooth term to be plotted
#' @param groupCovs (character)  name of group variable to plot by, if NULL (default) then there are no groups in plot
#' @param orderedAsFactor if TRUE then the model is refitted with ordered variables as factors.
#' @param rawOrFitted If FALSE (default) then only smooth terms are plotted; if rawOrFitted = "raw" then raw values are plotted against smooth; if rawOrFitted = "fitted" then fitted values are plotted against smooth
#' @param plotCI if TRUE (default) upper and lower confidence intervals are added at 2 standard errors above and below the mean
#'
#' @return Returns a ggplot object that can be visualized using the print() function
#'
#' @export
#' @examples
#'
#' data <- data.frame(x = rep(1:20, 2), group = rep(1:2, each = 20))
#' set.seed(1)
#' data$y <- (data$x^2)*data$group*3 + rnorm(40, sd = 200)
#' data$group <- ordered(data$group)
#'
#' gam <- mgcv::gam(y ~ s(x) + group, data=data)
#'
#' plot1 <- plotGAM(gamFit = gam, smooth.cov = "x", groupCovs = NULL,
#'                   rawOrFitted = "raw", plotCI=TRUE, orderedAsFactor = FALSE)
#' gam <- mgcv::gam(y ~ s(x) + group + s(x, by=group), data=data)
#' plot2 <- plotGAM(gamFit = gam, smooth.cov = "x", groupCovs = "group",
#'                              rawOrFitted = "raw", orderedAsFactor = FALSE)





plotGAM <- function(gamFit, smooth.cov ,
                    groupCovs = NULL, orderedAsFactor = T,
                    rawOrFitted = F, plotCI = T) {



  if (missing(gamFit)) {
    stop("gamFit is missing")}
  if (missing(smooth.cov)) {
    stop("smooth.cov is missing")}
  if (class(smooth.cov) != "character") {
    stop("smooth.cov must be a character")}
  if (class(groupCovs) != "character" & (!is.null(groupCovs))) {
    stop("groupCovs must be a character")}
  if (!(rawOrFitted == F | rawOrFitted == "raw" | rawOrFitted == "fitted")) {
    stop("Wrong input for rawOrFitted")
  }



  fit <- fitted <- group <- se.fit <- NULL

  if (base::class(gamFit)[1] != "gam") {
    stop("gamFit is not a object of type gam")
  }

  if (base::is.null(groupCovs)) {

    gam <- gamFit
    temp.data <- gam$model

    old.formula <- gam$formula
    old.covariates <- base::as.character(old.formula)[3]
    old.covariates <- base::gsub(" ", "", old.covariates)
    main.smooth <- base::paste0("s(", smooth.cov)
    interaction.smooth <- base::paste0("s(", smooth.cov, ",by=")

    if (!base::grepl(pattern = main.smooth, old.covariates, fixed=T)) {
      base::stop("smooth.cov is not included in formula for original model fit")
    }

    if (base::grepl(pattern = interaction.smooth, old.covariates, fixed=T)) {
      base::stop("smooth.cov has an interaction term in original model fit, pass groupCovs as an argument")
    }

    plot.df = base::data.frame( x = base::seq(min(temp.data[smooth.cov]),
                                              max(temp.data[smooth.cov]),
                                              length.out=200))
    base::names(plot.df) <- smooth.cov

    for (i in base::names(temp.data)[-1]) {
      if (i != smooth.cov) {
        if (base::any(base::class(temp.data[i][,1])[1] == c("numeric", "integer","boolean"))) {
          plot.df[, base::dim(plot.df)[2] + 1] <- base::mean(temp.data[i][,1])
          base::names(plot.df)[base::dim(plot.df)[2]] <- i
        }
        else if (base::any(base::class(temp.data[i][,1])[1] == c("character", "factor","ordered"))) {
          plot.df[, base::dim(plot.df)[2] + 1] <- temp.data[i][1,1]
          base::names(plot.df)[base::dim(plot.df)[2]] <- i
        }
        else {
          stop(base::paste("Unrecognized data type for",i,"please refit cluster model with different datatype"))
        }
      }
    }


    for (i in 1:dim(plot.df)[2]) {
      if (class(plot.df[,i])[1] == "ordered" |  class(plot.df[,i])[1] == "factor") {
        warning("There are one or more factors in the model fit, please consider plotting by group since plot might be unprecise")
      }
    }


    plot.df = base::cbind(plot.df, base::as.data.frame(mgcv::predict.gam(gam, plot.df, se.fit = TRUE)))

    if (plotCI) {
      plot <- (ggplot2::ggplot(data=plot.df, ggplot2::aes(x=plot.df[,1]))
               + ggplot2::geom_line(ggplot2::aes(y=fit), size=1)
               + ggplot2::geom_ribbon(data=plot.df, ggplot2::aes(ymax = fit+1.96*se.fit, ymin = fit-1.96*se.fit, linetype=NA), alpha = .2)
               + ggplot2::ggtitle(base::paste0(base::as.character(old.formula)[2], "  vs. s(",smooth.cov,")"))
               + ggplot2::ylab(base::paste0("Predicted ", base::as.character(old.formula)[2]))
               + ggplot2::xlab(base::paste0("s(",smooth.cov,")"))
               + ggplot2::theme_bw())
    }
    else {
      plot <- (ggplot2::ggplot(data=plot.df, ggplot2::aes(x=plot.df[,1]))
               + ggplot2::geom_line(ggplot2::aes(y=fit), size=1)
               + ggplot2::ggtitle(base::paste0(base::as.character(old.formula)[2], "  vs. s(",smooth.cov,")"))
               + ggplot2::ylab(base::paste0("Predicted ", base::as.character(old.formula)[2]))
               + ggplot2::xlab(base::paste0("s(",smooth.cov,")"))
               + ggplot2::theme_bw())
    }


    if (rawOrFitted == "raw") {
      plot <- plot + ggplot2::geom_point(data = temp.data, ggplot2::aes(x=purrr:as_vector(temp.data[smooth.cov]), y=temp.data[,1]))
      base::return(plot)
    }
    else if (rawOrFitted == "fitted") {
      temp.data$fitted <- gam$fitted.values
      plot <- plot + ggplot2::geom_point(data = temp.data, ggplot2::aes(x=purrr:as_vector(temp.data[smooth.cov]), y=fitted))
      base::return(plot)
    }
    else {
      base::return(plot)
    }

  }
  else if (base::is.character(groupCovs)) {

    gam <- gamFit
    temp.data <- gam$model

    old.formula <- gam$formula
    old.covariates <- base::as.character(old.formula)[3]
    old.covariates <- base::gsub(" ", "", old.covariates)
    main.smooth <- base::paste0("s(", smooth.cov)
    interaction.smooth <- base::paste0("s(", smooth.cov, ",by=")

    if(regexpr(groupCovs, as.character(old.covariates))[1] == 1) {
      old.covariates <- paste0("+", old.covariates)
    }

    main.group <- base::paste0("+", groupCovs)

    if (!base::grepl(pattern = main.group, old.covariates, fixed=T)) {
      base::stop("Error groupCovs is not included as a main effect in original model, include in original fit")
    }

    temp.fm <- strsplit(old.covariates, split='+', fixed = T)[[1]]
    for (i in 1:length(temp.fm)) {
      if (grepl(main.smooth, temp.fm[i], fixed=T)) {
        temp.term <- strsplit(temp.fm[i], split=',', fixed = T)[[1]]
        if (length(temp.term) > 1) {
          for (j in 2:length(temp.term)) {
            if (grepl('=',temp.term[j])) {
              if(grepl('by',temp.term[j])) {
                order <- 1:length(temp.term)
                order <- setdiff(order, j)
                  if (all(order == 1)) {
                    order <- c(1,2)
                    temp.term <- temp.term[order]
                  }
                else {
                    order <- c(order, NA)
                    order <- order[c(1, length(order), 2:(length(order) - 1))]
                    order[2] <- j
                    temp.term <- temp.term[order]
                    temp.term <- gsub(pattern=")", replacement = "", x = temp.term)
                    temp.term[length(temp.term)]<- paste0(temp.term[length(temp.term)], ")")
                  }
                }
            }
            else {
              base::stop("spline containing smooth.cov has more than one variable, refit with only one")
            }
          }
        }
        temp.fm[i] <- toString(temp.term)
      }
    }

    old.covariates <- paste0(temp.fm, collapse="+")
    old.covariates <- base::gsub(" ", "", old.covariates)

    if (!base::grepl(pattern = main.smooth, old.covariates, fixed=T) &  base::class(temp.data[groupCovs][,1])[1] == "ordered") {
      base::stop("smooth.cov is not included in formula for original model fit, model might be incorrectly parametrized")
    }
    else if (!base::grepl(pattern = interaction.smooth, old.covariates, fixed=T)) {
      base::warning("smooth.cov has no interaction term in original model fit")


      foo <-  base::seq(min(temp.data[smooth.cov]), max(temp.data[smooth.cov]), length.out=200)
      plot.df = base::data.frame( x = rep(foo, nlevels(as.factor(temp.data[groupCovs][,1]))))
      base::names(plot.df) <- smooth.cov

      plot.df$temp <- rep(levels(as.factor(temp.data[groupCovs][,1])), each = 200)
      base::names(plot.df)[2] <- groupCovs

      for (i in base::names(temp.data)[-1]) {
        if (!(i == smooth.cov | i == groupCovs)) {
          if (base::any(base::class(temp.data[i][,1])[1] == c("numeric", "integer","boolean"))) {
            plot.df[, base::dim(plot.df)[2] + 1] <- base::mean(temp.data[i][,1])
            base::names(plot.df)[base::dim(plot.df)[2]] <- i
          } else if (base::any(base::class(temp.data[i][,1])[1] == c("character", "factor","ordered"))) {
            plot.df[, base::dim(plot.df)[2] + 1] <- temp.data[i][1,1]
            base::names(plot.df)[base::dim(plot.df)[2]] <- i
          } else {
            stop(base::paste("Unrecognized data type for",i,"please refit cluster model with different datatype"))
          }
        }
      }

      plot.df = base::cbind(plot.df, base::as.data.frame(mgcv::predict.gam(gam, plot.df, se.fit = TRUE)))
      plot.df$group <- plot.df[,2]

      if (plotCI) {
        plot <- (ggplot2::ggplot(data=plot.df, ggplot2::aes(x=plot.df[,1]))
                 + ggplot2::geom_line(ggplot2::aes(y=fit, col=factor(group)), size=1)
                 + ggplot2::geom_ribbon(data=plot.df, ggplot2::aes(ymax = fit+1.96*se.fit, ymin = fit-1.96*se.fit, col=factor(group), linetype=NA), alpha = .2)
                 + ggplot2::ggtitle(base::paste0(base::as.character(old.formula)[2], " vs. s(",smooth.cov,")"))
                 + ggplot2::ylab(base::paste0("Predicted ", base::as.character(old.formula)[2]))
                 + ggplot2::xlab(base::paste0("s(",smooth.cov,")"))
                 + ggplot2::scale_colour_discrete(name = groupCovs)
                 + ggplot2::theme_bw())

      } else {
        plot <- (ggplot2::ggplot(data=plot.df, ggplot2::aes(x=plot.df[,1]))
                 + ggplot2::geom_line(ggplot2::aes(y=fit, col=factor(group)), size=1)
                 + ggplot2::ggtitle(base::paste0(base::as.character(old.formula)[2], " vs. s(",smooth.cov,")"))
                 + ggplot2::ylab(base::paste0("Predicted ", base::as.character(old.formula)[2]))
                 + ggplot2::xlab(base::paste0("s(",smooth.cov,")"))
                 + ggplot2::scale_colour_discrete(name = groupCovs)
                 + ggplot2::theme_bw())
      }

      if (rawOrFitted == "raw") {
        plot <- plot + ggplot2::geom_point(data = temp.data, ggplot2::aes(x=purrr:as_vector(temp.data[smooth.cov]), y=temp.data[,1], col=temp.data[groupCovs][,1]))
        base::return(plot)
      }
      else if (rawOrFitted == "fitted") {
        temp.data$fitted <- gam$fitted.values
        plot <- plot + ggplot2::geom_point(data = temp.data, ggplot2::aes(x=purrr:as_vector(temp.data[smooth.cov]), y=fitted, col=temp.data[groupCovs][,1]))
        base::return(plot)
      }
      else {
        base::return(plot)
      }

    }

    else {

      main.group <- base::paste0("+", groupCovs)
      if (!base::grepl(pattern = main.group, old.covariates, fixed=T)) {
        base::stop("groupCovs is not included as main effect in original model, please include in original fit")
      }

      if (base::class(temp.data[groupCovs][,1])[1] == "ordered" & orderedAsFactor) {
        base::warning("groupCovs is an ordered variable in original fit, reffiting model with groupCovs as a factor")
        temp.data[groupCovs] <- as.factor(as.character(temp.data[groupCovs][,1]))

        textCall <- gam$call
        temp.fm.updated <- strsplit(old.covariates, split='+', fixed = T)[[1]]
        index.Main <- setdiff(grep(pattern = main.smooth, x = temp.fm.updated, fixed = T),
                              grep(pattern = interaction.smooth, x = temp.fm.updated, fixed = T))
        temp.fm.updated <- temp.fm[-index.Main]
        updatedCovs <- paste0(temp.fm.updated, collapse="+")
        updatedCovs <- base::gsub(" ", "", updatedCovs)
        updatedformula <- stats::as.formula(paste(old.formula[2], old.formula[1], updatedCovs))
        new.data <- temp.data

        gam <- stats::update(gam, formula = updatedformula, data = new.data)
      }

      temp.data <- gam$model
      old.formula <- gam$formula
      old.covariates <- base::as.character(old.formula)[3]
      old.covariates <- base::gsub(" ", "", old.covariates)
      main.smooth <- base::paste0("s(", smooth.cov, ")")
      interaction.smooth <- base::paste0("s(", smooth.cov, ",by=")


      foo <-  base::seq(min(temp.data[smooth.cov]), max(temp.data[smooth.cov]), length.out=200)
      plot.df = base::data.frame( x = rep(foo, nlevels(temp.data[groupCovs][,1])))
      base::names(plot.df) <- smooth.cov

      plot.df$temp <- rep(levels(temp.data[groupCovs][,1]), each = 200)
      base::names(plot.df)[2] <- groupCovs

      for (i in base::names(temp.data)[-1]) {
        if (!(i == smooth.cov | i == groupCovs)) {
          if (base::any(base::class(temp.data[i][,1])[1] == c("numeric", "integer","boolean"))) {
            plot.df[, base::dim(plot.df)[2] + 1] <- base::mean(temp.data[i][,1])
            base::names(plot.df)[base::dim(plot.df)[2]] <- i
          } else if (base::any(base::class(temp.data[i][,1])[1] == c("character", "factor","ordered"))) {
            plot.df[, base::dim(plot.df)[2] + 1] <- temp.data[i][1,1]
            base::names(plot.df)[base::dim(plot.df)[2]] <- i
          } else {
            stop(base::paste("Unrecognized data type for",i,"please refit cluster model with different datatype"))
          }
        }
      }

      plot.df = base::cbind(plot.df, base::as.data.frame(mgcv::predict.gam(gam, plot.df, se.fit = TRUE)))
      plot.df$group <- plot.df[,2]

      if (plotCI) {
        plot <- (ggplot2::ggplot(data=plot.df, ggplot2::aes(x=plot.df[,1]))
                 + ggplot2::geom_line(ggplot2::aes(y=fit, col=factor(group)), size=1)
                 + ggplot2::geom_ribbon(data=plot.df, ggplot2::aes(ymax = fit+1.96*se.fit, ymin = fit-1.96*se.fit, col=factor(group), linetype=NA), alpha = .2)
                 + ggplot2::ggtitle(base::paste0(base::as.character(old.formula)[2], "  vs. s(",smooth.cov,")"))
                 + ggplot2::ylab(base::paste0("Predicted ", base::as.character(old.formula)[2]))
                 + ggplot2::xlab(base::paste0("s(",smooth.cov,")"))
                 + ggplot2::scale_colour_discrete(name = groupCovs)
                 + ggplot2::theme_bw())
      } else {
        plot <- (ggplot2::ggplot(data=plot.df, ggplot2::aes(x=plot.df[,1]))
                 + ggplot2::geom_line(ggplot2::aes(y=fit, col=factor(group)), size=1)
                 + ggplot2::ggtitle(base::paste0(base::as.character(old.formula)[2], "  vs. s(",smooth.cov,")"))
                 + ggplot2::ylab(base::paste0("Predicted ", base::as.character(old.formula)[2]))
                 + ggplot2::xlab(base::paste0("s(",smooth.cov,")"))
                 + ggplot2::scale_colour_discrete(name = groupCovs)
                 + ggplot2::theme_bw())
      }


      if (rawOrFitted == "raw") {
        plot <- plot + ggplot2::geom_point(data = temp.data, ggplot2::aes(x= purrr:as_vector(temp.data[smooth.cov]), y=temp.data[,1], col=temp.data[groupCovs][,1]))
        base::return(plot)
      } else if (rawOrFitted == "fitted") {
        temp.data$fitted <- gam$fitted.values
        plot <- plot + ggplot2::geom_point(data = temp.data, ggplot2::aes(x= purrr:as_vector(temp.data[smooth.cov]), y=fitted, col=temp.data[groupCovs][,1]))
        base::return(plot)
      } else {
        base::return(plot)
      }
    }

  }
  else {
    stop("groupCovs is not a character or null please correct")
  }

}


