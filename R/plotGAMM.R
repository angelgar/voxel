#' GAMM plotting using ggplot2
#'
#' @family Plotting
#' @param gammFit fitted gam model as produced by gamm4::gamm()
#' @param smooth.cov (character) name of smooth term to be plotted
#' @param groupCovs (character)  name of group variable to plot by, if NULL (default) then there are no groups in plot
#' @param orderedAsFactor Disabled
#' @param rawOrFitted If FALSE (default) then only smooth terms are plotted; if rawOrFitted = "raw" then raw values are plotted against smooth; if rawOrFitted = "fitted" then fitted values are plotted against smooth
#' @param plotCI if TRUE (default) upper and lower confidence intervals are added at 2 standard errors above and below the mean
#' @param grouping (character) Name of variable that you want to use as the group argument in ggplot2::aes(), useful for better visualization of longitudinal data, (default is NULL)
#'
#' @return Returns a ggplot object that can be visualized using the print() function
#'
#' @export
#' @examples
#'
#' set.seed(1)
#' data <- data.frame(x = (seq(.25,25, .25) +rnorm(100)), group = rep(1:2, 5), z=rnorm(100),
#'               index.rnorm = rep(rnorm(50, sd = 50), 2), index.var = rep(1:50, 2))
#' data$y <- (data$x)*data$group*10 + rnorm(100, sd = 700) + data$index.rnorm + data$z
#' data$group <- ordered(data$group)
#'
#'
#' gamm <- gamm4::gamm4(y ~ + s(x) + s(x, by=group) + z + group, data=data, random = ~ (1|index.var))
#'
#'
#' plot <- plotGAMM(gammFit <- gamm, smooth.cov <- "x", groupCovs = "group",
#'                     plotCI <- T, rawOrFitted = "raw", grouping = "index.var")
#'
#' plot2 <- plotGAMM(gammFit <- gamm, smooth.cov <- "x", groupCovs = "group",
#'                   plotCI <- T, rawOrFitted = "fitted", grouping = "index.var")





plotGAMM <- function(gammFit, smooth.cov , groupCovs = NULL, orderedAsFactor = F, rawOrFitted = F, plotCI = T, grouping = NULL) {

  if (missing(gammFit)) {
    stop("gammFit is missing")}
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

  if (!(base::class(gammFit$gam) == "gam" & base::class(gammFit$mer) == "lmerMod")) {
    stop("gamFit is not a object created by gamm4()")
  }



  if (base::is.character(rawOrFitted)) {
    if (((rawOrFitted == "raw") & base::is.null(grouping))) {
      base::warning("grouping variable is missing will plot raw values as points rather than trajectories")
    } else if (((rawOrFitted == "fitted") & base::is.null(grouping))) {
      base::warning("grouping variable is missing will plot fitted values as points rather than trajectories")
    }
  }

  if (base::is.null(rawOrFitted) & !(is.null(grouping))) {
    warning("Ignoring grouping variable and defaulting to no raw or fitted values in plot")
  }


  if (base::is.null(groupCovs)) {

    gam <- gammFit$gam
    temp.data <- gam$model

    if (is.character(grouping)) {
      if (!any(names(temp.data) == grouping) ) {
        stop("grouping variable is not in model fit, please update model arguments")
      }
    }

    if (is.character(grouping)) {
      temp.data[grouping] <- as.factor(temp.data[grouping][,1])
    }


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

    plot.df = base::data.frame( x = base::seq(min(temp.data[smooth.cov]), max(temp.data[smooth.cov]), length.out=200))
    base::names(plot.df) <- smooth.cov

    for (i in names(gam$var.summary)) {
      if (i != smooth.cov) {
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
    } else {
      plot <- (ggplot2::ggplot(data=plot.df, ggplot2::aes(x=plot.df[,1]))
               + ggplot2::geom_line(ggplot2::aes(y=fit), size=1)
               + ggplot2::ggtitle(base::paste0(base::as.character(old.formula)[2], "  vs. s(",smooth.cov,")"))
               + ggplot2::ylab(base::paste0("Predicted ", base::as.character(old.formula)[2]))
               + ggplot2::xlab(base::paste0("s(",smooth.cov,")"))
               + ggplot2::theme_bw())
    }

    if (base::is.character(rawOrFitted)) {
      if (rawOrFitted == "raw" & base::is.null(grouping)) {
        plot <- plot + ggplot2::geom_point(data = temp.data, ggplot2::aes(x = as_vector(temp.data[smooth.cov]), y=temp.data[,1]))
        base::return(plot)

      } else if (rawOrFitted == "fitted" & base::is.null(grouping)) {
        temp.data$fitted <- gam$fitted.values
        plot <- plot + ggplot2::geom_point(data = temp.data, ggplot2::aes(x = as_vector(temp.data[smooth.cov]), y=fitted))
        base::return(plot)

      } else if (rawOrFitted == "raw" & base::is.character(grouping)) {
        plot <- plot + ggplot2::geom_line(data = temp.data, ggplot2::aes(x = as_vector(temp.data[smooth.cov]), y=temp.data[,1], group=factor(temp.data[grouping][,1])), alpha=.2)
        base::return(plot)

      } else if (rawOrFitted == "fitted" & base::is.character(grouping)) {
        temp.data$fitted <- gam$fitted.values
        plot <- plot + ggplot2::geom_line(data = temp.data, ggplot2::aes(x = as_vector(temp.data[smooth.cov]), y=fitted, group=factor(temp.data[grouping][,1])), alpha=.2)
        base::return(plot)

      }
    } else {
      base::return(plot)
    }

  }

  else if (base::is.character(groupCovs)) {

    gam <- gammFit$gam
    temp.data <- gam$model

    if (is.character(grouping)) {
      if (!any(names(temp.data) == grouping) ) {
        stop("grouping variable is not in model fit, please update model arguments")
      }
    }

    if (is.character(grouping)) {
      temp.data[grouping] <- as.factor(temp.data[grouping][,1])
    }


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
                } else {
                  order <- c(order, NA)
                  order <- order[c(1, length(order), 2:(length(order) - 1))]
                  order[2] <- j
                  temp.term <- temp.term[order]
                  temp.term <- gsub(pattern=")", replacement = "", x = temp.term)
                  temp.term[length(temp.term)]<- paste0(temp.term[length(temp.term)], ")")
                }
              }
            } else {
              base::stop("spline containing smooth.cov has more than one variable, refit with only one")
            }
          }
        }
        temp.fm[i] <- toString(temp.term)
      }
    }

    old.covariates <- paste0(temp.fm, collapse="+")
    old.covariates <- base::gsub(" ", "", old.covariates)


    if (!base::grepl(pattern = main.smooth, old.covariates, fixed=T) &
        base::class(temp.data[groupCovs][,1])[1] == "ordered") {
      base::stop("smooth.cov is not included in formula for original model fit, model might be incorrectly parametrized")
    } else if (!base::grepl(pattern = interaction.smooth, old.covariates, fixed=T)) {
      base::warning("smooth.cov has no interaction term in original model fit")


      foo <-  base::seq(min(temp.data[smooth.cov]), max(temp.data[smooth.cov]), length.out=200)
      plot.df = base::data.frame( x = rep(foo, nlevels(as.factor(temp.data[groupCovs][,1]))))
      base::names(plot.df) <- smooth.cov

      plot.df$temp <- rep(levels(as.factor(temp.data[groupCovs][,1])), each = 200)
      base::names(plot.df)[2] <- groupCovs

      for (i in names(gam$var.summary)) {
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


      if (base::is.character(rawOrFitted)) {
        if (rawOrFitted == "raw" & base::is.null(grouping)) {
          plot <- plot +
                  ggplot2::geom_point(data = temp.data,
                                      ggplot2::aes(x = as_vector(temp.data[smooth.cov]),
                                                   y=temp.data[,1],
                                                   col=factor(temp.data[groupCovs][,1])))
          base::return(plot)

        } else if (rawOrFitted == "fitted" & base::is.null(grouping)) {
          temp.data$fitted <- gam$fitted.values
          plot <- plot +
                  ggplot2::geom_point(data = temp.data,
                                      ggplot2::aes(x = as_vector(temp.data[smooth.cov]),
                                                   y=fitted,
                                                   col=factor(temp.data[groupCovs][,1])))
          base::return(plot)

        } else if (rawOrFitted == "raw" & base::is.character(grouping)) {
          plot <- plot +
                  ggplot2::geom_line(data = temp.data,
                                     ggplot2::aes(x = as_vector(temp.data[smooth.cov]),
                                                  y=temp.data[,1],
                                                  group=factor(temp.data[grouping][,1]),
                                                  col=factor(temp.data[groupCovs][,1])),
                                     alpha=.4)
          base::return(plot)

        } else if (rawOrFitted == "fitted" & base::is.character(grouping)) {
          temp.data$fitted <- gam$fitted.values
          plot <- plot +
                  ggplot2::geom_line(data = temp.data,
                                     ggplot2::aes(x = as_vector(temp.data[smooth.cov]),
                                                  y=fitted,
                                                  group=factor(temp.data[grouping][,1]),
                                                  col=factor(temp.data[groupCovs][,1])),
                                     alpha=.4)
          base::return(plot)
        }
      } else {
        base::return(plot)
      }

    } else {

      main.group <- base::paste0("+", groupCovs)
      if (!base::grepl(pattern = main.group, old.covariates, fixed=T)) {
        base::stop("groupCovs is not included as main effect in original model, please include in original fit")
      }

      orderedAsFactor <- F


      if (base::class(temp.data[groupCovs][,1])[1] == "ordered" & orderedAsFactor) {
        base::warning("groupCovs is an ordered variable in original fit, reffiting model with groupCovs as a factor")
        temp.data[groupCovs] <- as.factor(as.character(temp.data[groupCovs][,1]))

        textCall <- gam$call
        updatedCovs <- base::gsub(main.smooth, replacement = "", old.covariates, fixed = T)
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

      for (i in names(gam$var.summary)) {
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


      if (base::is.character(rawOrFitted)) {
        if (rawOrFitted == "raw" & base::is.null(grouping)) {
          plot <- plot + ggplot2::geom_point(data = temp.data, ggplot2::aes(x = as_vector(temp.data[smooth.cov]), y=temp.data[,1], col=factor(temp.data[groupCovs][,1])))
          base::return(plot)

        } else if (rawOrFitted == "fitted" & base::is.null(grouping)) {
          temp.data$fitted <- gam$fitted.values
          plot <- plot + ggplot2::geom_point(data = temp.data, ggplot2::aes(x = as_vector(temp.data[smooth.cov]), y=fitted, col=factor(temp.data[groupCovs][,1])))
          base::return(plot)

        } else if (rawOrFitted == "raw" & base::is.character(grouping)) {
          plot <- plot + ggplot2::geom_line(data = temp.data, ggplot2::aes(x = as_vector(temp.data[smooth.cov]), y=temp.data[,1], group=factor(temp.data[grouping][,1]), col=factor(temp.data[groupCovs][,1])), alpha=.4)
          base::return(plot)

        } else if (rawOrFitted == "fitted" & base::is.character(grouping)) {
          temp.data$fitted <- gam$fitted.values
          plot <- plot + ggplot2::geom_line(data = temp.data, ggplot2::aes(x = as_vector(temp.data[smooth.cov]), y=fitted, group=factor(temp.data[grouping][,1]), col=factor(temp.data[groupCovs][,1])), alpha=.4)
          base::return(plot)

        }
      } else {
        base::return(plot)
      }

    }
  }

  else {
    stop("groupCovs is not a character or null please correct")
  }
}

