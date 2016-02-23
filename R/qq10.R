#
# functions to draw QQ plots on -log10 scale
#

qq10.new <- function (pmin) {
    plot.new()
    plot.window(c(0, -log10(pmin)), c(0, -log10(pmin)))
    axis(1)
    axis(2, las = 1)
    title(xlab = expression(-log[10](plain(expected) ~ ~paste(italic(P), 
        "-value"))), line = 2)
    title(ylab = expression(-log[10](plain(observed) ~ ~paste(italic(P), 
        "-value"))), line = 2)
    box()
}

qq10.envelope <- function (n, pmin, alpha = 0.01, col = "grey") {
    pre <- -log10(pmin) + 2
    pp <- 10^(-seq(from = 0.01, to = pre, by = 0.01))
    qqenv <- data.frame(x = c(0), y = c(0))
    ppq <- pp + sqrt(pp * (1 - pp)/n) * qnorm(alpha/2)
    ppq <- ifelse(0 < ppq & ppq < 1, ppq, NA)
    qqenv <- rbind(qqenv, na.omit(data.frame(x = -log10(ppq), 
        y = -log10(pp))))
    qqenv <- rbind(qqenv, data.frame(x = c(pre), y = max(-log10(pp[!is.na(ppq)]))))
    qqenv <- rbind(qqenv, data.frame(x = c(pre), y = c(pre)))
    pp <- 10^(-seq(from = pre, to = 0.01, by = -0.01))
    ppq <- pp + sqrt(pp * (1 - pp)/n) * qnorm(1 - alpha/2)
    ppq <- ifelse(0 < ppq & ppq < 1, ppq, NA)
    qqenv <- rbind(qqenv, na.omit(data.frame(x = -log10(ppq), 
        y = -log10(pp))))
    polygon(qqenv, col = col, border = col)
}

qq10.points <- function (p, ...) {
  pe <- (rank(p, na.last = "keep")-0.5) / length(na.omit(p))
  points(-log10(pe), -log10(p), ...)
}

qq10 <- function (p, pmin = NULL, alpha = 0.01, ...) {
  p <- p[!is.na(p)]
  stopifnot(length(p) > 0)
  if (is.null(pmin)) pmin <- min(c(p, 0.5/length(p)))
  stopifnot(pmin > 0)
  qq10.new(pmin)
  if (all(!is.na(alpha))) qq10.envelope(length(p), pmin, alpha = alpha, col = "grey")
  qq10.points(p, ...)
  box()
}

gclambda <- function(p, df = 1, prob = 0.5) {
  return(unname(qchisq(quantile(p, probs = prob, na.rm = TRUE), df = df, lower.tail = FALSE) /
                qchisq(prob, df = df, lower.tail = FALSE)))
}
