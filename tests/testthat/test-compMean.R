### test-compMean.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jan  3 2023 (09:44) 
## Version: 
## Last-Updated: jan  3 2023 (10:30) 
##           By: Brice Ozenne
##     Update #: 2
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

d.test <- structure(list(Y1 = c(-0.59, 0.03, -1.52, -1.36, 1.18, -0.93, 
                                1.32, 0.62, -0.05, -1, 0.66, 1.65, -1.26, 0.88, 0.74, 0.01, -0.36, 
                                -0.53, -1.04, -1.27, 1.34, -1.04, 1.18, -2.8, 1.04, 0.22, 0.63, 
                                -0.96, -1.66, 2.07, 0.33, -0.82, 0.35, -0.91, -2.23, 0.23, 0.97, 
                                0.62, 0.07, 1.02),
                         Y2 = c(-0.44, 0.24, 0.6, -0.12, -2.07, 0.59, 
                                0.49, -1.01, 1.27, 1.12, -1.55, -0.73, -0.8, 2.07, 2.05, 0.79, 
                                0.18, -0.68, 1.46, 0.44, 0.95, -1.34, 0.41, 0.93, -0.93, 0.61, 
                                -0.76, 0.26, -0.1, -1.57, -0.03, -0.88, -0.27, 1.37, -0.52, 0.1, 
                                0.56, 1.22, -0.55, 1.51),
                         group = c("C", "C", "C", "C", "C", 
                                   "C", "C", "C", "C", "C", "T", "T", "T", "T", "T", "T", "T", "T", 
                                   "T", "T", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "T", 
                                   "T", "T", "T", "T", "T", "T", "T", "T", "T"),
                         time = c("D1", 
                                  "D1", "D1", "D1", "D1", "D1", "D1", "D1", "D1", "D1", "D1", "D1", 
                                  "D1", "D1", "D1", "D1", "D1", "D1", "D1", "D1", "W1", "W1", "W1", 
                                  "W1", "W1", "W1", "W1", "W1", "W1", "W1", "W1", "W1", "W1", "W1", 
                                  "W1", "W1", "W1", "W1", "W1", "W1")),
                    row.names = c(NA, -40L), class = c("data.frame"))

d.test2 <- d.test[c(13L, 9L, 32L, 31L, 11L, 15L, 18L, 30L, 20L, 12L, 26L, 37L, 
                    38L, 17L, 4L, 8L, 40L, 29L, 2L, 24L, 1L, 7L, 35L, 3L, 5L, 23L, 
                    6L, 36L, 16L, 10L, 21L, 19L, 22L, 14L, 25L, 34L, 28L, 27L, 39L, 
                    33L),]

res.work <- compMean(Y = cbind(Y1 = d.test$Y1, Y2 = d.test$Y2), 
                     group = d.test$group, time = d.test$time,
                     permutation.type = "group")

res.fail <- compMean(Y = cbind(Y1 = d.test2$Y1, Y2 = d.test2$Y2), 
                     group = d.test2$group, time = d.test2$time,
                     permutation.type = "group")
## previous error message
## Y[(group == Ugroup[1]) * (time == Utime[1]) == 1, , drop = FALSE] : 
##   (subscript) logical subscript too long
## In addition: Warning message:
## In (group == Ugroup[1]) * (time == Utime[1]) :
##   longer object length is not a multiple of shorter object length




##----------------------------------------------------------------------
### test-compMean.R ends here
