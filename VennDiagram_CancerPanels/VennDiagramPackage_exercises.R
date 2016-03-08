# VennDiagram R package practice
# Taken from http://mfcovington.github.io/r_club/exercises/2013/05/15/venn-euler-demo/
# March 1, 2015

library(VennDiagram)

# Exercise 1
grid.newpage()
venn.plot <- draw.pairwise.venn(100, 70, 30, c("First", "Second"))
grid.draw(venn.plot)

# Exercise 2
grid.newpage()
venn.plot <- draw.pairwise.venn(area1        = 100,
                                area2        = 70,
                                cross.area   = 68,
                                scaled       = T,
                                category     = c("First", "Second"),
                                fill         = c("blue", "red"),
                                alpha        = 0.2,
                                lty          = "blank",
                                cex          = 2,
                                cat.cex      = 2,
                                cat.pos      = c(285, 105),
                                cat.dist     = 0.09,
                                cat.just     = list(c(-1, -1), c(1, 1)),
                                ext.pos      = 30,
                                ext.dist     = -0.05,
                                ext.length   = 0.85,
                                ext.line.lwd = 2,
                                ext.line.lty = "dashed")
grid.draw(venn.plot)

# Exercise 3
grid.newpage()
venn.plot <- draw.triple.venn(65, 75, 85, 35, 15, 25, 5, c("First", "Second", "Third"))
grid.draw(venn.plot)

# Exercise 4
grid.newpage()
venn.plot <- draw.triple.venn(area1    = 65,
                              area2    = 75,
                              area3    = 85,
                              n12      = 35,
                              n23      = 15,
                              n13      = 25,
                              n123     = 5,
                              category = c("First", "Second", "Third"),
                              cat.pos  = c(0, 40, 250),
                              cat.dist = c(0.05, 0.05, 0.05),
                              fill     = c("blue", "red", "green"),
                              alpha    = 0.3,
                              lty      = "blank",
                              cex      = 2,
                              cat.cex  = 2,
                              cat.col  = c("blue", "red", "green"))
grid.draw(venn.plot)