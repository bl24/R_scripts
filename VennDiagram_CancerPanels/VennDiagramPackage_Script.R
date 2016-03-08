install.packages('VennDiagram')
library(VennDiagram)

# Alex's version
# Written for 3 datasets
venn.diagram(
	x = list(
		'Foundation\nMedicine\n(315 genes)' = c(1:3,4:26,27:117,135:332),
		'TruSight\nTumor\n(26 genes)' = c(1:3,4:26),
		'NGS\nGateway\n(131 genes)' = c(4:26,27:117,118:134)
		),
	filename = "VennDiagTest",
	col = "transparent",
	fill = c("red", "blue", "green"),
	alpha = 0.5,
	label.col = c("darkred", "white", "darkblue", "white", "white", "white", "darkgreen"),
	cex = 2.5,
	fontfamily = "serif",
	fontface = "bold",
	cat.default.pos = "text",
	cat.col = c("darkred", "darkblue", "darkgreen"),
	cat.cex = 2.0,
	cat.fontfamily = "serif",
	cat.dist = c(0.36, -0.16, 0.16),
	cat.pos = -10
	)

# My version/redo
# Written for 4 datasets

# Calculating the overlaps
calculate.overlap()

# To get a list of non-overlapping genes
data.frame(panel_overlaps$a9) # Change "a" coordinates where appropriate

# Drawing the diagram
grid.newpage()
venn.plot <- draw.quad.venn(area1    = 131, #a9
                              area2    = 315, #a14
                              area3    = 26, #a1
                              area4    = 35, #a3
                              n1234    = 15, #a6
                              n123     = 23, #a12
                              n124     = 33, #a11
                              n134     = 15, #a5
                              n234     = 15, #a7
                              n12      = 114, #a15
                              n13      = 23, #a4
                              n14      = 33, #a10
                              n23      = 26, #a13
                              n24      = 35, #a8
                              n34      = 15, #a2
                              category = c("NGS\nGateway\n(131 genes)", "Foundation\nMedicine\n(315 genes)", "TruSight Tumor\n(26 genes)", "Oncomine\n(35 genes)"),
                              cat.pos  = c(-0.8, 10, 165, 1),
                              cat.dist = c(0.16, 0.15, -0.11, 0.11),
                              fill     = c("blue", "red", "green", "yellow"),
                              alpha    = 0.2,
                              scaled   = T,
                              lty      = "blank",
                              cex      = 2,
                              cat.cex  = 2,
                              cat.col  = c("black", "black", "black", "black"))
grid.draw(venn.plot)
