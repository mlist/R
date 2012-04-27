p <- qplot(Transfection, CTblue, data=kIn, colour=InnerOuter, shape=TopBottom)
p <- p + facet_wrap(~cellLine) + opts(axis.text.x=theme_text(angle=45))

p <- qplot(plateColumn, plateRow, data=blanc) + geom_tile(aes(fill = CTblue), colour = "white") 
p <- p + scale_fill_gradient(low = "white", high = "steelblue")
p <- p + facet_wrap(~cellLine)
p <- p + scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0))

p <- qplot(Transfection, CTblue, data=kIn) + geom_boxplot(aes(fill = factor(InnerOuter))) 
p <- p + facet_wrap(~cellLine) + opts(axis.text.x=theme_text(angle=45))