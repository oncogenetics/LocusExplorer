library(ggplot2)

ggplot(data=data.frame(x=c(1,1,1,1),
                       y=4:1,
                       fill=c("black","black","green","purple"),
                       shape=c(2,1,15,15),
                       label=c("Typed","Imputed","Positive","Negative"),
                       stringsAsFactors = FALSE),
       
       aes(x=x,y=y,col=fill,fill=fill,shape=shape,label=label)) +
  geom_point(size=4) +
  geom_text(col="black",hjust=-0.2, vjust=0.5) +
  #annotate("text",x=1.5,y=4:1,) +
  #xlim(0,5) + 
  scale_fill_identity() +
  scale_color_identity() +
  scale_shape_identity() +
  theme_null() 
