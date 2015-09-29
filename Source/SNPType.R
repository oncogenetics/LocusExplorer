# Author: Tokhir Dadaev
# License: MIT + file LICENSE.txt

ggplot(data=ROIdatStats(),
       aes(x=BP, xend=BP, 
           y=TYPED-1,yend=TYPED, 
           colour=ifelse(TYPED==2,"#003333","#669999"))) +
  geom_segment() +
  scale_color_identity() +
  scale_y_continuous(breaks=c(0.5,1.5),
                     labels=udf_pad(c("Imputed","Typed")),
                     name="SNP") +  
  #testing plot alignment
  #geom_vline(xintercept=173000000,col="red") +
  xlim(c(zoomStart(),zoomEnd())) +
  udf_theme() +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())
