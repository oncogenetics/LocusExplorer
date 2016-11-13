

df <- data.frame(x1 = 2:4, x2 = 4:6, y1 = 0, y2 = 0, R2 = 1:3)
ggplot(mtcars, aes(wt, mpg)) +
  geom_point() + ylim(0,40) +
  geom_curve(aes(x = x1, y = y1, xend = x2, yend = y2, size = R2),
             alpha = 0.3, col = "red",
             curvature = -1, data = df)
  