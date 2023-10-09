#------------------- Simulation figs --------------------------
brewer.pal.info #See different palette options from Rcolorbrewer
cp5 <- brewer.pal(3, "Dark2") #Choose palette + number of colors
cp5[4] <- "black" #If not including average trendline on second plot
#cp5[4] <- "darkslateblue" #If including average trendline on second plot
#cp5[5] <- "black" #If including average trendline on second plot
scales::show_col(cp5)


ls_data <- results

x.ticks <- c(seq(1997, 2015, by = 1))
x.tick.labels <- c("", "", "", "2000", rep("", times = 4), "2005", rep("", times = 4), "2010", rep("", times = 4), "2015")

fig5a <- ls_data %>% 
  dplyr::filter(parameter == "Nfbt") %>%
  ggplot(aes(x = estimation.year, y = Q50, color = time_window)) + 
  #geom_point() + 
  geom_pointrange(aes(ymin = HPD2.5, ymax = HPD97.5),
                  position = position_dodge(width = 0.5)) +
  geom_smooth(method = "gam", se = FALSE, aes(linetype = time_window), show.legend = FALSE) +
  #facet_wrap(~years.lab) +
  theme_bw() +
  ylim(0, 75) +
  labs(x = "Estimation year", y = "Female abundance") +
  scale_color_manual(values = cp5, name = "Sample window") +
  #  ggtitle("a)") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 20),
        plot.title = element_text(size=25)) +
  scale_x_continuous(breaks = x.ticks,
                     labels = x.tick.labels)

fig5a
