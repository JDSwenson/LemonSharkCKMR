plot.Nfb <- function(data){
  Nfb.boxplot <- data %>% dplyr::filter(parameter %in% c("Nfbt")) %>% 
    ggplot() +
    geom_boxplot(aes(x=factor(estimation.year), 
                     fill = time_window, 
                     y=relative_bias), 
                 notch = TRUE, 
                 outlier.alpha = 0.1) +
    geom_hline(yintercept=0, 
               col="blue", 
               size=1.25, 
               linetype = "dashed") +
    labs(y="Relative bias", 
         x = "Year") +
    scale_y_continuous(limits = c(-50, 150),
                       breaks = c(-50, -25, 0, 25, 50, 75, 100, 125, 150)) +
    scale_fill_manual(values = cp5) +
    geom_vline(xintercept=7, col="darkslategrey", size=1.25, linetype = "longdash") +
    geom_vline(xintercept=12, col="darkslategrey", size=1.25, linetype = "longdash") +
    #facet_wrap(~parameter, scales = "free", ncol = 2) +
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          #axis.title.y = element_text(size = 20),
          axis.title.y = element_blank(),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 20),
          #legend.position = "none",
          plot.title = element_text(size=25),
          axis.title.y.right = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    guides(fill=guide_legend(title="Sample window")) +
    ggtitle("a) Nfb")
  
  Nfb.boxplot
}

plot.survival <- function(data){
  survival.boxplot <- data %>% dplyr::filter(parameter %in% c("survival")) %>% 
    ggplot() +
    geom_boxplot(aes(x=factor(estimation.year), fill = time_window, y=relative_bias),
                 notch = TRUE, 
                 outlier.alpha = 0.1) +
    geom_hline(yintercept=0, col="blue", size=1.25, linetype = "dashed") +
    labs(y="Relative bias", x = "Year") +
    scale_y_continuous(limits = c(-50, 50), breaks = c(-50, -25, 0, 25, 50),
                       sec.axis = sec_axis(~. *1, name = "abundance")
    ) +
    scale_fill_manual(values = cp5) +
    geom_vline(xintercept=7, col="darkslategrey", size=1.25, linetype = "longdash") +
    geom_vline(xintercept=12, col="darkslategrey", size=1.25, linetype = "longdash") +
    #facet_wrap(~parameter, scales = "free", ncol = 2) +
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          #axis.title = element_text(size = 20),
          axis.title.y = element_blank(),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 20),
          plot.title = element_text(size=25),
          axis.title.y.right = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank()) +
    guides(fill=guide_legend(title="Sample window")) +
    ggtitle(expression(paste("c) ", phi)))
  
  survival.boxplot
}

plot.lambda <- function(data){
  lambda.boxplot <- data %>% dplyr::filter(parameter %in% c("lambda")) %>% 
    ggplot() +
    geom_boxplot(aes(x=factor(estimation.year), fill = time_window, y=relative_bias),
                 notch = TRUE, 
                 outlier.alpha = 0.1) +
    geom_hline(yintercept=0, col="blue", size=1.25, linetype = "dashed") +
    labs(y="Relative bias", x = "Year") +
    scale_y_continuous(limits = c(-25, 50), breaks = c(-25, -25, 0, 25, 50),
                       sec.axis = sec_axis(~. *1, name = "abundance")
    ) +
    geom_path(data = ls_dataTrend_allSibs.df, #Include if we want a trend line of the average realized abundance
              aes(x = factor(year),
                  y = Nfb.trend.trans,
                  color = type,
                  group = 1),
              size = 1.25,
              linetype = "solid") +
    scale_color_manual(values = "firebrick4", name = element_blank()) +
    # facet_wrap(~parameter,
    #            scales = "free") +
    scale_fill_manual(values = cp5) +
    #facet_wrap(~parameter, scales = "free", ncol = 2) +
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          #axis.title = element_text(size = 20),
          axis.title.y = element_blank(),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 20),
          plot.title = element_text(size=25),
          axis.title.y.right = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank()) +
    geom_vline(xintercept=7, col="darkslategrey", size=1.25, linetype = "longdash") +
    geom_vline(xintercept=12, col="darkslategrey", size=1.25, linetype = "longdash") +
    guides(fill=guide_legend(title="Sample window")) +
    ggtitle(expression(paste("b) ", lambda)))
  
  lambda.boxplot
}

plot.single.iter <- function(data){
  single.iter.plot <- data %>% dplyr::filter(iteration == iter,
                                                               parameter %in% params) %>% 
    ggplot(aes(x = estimation.year, 
               y = estimate,
               color = time_window)) +
    geom_pointrange(aes(ymin = HPD2.5, 
                        ymax = HPD97.5),
                    position = position_dodge(width = 0.5)) +
    geom_smooth(method = "loess", 
                se = FALSE, 
                aes(linetype = time_window)) +
    geom_line(data = ls_truthSim.abund.Plot_allSibs,
              aes(y = truth, 
                  color = type),
              size = 1) +
    scale_x_continuous(breaks = seq(from = 74, to=90, by=1)) +
    # geom_path(data = ls_dataTrend.df, #Include if we want a trend line of the average realized abundance
    #           aes(x = year,
    #               y = Nfb.trend.trans,
    #               color = type,
    #               group = 1),
    #           size = 1.25,
    #           linetype = "solid") +
    #  facet_wrap(~parameter, 
    #             scales = "free") +
    theme_bw() +
    #  labs(x = "Estimation year", y = "Female abundance") +
    scale_color_manual(name = "Sample window", values = cp5) +
    #  ggtitle("a)") +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 20),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 20),
          plot.title = element_text(size=25)) +
    xlab("Year") +
    ylab("Female abundance") +
    guides(linetype = "none") #Can decide what to keep or remove from legend with this 
  
  single.iter.plot
}

plot.psi <- function(data){
  psi.boxplot <- data %>% dplyr::filter(parameter %in% c("psi")) %>% 
    ggplot() +
    geom_boxplot(aes(x=factor(estimation.year), fill = time_window, y=Q50),
                 notch = TRUE, 
                 outlier.alpha = 0.1) +
    geom_hline(yintercept=.83, col="blue", size=1.25, linetype = "dashed") +
    labs(y="Relative bias", x = "Year") +
    # scale_y_continuous(limits = c(-50, 50), breaks = c(-50, -25, 0, 25, 50),
    #                    sec.axis = sec_axis(~. *1, name = "abundance")
    # ) +
    scale_fill_manual(values = cp5) +
    geom_vline(xintercept=7, col="darkslategrey", size=1.25, linetype = "longdash") +
    geom_vline(xintercept=12, col="darkslategrey", size=1.25, linetype = "longdash") +
    #facet_wrap(~parameter, scales = "free", ncol = 2) +
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          #axis.title = element_text(size = 20),
          axis.title.y = element_blank(),
          strip.text.x = element_text(size = 15),
          strip.text.y = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 20),
          plot.title = element_text(size=25),
          axis.title.y.right = element_blank(),
          axis.text.y.right = element_blank(),
          axis.ticks.y.right = element_blank()) +
    guides(fill=guide_legend(title="Sample window")) +
    ggtitle(expression(paste("c) ", phi)))
  
  psi.boxplot
}
