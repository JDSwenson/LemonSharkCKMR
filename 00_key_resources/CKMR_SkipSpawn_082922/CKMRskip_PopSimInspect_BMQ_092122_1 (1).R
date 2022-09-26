library(tidyverse)
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5))
library(gridExtra)

simdir = "D:\\Google Drive\\Research\\Manuscripts\\Swenson CKMR\\Simulations\\Population.simulations"

simfiles = list.files(simdir, pattern = "sample.info.*ALL")
simfiles.psi = sapply(strsplit(simfiles[], "_"), `[`, 6)

plots.out = list()
i=1
for(i in seq(length(simfiles)))
{
  simfile = simfiles[i]
    
  sample = readRDS(paste0(simdir,"/",simfile))
  
  sample.sub = sample %>% filter(birth.year == 90)
  
  sample.sub$birth.year.mom = sample$birth.year[match(sample.sub$mother.x, sample$indv.name)]
  sample.sub$repro.strategy.mom = sample$repro.strategy[match(sample.sub$mother.x, sample$indv.name)]
  sample.sub$mom = sample$indv.name[match(sample.sub$mother.x, sample$indv.name)]
  sample.sub$repro.cycle.mom = sample$repro.cycle[match(sample.sub$mother.x, sample$indv.name)]
  
  
  plots.out[[i]] = ggplot(sample.sub)+
    geom_histogram(aes(x = birth.year.mom, fill = repro.strategy.mom), binwidth = 1, col = 1)+

    scale_x_continuous(breaks = seq(100))+
    theme(axis.text.x = element_text(size = 6))+
    ggtitle(paste0(simfiles.psi[i],"\ncohort year: ", max(sample$birth.year),", capture year: ",max(sample$capture.year)))
}

n <- length(plots.out)
nCol <- floor(sqrt(n))

do.call("grid.arrange", c(plots.out, ncol=nCol))

png("D:\\Google Drive\\Research\\Manuscripts\\Swenson CKMR\\Simulations\\CKMRskip_PopSimInspect_BMQ_Plot_092122_1.png",
    res = 600, width = 14, height = 6, units = 'in')
do.call("grid.arrange", c(plots.out, ncol=nCol))
dev.off()


