library(gRbase)
library(gtools)

v1 <- c(1:4)

combinations(samples.df.all$indv.name, n=length(samples.df.all$indv.name), r=2, repeats.allowed = FALSE) %>% 
  as_tibble()

c.df <- crossing(v1 = v1, v2 = v1) %>% 
  mutate(c1 = paste0(v1, "_", v2),
         c2 = paste0(v2, "_", v1))


data.frame(t(fastcombn(v1, m=2))) %>% 
  as_tibble()

data.frame(t(fastcombn(filtered.samples.HS.df$indv.name, m=2))) %>% 
  as_tibble()
