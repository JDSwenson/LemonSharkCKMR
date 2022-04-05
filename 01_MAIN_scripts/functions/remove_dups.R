split.dups <- function(samples){

dups.first.capture <- samples %>% group_by(indv.name) %>% 
  filter(n() > 1) %>% 
  distinct(indv.name, .keep_all = TRUE) %>% 
  ungroup() 
  

dups.later.capture <- samples %>% group_by(indv.name) %>% 
  filter(n() > 1) %>% 
  anti_join(dups.first.capture, by = c("indv.name", "capture.year")) %>% 
  ungroup()

#capture.yrs <- cbind(kept = dups.kept$capture.year, removed = dups.removed$capture.year)

samples.save.first <- samples %>% anti_join(dups.later.capture, by = c("indv.name", "capture.year")) %>% 
  dplyr::arrange(birth.year) %>% 
  as_tibble()

samples.save.second <- samples %>% anti_join(dups.first.capture, by = c("indv.name", "capture.year")) %>%
  dplyr::arrange(birth.year) %>% 
  as_tibble()

return(list(samples.save.first, samples.save.second))  
}