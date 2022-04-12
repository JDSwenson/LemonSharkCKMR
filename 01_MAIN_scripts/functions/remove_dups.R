split.dups <- function(samples){

  #Identify duplicates
dups.first.capture <- samples %>% group_by(indv.name) %>% 
  filter(n() > 1) %>% 
  distinct(indv.name, .keep_all = TRUE) %>% 
  ungroup() 
  
dups.later.capture <- samples %>% group_by(indv.name) %>% 
  filter(n() > 1) %>% 
  anti_join(dups.first.capture, by = c("indv.name", "capture.year")) %>% 
  ungroup()

#Save the first capture instance of recaptured individuals and sort so older individual comes first (for pairwise comparison matrix)
samples.save.first <- samples %>% anti_join(dups.later.capture, by = c("indv.name", "capture.year")) %>% 
  dplyr::arrange(birth.year) %>% 
  as_tibble()

#Save the last capture instance of recaptured individuals and sort so older individual comes first (for pairwise comparison matrix)
samples.save.second <- samples %>% anti_join(dups.first.capture, by = c("indv.name", "capture.year")) %>%
  dplyr::arrange(birth.year) %>% 
  as_tibble()

return(list(samples.save.first, samples.save.second))  
}