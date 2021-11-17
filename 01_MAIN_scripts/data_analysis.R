#Check autocorrelation of consecutive posterior estimates
#Next: run acf function
#Then: re-run with more thinning
post.all <- read_csv("02_model_validation/results/Lemon_Shark_life_history/model_validation/posterior.samples_thin1.csv")

post.thin <- read_csv("02_model_validation/results/Lemon_Shark_life_history/model_validation/posterior.samples_thin5.csv")

head(post.all) # results when thinning = 1
nrow(post.all)
allacf <- acf(post.all$Nf[1:5000], lag.max = 10)
h <- hist(post.all$Nf)


head(post.thin) # results when thinning = 5
nrow(post.thin)
thinacf <- acf(post.thin$Nf[1:5000], lag.max = 10)



#Violin plot
RB_violin <- ggplot(data=post.all, aes(x=factor(iter))) +
  geom_violin(aes(y=Nf)) +
  ylim(-50, 160) +
  geom_hline(yintercept=0, col="black", size=1.25) +
  annotate("rect", xmin=0, xmax=Inf, ymin=-20, ymax=20, alpha=.5, col="red") +
  labs(x="Iteration", y="Relative bias", title="Relative bias by iteration") +
  scale_fill_brewer(palette="Set2") +
  font("title", size = 10, face = "bold")
