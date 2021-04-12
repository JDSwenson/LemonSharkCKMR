02/08/2021
On 2/3/2021 I ran fishSim on a loop and fit a CKMR model when ages were correctly assigned or misassigned by a small amount.
On 2/5/2021 I ran fishSim using a truncated poisson distribution in the call to altMate. The estimates were still very positively biased.
Today (2/8/2021) I realized that there is some slop in my model, where there is a non-zero probability for pairs of animals that could not possibly be related. I also realize that by estimating abundance for the first POSITIVE comparison, rather than the year of the oldest animal, there are a bunch of negative comparisons that were receiving non-zero probabilities.
