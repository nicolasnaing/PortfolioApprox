# PortfolioApprox
The broader question for this project was “Given data for 3 assets, how can we find the allotment for each asset that maximizes return?” We can capture the nature of these assets through a 1x3 mean vector mu and 3x3 matrix sigma for the pairwise covariances of the 3 assets. Using a 1x3 vector for allotments of the assets, the return is computed as mu^t w -gamma w^T sigma w, where gamma is a number from 0 to 1 that serves as a risk aversion factor (high gamma -> high risk aversion -> avoid assets with high variance). 

By choosing distributions for the 3 assets, we can find the exact value of mu and sigma, which have their weights that maximize return (W1). The goal is to take data for the 3 assets, approximations mu,sigma, (called mu2, sigma2) and find the weights (W2) that maximize return for mu2,sigma2. I found these weights by using a formula from a book Prof. Cornuejols wrote on optimization. I then computed two sets of returns: The first return is computed with mu,sigma, W1 and the second is computed with mu,sigma, W2. The difiference between these two values is the main source of data, and it is better to have this difference be as small as possible. 

My work specifically focused on how multiplying the covariance matrix by a constant from 0 to 1 affected the difference between the two computed returns. This multiplication occurred after building mu2, sigma2, so if we call the adjusted sigma sigma3, then I used the formula on mu2 and sigma3 to find a set of weights W3. I then computed two returns again: The first return is computed with mu,sigma, W1 and the second is computed with mu,sigma, W3.

I conducted this research from February 2020 – May 2020 with Prof. Cornuejols and used MATLAB to collect my data.

# List of file names (The pdf file ShrinkageCov.pdf contains the full explanation for the work and the "questions" that the below descriptions refer to
portfolioapprox is tex files

portfolioapprox1 is code for first question

portfolioapprox2 is code for first part of second question (only tested .5 .6 .7 .8 .9)

unifgroupresults is code for second part of second question, first set of rvs with 0.00 to 1.00 coefficients

correlated is code for third part of second question, second set of rvs (uncorrelated)

results1 contains results for first question 

results2 is for first + second part of second question

results3 are for third part of second question
