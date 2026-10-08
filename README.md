# A Latent Editing Mixture Model for Clustering DNA Reads

This repository contains the reproducibility materials for the article Dallari, Redivo, Viroli (2026+) *A Latent Editing Mixture Model for Clustering DNA Reads*. 
It provides the code needed to reproduce the computational results and figures in the paper.

The proposed model is implemented in the function `ledm_em` (contained in the file `functions/ledm-em.R`).
It requires the following inputs:

- `X`: a list of integer vectors where the nucleotides (A, C, G, T) are coded as (1, 2, 3, 4).
- `p`: length of the templates/references.
- `K`: number of clusters.

The function `sample_ledm` (also contained in the file `functions/ledm-em.R`) allows to simulate from the proposed
data generating process. A minimal example is provided in the following:

```r
p <- 50
K <- 5
data <- sample_ledm(nk = rep(10, K), p = p, theta = c(0.9, 0.03, 0.05, 0.02))
out_ledm <- ledm_em(X = data$X, p = p, K = K)
```