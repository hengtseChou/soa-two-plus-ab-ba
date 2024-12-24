# Optimal alpha-beta and beta-alpha designs for SOA 2+

This repo contains the search result and necessay data and code to reproduce them.

For search results, the recommanded way to run them is to execute these R codes in a terminal with `Rscript`:

```
# 16 runs
cd '16 runs'
Rscript 16runs_data_prep.R
Rscript 16runs_search.R

# 32 runs
cd '32 runs'
Rscript 32runs_data_prep.R
Rscript 32runs_search.R

# 81 runs
cd '81 runs'
Rscript 81runs_data_prep.R
Rscript 81runs_search.R
```

As for `greedy_experiments.R`, one may run it in RStudio to get the resulting plots.
