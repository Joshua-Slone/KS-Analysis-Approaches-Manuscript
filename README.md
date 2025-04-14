1. bootstrap_MI_for_cluster.R contains the bootstrap code and runs it before calling the multiple imputation (MI) functions in MI_region.R and MI_site.R
2. MI_region.R and MI_site.R are MI functions that differ only in their imputation of the KS indicator. MI_region.R imputes KS by region, and MI_site.R imputes KS by clinic site.
3. GR_IPW_EHR_prevalence.Rmd, GR_IPW_EHR_incidence.Rmd, and GR_IPW_EHR_mortality.Rmd contain the code for the generalized raking (GR), inverse probability weighted (IPW), and EHR estimators for each outcome, respectively.
