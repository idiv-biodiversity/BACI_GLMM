# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# BACI (Before-After_Control_Impact) design analysis using GLMM frame (Generalized linear mixed model).
# References are provided in comments and the full list is at the end of this script.
# Many thanks to Dr. Dylan Craven for his valuable advice about GLMM.
#
# Note also the repository https://github.com/valentinitnelav/Saldur-ROR-analysis
# There we made use of the more recent package glmmTMB
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# I) Attach packages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# List of packages to be used in the R session
.packages = c("lme4", "AICcmodavg", "MuMIn", "pbkrtest",
              "parallel", "data.table", "blmeco", "lsmeans")
# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])
# Attach packages
sapply(.packages, require, character.only=TRUE)
# Display metadata about the R session at run time
sessionInfo()
# R version 3.4.2 (2017-09-28)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 16.04.3 LTS
# 
# Matrix products: default
# BLAS: /usr/lib/openblas-base/libblas.so.3
# LAPACK: /usr/lib/libopenblasp-r0.2.18.so
# 
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] lsmeans_2.27-2    estimability_1.2  blmeco_1.1        MASS_7.3-47       data.table_1.10.4 pbkrtest_0.4-7   
# [7] MuMIn_1.40.0      AICcmodavg_2.1-1  lme4_1.1-14       Matrix_1.2-11    
# 
# loaded via a namespace (and not attached):
# [1] Rcpp_0.12.13     compiler_3.4.2   nloptr_1.0.4     plyr_1.8.4       tools_3.4.2      nlme_3.1-131     lattice_0.20-35 
# [8] mvtnorm_1.0-6    coda_0.19-1      raster_2.5-8     stats4_3.4.2     grid_3.4.2       reshape_0.8.7    survival_2.41-3 
# [15] VGAM_1.0-4       arm_1.9-3        sp_1.2-5         multcomp_1.4-7   TH.data_1.0-8    minqa_1.2.4      codetools_0.2-15
# [22] splines_3.4.2    abind_1.4-5      unmarked_0.12-2  xtable_1.8-2     sandwich_2.4-0   zoo_1.8-0    


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# II) GLMM fitting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read data
BACI.dt <- fread("Data/BACI.dt.csv", stringsAsFactors = TRUE)
BACI.dt[, Year_F := factor(Year_F)] # Year_F needs to be factor


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# II.a) Specify fixed and random effects and test for random effects structure.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The general model formula with fixed effects is 
# Predation_ratio ~ Period_BA + SiteClass_CI + Period_BA:SiteClass_CI
# Note the interaction term Period_BA:SiteClass_CI which is the "BACI effect",
# Testing for its statistical significance is equivalent to testing for an environmental impact.
# See details in Schwarz (2015).
# To the fixed effects model, random effects are added: SITE and YEAR.
# Below, it was tested if SITE:YEAR interaction should be kept in the model.

Cand.set <- list() # create an empty list to be populated with models
# Fit model without the interaction of random effects (Site and Year)
Cand.set[[1]] <- glmer(CL/NOTAB ~ Period.BA_F*SiteClass.CI_F+
                           (1|Site_F)+(1|Year_F), 
                       data = BACI.dt, weights = NOTAB, family = binomial)
# Model with interaction in the random effects structure
Cand.set[[2]] <- glmer(CL/NOTAB ~ Period.BA_F*SiteClass.CI_F+(1|Site_F)+
                           (1|Year_F)+(1|Site_F:Year_F), 
                       data = BACI.dt, weights = NOTAB, family = binomial)
# AICc comparison
AIC.res.table <- aictab(cand.set = list(Cand.set[[1]], Cand.set[[2]]), 
                        modnames = paste0("Cand.set_", c(1,2)), 
                        second.ord = TRUE)
AIC.res.table
##            K     AICc Delta_AICc AICcWt Cum.Wt       LL
## Cand.set_2 7 13372.71       0.00      1      1 -6679.34
## Cand.set_1 6 14696.73    1324.02      0      1 -7342.35
# Choose random effects structure with lowest AICc.
# Conclusion: 
# AICc suggests to include the term (1|Site_F:Year_F) in the final model,
# that is, the interaction of random effects.


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# II.b) Testing for overdispersion in the binomial mixed model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# "Overdispersion: the occurrence of more variance in the data than predicted by a statistical model." (Bolker, 2009)
# In case of overdispersed then fit model with observation-level random effects 
# see Harrison XA (2014) or Harrison XA (2015)

# Full model without overdispersion control
model.full.no.disp.ctrl <- Cand.set[[2]]
# Full model with overdispersion control
model.full.disp.ctrl <- update(Cand.set[[2]],.~.+(1|obs))

# Measure overdispersion in the two binomial glmer-models
blmeco::dispersion_glmer(model.full.no.disp.ctrl)
## [1] 1.390638 
# If the value is between 0.75 and 1.4, there may not be an overdispersion
blmeco::dispersion_glmer(model.full.disp.ctrl)
## [1] 1.010785
# The scale parameter dropped from 1.4 to 1.0;
# accounting for overdispersion can be justified.


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# II.c) Check model assumptions (visually)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
qqnorm(residuals(model.full.disp.ctrl), 
       main = "Q-Q plot - residuals")
qqline(residuals(model.full.disp.ctrl), col="red")

# inspecting the random effects (see also Bolker, 2009 - supp 1)
qqnorm(unlist(ranef(model.full.disp.ctrl)), 
       main = "Q-Q plot, random effects")
qqline(unlist(ranef(model.full.disp.ctrl)), col="red")

# fitted vs residuals
scatter.smooth(fitted(model.full.disp.ctrl), 
               residuals(model.full.disp.ctrl, type="pearson"),
               main="fitted vs residuals",
               xlab="Fitted Values", ylab="Residuals")
abline(h=0, col="red")

# fitted vs observed
scatter.smooth(fitted(model.full.disp.ctrl), BACI.dt$PropCL,
               xlab="Fitted Values", ylab="Observed PropCL")
abline(0,1, col="red")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# II.d) Test for statistical significance of BACI interaction term
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Using a parametric bootstrap comparison between nested models
# (see Halekoh & Højsgaard 2014 and Bolker, 2009 - supp 1)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# II.d.1) For model with all data points
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model without the BACI interaction term
model.no.interaction <- glmer(CL/NOTAB ~ Period.BA_F+SiteClass.CI_F+
                                  (1|Site_F)+(1|Year_F)+(1|Site_F:Year_F)+(1|obs), 
                              data = BACI.dt, weights = NOTAB, family = binomial)

# Calculate reference distribution of likelihood ratio statistic
# Warning: is computationally slow (can take 3-4 hours for 100 simulations)
refdist.pb.100.interaction <- PBrefdist(largeModel = model.full.disp.ctrl, 
                                        smallModel = model.no.interaction, 
                                        nsim = 100, seed = 2017)
# Model comparison test using the reference distribution from above
compar.interaction.100 <- PBmodcomp(largeModel = model.full.disp.ctrl, 
                                    smallModel = model.no.interaction,
                                    ref = refdist.pb.100.interaction)
compar.interaction.100
##           stat df   p.value 
## LRT    34.832  1 3.595e-09 *** # likelihood ratio test statistic
## PBtest 34.832      0.01389 *   # parametric bootstrap test (is more conservative)
# Both likelihood ratio and parametric bootstrap tests give significant p-values.
# Therefore, the interaction term (the BACI effect) is statistically significant
# (there is an environmental impact)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# II.d.2) For model only with single point ("Snapshot") method of estimating seed predation 
# This is to indicate that method of estimating pre-dispersal seed predation 
# does not influence the statistical significance of the BACI effect proven above.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
BACI.dt2 <- BACI.dt[Method_F == "Snapshot"] # subset
# Fit model for single point ("Snapshot") case - model without overdispersion control
model.snap <- glmer(CL/NOTAB ~ Period.BA_F*SiteClass.CI_F+
                        (1|Site_F)+(1|Year_F)+(1|Site_F:Year_F), 
                    data = BACI.dt2, weights = NOTAB, family = binomial)
# Account for overdispersion
model.snap.disp.ctrl <- update(model.snap,.~.+(1|obs))
blmeco::dispersion_glmer(model.snap)
## [1] 1.393408
blmeco::dispersion_glmer(model.snap.disp.ctrl)
## [1] 1.008827
# Here also, the scale parameter dropped from 1.4 to 1.0; 
# accounting for overdispersion can be justified.

# Model without the BACI interaction term
model.no.interaction.snap <- glmer(CL/NOTAB ~ Period.BA_F+SiteClass.CI_F+
                                       (1|Site_F)+(1|Year_F)+(1|Site_F:Year_F)+(1|obs), 
                                   data = BACI.dt2, weights = NOTAB, family = binomial)
# Calculate reference distribution of likelihood ratio statistic
# Warning: is computationally slow (can take 2-3 hours for 100 simulations)
refdist.pb.100.interaction.snap <- PBrefdist(largeModel = model.snap.disp.ctrl, 
                                             smallModel = model.no.interaction.snap, 
                                             nsim = 100, seed = 2017)
# Model comparison test using the reference distribution from above
compar.interaction.100.snap <- PBmodcomp(largeModel = model.snap.disp.ctrl, 
                                         smallModel = model.no.interaction.snap,
                                         ref = refdist.pb.100.interaction.snap)
compar.interaction.100.snap
##          stat df   p.value    
## LRT    34.903  1 3.466e-09 ***
## PBtest 34.903      0.01333 * 
# Also when using only observations corresponding to single point ("Snapshot") method of estimating seed predation,
# both likelihood ratio and parametric bootstrap tests give significant p-values for BACI effect.


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# III) Extract coefficients from final model 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
final.model <- model.full.disp.ctrl
final.model.noIntercept <- update(model.full.disp.ctrl, . ~ . -1)
# Remove the intercept for extracting group means and their confidence intervals.
# see Schielzeth H (2010)

# Note that the values from summary() are on original (logit) scale.
# To get the real proportions (Least-squares means/predicted marginal means/treatment means) one can do:
estimates <- lsmeans::lsmeans(final.model.noIntercept, ~ SiteClass.CI_F:Period.BA_F, type = "response")
# https://stats.stackexchange.com/questions/192062/issue-calculating-adjusted-means-for-glmer-model
# Confidence level used: 0.95 
estimates
# SiteClass.CI_F Period.BA_F       prob         SE df   asymp.LCL asymp.UCL
## Control        After       0.52995542 0.09008968 NA 0.356894394 0.6961011
## Impact         After       0.02982617 0.01797237 NA 0.009018536 0.0940835
## Control        Before      0.73204824 0.08241147 NA 0.545271980 0.8615822
## Impact         Before      0.77994713 0.11021259 NA 0.501690025 0.9258043

# As indicated in Schwarz CJ (2015), the BACI effect is computed as 
# BACI = avgCA-avgCB-(avgIA-avgIB)
est <- predict(ref.grid(final.model.noIntercept), type = "response") # get the estimates only (without CI)
names(est) <- c("CA","CB","IA","IB"); est # give names to the estimates vector
baci <- est["CA"]-est["CB"]-(est["IA"]-est["IB"])
baci # 0.5480281
# One can also get the BACI effect like 
contrast(regrid(estimates), list(baci=c(1,-1,-1,1)))
## contrast  estimate        SE df z.ratio p.value
## baci     0.5480281 0.1236863 NA   4.431  <.0001
## Confidence level used: 0.95
# Or with asymptotic CI-s
confint(contrast(regrid(estimates), list(baci=c(1,-1,-1,1)))) 
## contrast  estimate        SE df asymp.LCL asymp.UCL
## baci     0.5480281 0.1236863 NA 0.3056074 0.7904489
## Confidence level used: 0.95 
# https://stats.stackexchange.com/questions/241523/testing-for-pairwise-proportion-differences

# Calculate conditional and marginal coefficient of determination
MuMIn::r.squaredGLMM(final.model)  
## R2m       R2c 
## 0.3305502 0.4736158
# Rm2 = represents the variance explained by fixed factors (Marginal R_GLMM²)
# R2c = variance explained by both fixed and random factors, i.e. the entire model (Conditional R_GLMM²)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# REFERENCES
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Bolker BM et al. (2009) Generalized linear mixed models: a practical guide for ecology and evolution.
#   Trends in ecology & evolution 24:127–135
#   at http://www.sciencedirect.com/science/article/pii/S0169534709000196

# Halekoh U, Højsgaard S (2014) A kenward-roger approximation and parametric bootstrap methods for
#   tests in linear mixed models–the R package pbkrtest. Journal of Statistical Software 59:1–32
#   at http://www.jstatsoft.org/v59/i09/

# Harrison XA (2014) Using observation-level random effects to model overdispersion 
#   in count data in ecology and evolution. PeerJ 2:e616 https://doi.org/10.7717/peerj.616

# Harrison XA (2015) A comparison of observation-level random effect and Beta-Binomial models for
#   modelling overdispersion in Binomial data in ecology & evolution. PeerJ 3:e1114
#   at https://peerj.com/articles/1114.pdf

# Schwarz CJ (2015) Analysis of BACI experiments. In Course Notes for Beginning and Intermediate Statistics.
#   at http://people.stat.sfu.ca/~cschwarz/Stat-650/Notes/PDFbigbook-R/R-part013.pdf
#   Note that this is not maintained by the stat.sfu.ca anymore. Use https://web.archive.org/ to retrieve it. For example
#   https://web.archive.org/web/20190107174420/http://people.stat.sfu.ca/~cschwarz/Stat-650/Notes/PDFbigbook-R/R-part013.pdf

# Schielzeth H (2010) Simple means to improve the interpretability of regression coefficients. 
#   Methods in Ecology and Evolution, 1(2), 103-113.
#   at https://besjournals.onlinelibrary.wiley.com/doi/pdfdirect/10.1111/j.2041-210X.2010.00012.x
