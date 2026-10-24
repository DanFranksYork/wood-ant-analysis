# Statistical methods
#Author: Daniel Franks

We used Bayesian multilevel models to estimate how genetic differentiation changes with distance within and between colonies, and how distance, genetic differentiation and food-income difference are associated with trail presence and investment. Shared nest and colony effects accounted for dependence among pairs. Partial pooling allows nests, colonies and colony pairs to differ, while drawing less precisely estimated effects towards their population means. We fitted models separately and interpreted posterior predictions as adjusted associations. We considered causal limitations but did not conduct a formal causal or confounding sensitivity analysis.

## Data and notation

Fifteen of 196 nests had no genotype samples. The remaining 181 nests contributed 1,448 ants, with records for 11 diploid microsatellite loci per ant. We combined each marker's two allele calls into one diploid genotype and processed the genotypes using the R package adegenet. We treated `000` as a failed allele call (2,242 of 31,856 allele calls; 7.0%). At each locus, we treated an ant's genotype as missing if either allele call had failed; this affected 1,122 of 15,928 genotypes (7.0%). We excluded one ant with no usable genetic data. Missing genotypes were not imputed; ants contributed data only at loci where they had usable genotypes.

All analyses used the same screened genotypes. Restricting to known colonies retained 847 ants from 106 nests, excluding 600 ants from 75 nests of unknown colony identity. We removed 55 ants missing more than 30% of the 11 loci, then loci with more than 5% missing genotypes among the remaining ants (FE16 and FE21). One nest with only three remaining ants fell below the four-ant minimum and was excluded. These sequential screens retained 789 ants from 105 nests in ten colonies, with nine loci and 6–8 ants per nest; thresholds were not reapplied iteratively.

We calculated pairwise FST from diploid genotypes using Latter’s estimator in the R package hierfstat. We also calculated Weir–Cockerham and Nei estimates to assess agreement. Validation covered allele counts, identifiers, joins, matrix symmetry, unique unordered pairs and invariance to allele order, individual order and exchanging endpoints.

Pair $i=(a_i,b_i)$ joins two distinct nests; $c(j)$ denotes nest $j$'s colony. Within colonies, geographic distance $D_i$ was the Euclidean separation between local map coordinates, converted from centimetres to metres. Between colonies, it was calculated from GPS coordinates on the WGS84 ellipsoid, in metres. We excluded one of 87 internest and two of 97 foraging trail records because their nest identifiers could not be matched. One further internest trail involved the nest excluded by genetic screening, leaving 85. Nest sizes were not imputed; the model-specific exclusions are given below. Internest and foraging surveys were complete for monitored colonies, so unlisted trails represented absences and no listed foraging trail represented zero observed income. Missing recorded strengths remained unknown.

Pair-level exclusions were applied after FST estimation and did not remove genetic data from that calculation. The genetic model retained all 5,460 pairs, comprising 662 within-colony and 4,798 between-colony pairs; none was excluded for missing FST or distance. Both trail models excluded the 4,798 between-colony pairs. Of the remaining 662 pairs, presence excluded 17 because one nest lacked the population estimate needed for food difference, leaving 645. Investment also excluded 577 pairs without trails, then one of 85 observed trails because the same missing estimate prevented calculation of food difference and investment, leaving 84. No additional pairs were excluded for missing or invalid measurements. Exclusions overlap across analyses and should not be summed.

Exclusions for missing genetic data or model predictors are simplest to justify when missingness is unrelated to observed or unobserved characteristics (missing completely at random), which was not established. Selective missingness could bias FST or regression estimates.

Below, Normal distributions are parameterised by mean and standard deviation (SD); HalfNormal denotes the corresponding zero-centred normal restricted to positive values; Exponential uses a rate; and MVNormal uses a covariance matrix. All logarithms are natural. Observations were assumed conditionally independent within each model, given predictors and model parameters. Group-effect vectors are independent across groups and grouping factors conditional on their population parameters. Priors on population-level coefficients and scale parameters are independent unless stated otherwise.

## Analysis of molecular variance (AMOVA)

We partitioned genetic differences among diploid ants into components between colonies, between nests within colonies and between ants within nests, using the screened genotypes and the R package poppr 2.9.8 ([Kamvar et al., 2014](https://doi.org/10.7717/peerj.281)). Within-ant heterozygosity was not a separate component. Allele-frequency distances used available allele entries for each pair, rescaled for the number available, followed by a quasi-Euclidean correction. We tested variance components with 999 hierarchical permutations using ade4 1.7-24 and a fixed random seed. This analysis describes genetic structure rather than directly testing genetic associations with trails.

## Genetic differentiation (Model 1)

The genetic analysis included 5,460 pairs from 105 nests in ten known colonies: 662 within-colony and 4,798 between-colony pairs. Let $\widehat F_i$ denote the Latter estimate and $S_i=1$ for nests in the same colony and $0$ otherwise. All retained estimates lay strictly inside $(0,1)$; none was clipped or shifted. Let $g(i)$ identify the unordered colony combination, including within-colony combinations, and $\bar D_g$ its mean pairwise distance among known-colony pairs before genetic screening. We defined
\[
x_i=\frac{D_i-50}{100},\qquad
 h_i=(1-S_i)\frac{D_i-\bar D_{g(i)}}{100}.
\]

We modelled the logit-transformed estimates with separate within- and between-colony distance slopes and residual SDs:
\[
\begin{aligned}
\operatorname{logit}(\widehat F_i) &\sim \operatorname{Normal}(\mu_i,\sigma_i),\\
\mu_i &= \alpha+\delta S_i+\beta_Bx_i(1-S_i)+\beta_Wx_iS_i\\
      &\quad+u_{a_i}+u_{b_i}+v_{c(a_i)}+v_{c(b_i)}+q_{g(i)}+r_{g(i)}h_i,\\
\log\sigma_i &= \gamma_0+\gamma_S S_i,\\
u_j &\sim \operatorname{Normal}(0,\tau_u),\\
v_k &\sim \operatorname{Normal}(0,\tau_v),\\
q_g &\sim \operatorname{Normal}(0,\tau_q),\\
r_g &\sim \operatorname{Normal}(0,\tau_r),\\
\alpha &\sim \operatorname{Normal}(\operatorname{logit}(0.1),0.8),\\
\delta &\sim \operatorname{Normal}(0,0.7),\\
\beta_B &\sim \operatorname{Normal}(0,0.1),\\
\beta_W &\sim \operatorname{Normal}(0,1),\\
\gamma_0 &\sim \operatorname{Normal}(\log 0.3,0.6),\\
\gamma_S &\sim \operatorname{Normal}(0,0.5),\\
\tau_u,\tau_v,\tau_q &\sim \operatorname{HalfNormal}(0,0.3),\\
\tau_r &\sim \operatorname{HalfNormal}(0,0.5).
\end{aligned}
\]

Here $u_j$ and $v_k$ are multiple-membership random intercepts for nests and colonies. Each nest has one shared effect, whichever endpoint it occupies, and each pair adds its two nest effects. Colony effects follow the same rule, so within-colony pairs receive the same colony effect twice. All 55 colony combinations have baseline deviations $q_g$ beyond these endpoint effects; only the 45 between-colony combinations contribute distance-slope deviations $r_g$. Centring distance places each between-colony combination's baseline deviation at its mean distance before genetic screening. Pair intercept and slope deviations are independent conditional on their SDs.

The population slopes $\beta_W$ and $\beta_B$ describe changes in expected logit FST per 100 m. Each between-colony combination has slope $\beta_B+r_g$, while the within-colony slope is common across colonies. At 50 m with group deviations zero, $\alpha$ is the between-colony mean logit FST and $\delta$ the within-minus-between difference.

The genetic and trail models treat estimated FST values as known; FST estimation uncertainty is not propagated.

## Trail presence and investment

Both trail analyses used within-colony pairs. Let $T_i$ indicate trail presence, $Q_j$ be summed outgoing foraging strength (ants/cm), and $N_j$ be estimated nest population (ants). Food-income difference $H_i$ and positive trail investment $I_i$ were
\[
H_i=\left|\frac{Q_{a_i}}{N_{a_i}}-\frac{Q_{b_i}}{N_{b_i}}\right|,
 \qquad
 I_i=\frac{\lambda_i L_i}{(N_{a_i}+N_{b_i})/2},
\]

where $\lambda_i$ is internest trail strength (ants/cm) and $L_i$ is straight-line separation from local coordinates (cm), so $D_i=L_i/100$ for these within-colony pairs. Thus $H_i$ has units cm$^{-1}$ and investment is dimensionless, without a fixed upper bound of one. Both endpoint sizes had to be finite and positive. No retained investment response used substituted local coordinates.

To retain the prior scales, we centred local distance using the mean GPS distance and divided by its sample SD in the eligible trail-presence dataset before genetic screening. Genetic and food predictors retained their means and SDs from that same reference dataset; Model 1 retained its 50 m reference and between-colony pair centres. Only screened data entered the models. The resulting predictors were $z_{D,i}$ (distance), $z_{F,i}$ (Latter FST) and $z_{H,i}$ (food difference). Investment and robustness fits reused these scaling values without rounding. Presence used 645 pairs (104 nests; 84 trails), and investment used 84 trails (94 nests). All ten colonies were represented in both analyses.

Superscripts $P$ and $I$ distinguish presence and investment. Within each model, $u_j$ is a multiple-membership nest random intercept, summed across the two nests, and $k(i)=c(a_i)=c(b_i)$ is the pair's shared colony. The colony terms $c_{0,k}$ and $c_{D,k}$ are correlated random effects: baseline and log decline magnitude in Model 2, and baseline and distance slope in Model 3. Their covariance matrices are $\Sigma_P$ and $\Sigma_I$, with $2\times2$ correlation matrices $\mathbf R_P$ and $\mathbf R_I$.

### Trail presence (Model 2)

Trail presence followed
\[
\begin{aligned}
T_i &\sim \operatorname{Bernoulli}(p_i),\\
\operatorname{logit}(p_i) &= \alpha_P+b_{D,k(i)}^{(P)}z_{D,i}+b_{F,P}z_{F,i}+b_{H,P}z_{H,i}\\
 &\quad+u_{a_i}^{(P)}+u_{b_i}^{(P)}+c_{0,k(i)}^{(P)},\\
b_{D,k}^{(P)} &= -m_P\exp(c_{D,k}^{(P)}),\\
u_j^{(P)} &\sim \operatorname{Normal}(0,\tau_{u,P}),\\
\begin{bmatrix}c_{0,k}^{(P)}\\c_{D,k}^{(P)}\end{bmatrix}
 &\sim \operatorname{MVNormal}\!\left(\begin{bmatrix}0\\0\end{bmatrix},\Sigma_P\right),\\
\alpha_P &\sim \operatorname{Normal}(-2,1.5),\\
m_P &\sim \operatorname{HalfNormal}(0,3),\\
b_{F,P},b_{H,P} &\sim \operatorname{Normal}(0,1),\\
\tau_{u,P},\tau_{0,P},\tau_{D,P} &\sim \operatorname{Exponential}(2),\\
\Sigma_P &= \operatorname{diag}(\tau_{0,P},\tau_{D,P})\,
             \mathbf R_P\,\operatorname{diag}(\tau_{0,P},\tau_{D,P}),\\
\mathbf R_P &\sim \operatorname{LKJCorr}(2).
\end{aligned}
\]

Here $m_P$ is the shared decline magnitude. The positive multiplier $\exp(c_{D,k}^{(P)})$ lets the decline rate vary among colonies while keeping every colony's distance slope non-positive.

### Trail investment (Model 3)

For observed trails with positive investment, we used
\[
\begin{aligned}
I_i\mid T_i=1 &\sim \operatorname{Lognormal}(\eta_i^{(I)},\sigma_I),\\
\eta_i^{(I)} &= \alpha_I+b_{D,I}z_{D,i}+b_{F,I}z_{F,i}+b_{H,I}z_{H,i}\\
 &\quad+u_{a_i}^{(I)}+u_{b_i}^{(I)}+c_{0,k(i)}^{(I)}
       +c_{D,k(i)}^{(I)}z_{D,i},\\
u_j^{(I)} &\sim \operatorname{Normal}(0,\tau_{u,I}),\\
\begin{bmatrix}c_{0,k}^{(I)}\\c_{D,k}^{(I)}\end{bmatrix}
 &\sim \operatorname{MVNormal}\!\left(\begin{bmatrix}0\\0\end{bmatrix},\Sigma_I\right),\\
\alpha_I &\sim \operatorname{Normal}(\log(0.003),1.2),\\
b_{D,I} &\sim \operatorname{Normal}(0,0.25),\\
b_{F,I} &\sim \operatorname{Normal}(0,0.5),\\
b_{H,I} &\sim \operatorname{Normal}(0,0.15),\\
\tau_{u,I},\tau_{0,I} &\sim \operatorname{HalfNormal}(0,0.5),\\
\tau_{D,I} &\sim \operatorname{HalfNormal}(0,0.3),\\
\sigma_I &\sim \operatorname{HalfNormal}(0,0.75),\\
\Sigma_I &= \operatorname{diag}(\tau_{0,I},\tau_{D,I})\,
             \mathbf R_I\,\operatorname{diag}(\tau_{0,I},\tau_{D,I}),\\
\mathbf R_I &\sim \operatorname{LKJCorr}(2).
\end{aligned}
\]

The lognormal likelihood specifies the mean and SD of $\log I_i$. These associations concern existing trails; selecting them can introduce confounding through unmeasured causes of presence and investment. Investment also contains trail length and shares nest size with food difference, further limiting causal interpretation.

## Prior justification

We followed a Bayesian workflow: specify regularising priors, check prior predictions, fit models, and assess sampling and posterior predictions ([Gelman et al., 2020](https://arxiv.org/abs/2011.01808)). Prior calibration was done using biological expectations, measurement units and the predictor design. Except for the trail-presence distance effect, zero-centred slope and colony-contrast priors allowed either direction while favouring smaller associations. Prior predictive trend draws were plotted to sense-check choice of priors.

Genetic priors favoured modest differentiation: the intercept prior placed 95% of the reference conditional median FST between 0.023 and 0.35. Distance priors allowed comparable changes over 50 m within colonies and 500 m between colonies, reflecting their different spatial scales. Half-normal SD priors regularised nest, colony and colony-pair baselines; the wider pair-slope prior allowed local trends to depart from the overall between-colony trend. The residual prior allowed between-colony SDs of approximately 0.09–0.97 logit units and within/between SD ratios of 0.38–2.66, each with 95% prior probability, permitting substantial dispersion without favouring extreme values.

For trail presence, we assumed that greater distance cannot increase connection probability within a colony, reflecting the costs of maintaining longer connections. The half-normal prior on decline magnitude peaks at zero; scale 3 allowed weak through steep declines in prior predictions. Direction was assumed, so inference concerned the size of the decline. The baseline prior favoured sparse networks; exponential SD priors regularised nest and colony variation.

Investment priors allowed a small but uncertain allocation to one trail relative to mean nest population: the reference conditional median was centred on 0.003, with a 95% prior range of 0.00029–0.032. This is a broad biological judgement, not an empirical estimate. Half-normal nest, colony-intercept and residual SD priors allowed weak and strong trails while regularising extreme variation. Slope priors allowed moderate distance effects, greater uncertainty about genetics (a 95% prior range of 0.38–2.66-fold per predictor SD), and narrower food effects because food spanned more standardised units. LKJ(2) priors mildly favoured correlations near zero while allowing either sign.

## Fitting and model checks

We fitted models with the No-U-Turn Sampler in R 4.6.1 using brms 2.23.0 ([Bürkner, 2017](https://doi.org/10.18637/jss.v080.i01)), cmdstanr 0.9.0 and CmdStan 2.39.0. We ran four chains, initially with 2,000 warmup and 1,000 retained iterations per chain for genetic differentiation, 2,000 and 3,000 for trail presence, and 1,500 and 2,000 for investment. For genetic differentiation, sampling refinement increased the tree-depth limit to 14; the final run used 500 further warmup iterations and 2,000 retained iterations per chain with the learned metric.

We required $\widehat R<1.01$, bulk and tail effective sample sizes above 400, no divergent transitions or tree-depth-limit hits, and energy Bayesian fraction of missing information above 0.3 in every chain. We also checked that the implemented predictors and likelihoods were invariant to exchanging dyad endpoints.

Posterior predictive checks used 400 draws, including fitted group effects and response variation, to assess distributions, distance strata and colonies. Genetic checks also assessed within- and between-colony distributions and residual skewness. These were conditional checks for the sampled nests and colonies. Observed genetic residuals had a longer lower tail than the Gaussian model predicted, limiting predictions for unusually low FST.
