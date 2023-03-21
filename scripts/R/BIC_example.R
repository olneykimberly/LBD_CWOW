# devtools::install_github("GabrielHoffman/mvIC", repos=BiocManager::repositories())
library(mvIC)
#---------------------------
# Multivariate linear model
# Predict Sepal width and Length given Species
# Evaluate model fit
fit1 = lm( cbind(Sepal.Width, Sepal.Length) ~ Species, data=iris)
mvIC( fit1 )
# add Petal width and length
# smaller mvIC means better model
fit2 = lm( cbind(Sepal.Width, Sepal.Length) ~ Petal.Width + Petal.Length + Species, data=iris)
mvIC( fit2 )

#---------------------------
# Multivariate linear mixed model
library(lme4)
# Predict Sepal width and Length given Species
# Evaluate model fits separately
fit1 = lmer( Sepal.Width ~ (1|Species), data=iris, REML=FALSE)
fit2 = lmer( Sepal.Length ~ (1|Species), data=iris, REML=FALSE)
score = mvIC( list(fit1,fit2) )

# Use mvIC_fit() to fit multivariate regression
# Combine responses on *rows*
Y = with(iris, rbind(Sepal.Width, Sepal.Length))

# Evaluate multiple responses
# Can specify any formula of fixed or random effects
mvIC_fit( Y, ~ (1|Species), data=iris)

#---------------------------
# Forward stepwise regression
# Combine responses on *rows*
Y = with(iris, rbind(Sepal.Width, Sepal.Length))

# variables to consider in the model
variables = c("Petal.Length", "Petal.Width", "(1|Species)")

# fit forward stepwise regression starting with model: ~1. 
bestModel = mvForwardStepwise( Y, ~ 1, data=iris, variables=variables)

#-----
# Forward stepwise regression
# Combine responses on *rows*
Y = with(dge.filtered.norm, rbind(dge.filtered.norm$samples$sex_inferred, 
                                  dge.filtered.norm$samples$RIN, 
                                  dge.filtered.norm$samples$Age, dge.filtered.norm$samples$APOE_E4_allele_count,
                                  dge.filtered.norm$samples$PCT_CODING_BASES,
                                  dge.filtered.norm$samples$PCT_INTERGENIC_BASES, 
                                  dge.filtered.norm$samples$PCT_INTRONIC_BASES))

# variables to consider in the model
variables = c("(1|sex_inferred)", "RIN", "Age", "APOE_E4_allele_count", "PCT_CODING_BASES", "PCT_INTERGENIC_BASES", "PCT_INTRONIC_BASES", "(1|TYPE)")
# fit forward stepwise regression starting with model: ~1. 
bestModel = mvForwardStepwise( Y, ~ 1, data=dge.filtered.norm$samples, variables=variables)

#-----------------------------------------------
baseFormula_2 <- ~ (1 | TYPE) + PCT_CODING_BASES + RIN + (1 | flowcell_and_lane) + PCT_INTRONIC_BASES + PCT_INTERGENIC_BASES
# Combine responses on *rows*
Y = with(
  info,
  rbind(
    sex_inferred,
    APOE,
    scaled.info.df$APOE_E4_allele_count
  )
)

rownames(Y) <-
  c(
    "sex_inferred",
    "APOE",
    "APOE_E4_allele_count"
  )
# variables to consider in the model
# categorical variables must be modeled using (1|)
variables = c(
  "(1|sex_inferred)",
  "(1|APOE)",
  "APOE_E4_allele_count",
  "(1 | TYPE:sex_inferred)",
)

# fit forward stepwise regression starting
bestModel_voomcounts_2 = mvForwardStepwise(voomCounts,
                                         baseFormula_2,
                                         data = info,
                                         variables = variables)
bestModel_voomcounts_2



