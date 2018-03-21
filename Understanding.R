require(glmnet)
require(pROC)


G = read.csv('~/Downloads/understand.csv')
G_scaled = G
ind = sapply(G_scaled, is.numeric)
G_scaled[ind] = lapply(G_scaled[ind], scale)

# Note:
# The coefficients in the models were interpreted as follows
# Since the models are fit to scaled data
# Of the form (x - mu) / sd
# Multiplicative change associated with a one-unit increase of x:
#
# e^(B / sd)
#
# Because e^(B * (x + 1 - mu) / sd)
# = e^(B / sd) * e^(B * (x - mu) / sd)
#    factor          value

####

hist(G$AU)
hist(G$PE.spec..java.)

#####

# Some descriptive statistics for poster/paper
tab = table(G$AU > 0.6, G$PE.gen)
tab
ptab= prop.table(tab, margin=2)
ptab

barplot(ptab)

table(subset(G, PE.gen == 10)$AU)

boxplot(G$AU ~ G$JavaProfessional)

hist(subset(G$AU, G$JavaProfessional == TRUE), breaks="FD")

hist(subset(G$AU, G$JavaProfessional == FALSE), breaks="FD")

library(ggplot2)

levels(G$JavaProfessional) = c("No", "Yes")
ggplot(data=G, aes(x=G$AU, fill = factor(G$JavaProfessional, labels=list("No", "Yes")))) + geom_density(alpha=0.5) +
  xlab("% Understood") + ylab("Density") + theme(axis.title = element_text(size = 26)) +
  scale_fill_discrete(name=">5yrs Java Exp.") +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.text=element_text(size=16), legend.title=element_text(size=16)) +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    #plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  theme(axis.text=element_text(size=12),
        axis.title.y=element_text(size=14,),
        axis.title.x = element_text(size=14),
        legend.title=element_text(size=14),
        plot.title = element_text(size=16,hjust=1))
ggsave('UnderstandabilityPlot.png',width=8,height=5,units="in",bg = "transparent",dpi=300)

# PCA

ignore = colnames(G)[colSums(is.na(G)) > 0]
# Remove factors / irrelevant things
ignore = c(ignore, c("participant_id", "system_name", "snippet_signature", "developer_position", "Professional", "JavaProfessional", "Understood"))

pca = prcomp(G[, !(names(G) %in% ignore)], scale = TRUE)

par(mar=c(1, 5, 1, 1))
plot(pca, main="", cex.lab=1.2, cex.axis=1.2)
# Scree plot in the paper
screeplot(pca)

a1 = pca$rotation[,1]
sort(abs(a1))

sort(abs(pca$rotation[,1])) # Suggests Cyclomatic complexity or X.conditionals..dft.
sort(abs(pca$rotation[,2])) # Suggests Volume
sort(abs(pca$rotation[,3])) # Suggests X..words.max. or LOC
sort(abs(pca$rotation[,4])) # Suggests X.operators..avg. or String.identifiers..area.
sort(abs(pca$rotation[,5])) # Operators..Visual.X. or Line.length..max.
sort(abs(pca$rotation[,6])) # X.spaces..avg or Line.length..avg.

# New computed features
G_scaled$Understood = G$Understood = factor(G$AU > 0.6)
G_scaled$Professional = G$Professional = factor(G$PE.gen >= 5)
G_scaled$JavaProfessional = G$JavaProfessional = factor(G$PE.spec..java. >= 5)

# The model inspired by PCA
require(lme4)
lg.fit = glmer(Understood ~ 
                 Professional +
                 Cyclomatic.complexity + 
                 LOC + 
                 Volume + 
                 X.operators..avg. + 
                 Line.length..max. +
                 X.spaces..avg. + (1 | participant_id), data=G_scaled,  # Fitting on the scaled data
               glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e5)), family = binomial(link="logit"))

summary(lg.fit)
require(MuMIn)
r.squaredGLMM(lg.fit)

######## Stepwise
# Forward stepwise

terms = c("(1|participant_id)")
bestTerms = c("(1|participant_id)")

exclude = c("participant_id", "system_name", "snippet_signature", "developer_position", "AU", "Understood", "PBU", "TAU", "TNPU", "Identifiers.comments..area.", "Keywords.comments..area.",    "Numbers.comments..area.",     "Strings.comments..area.",    
            "Literals.comments..area.",    "Operators.comments..area.",   "Strings.numbers..area.",      "Literals.numbers..area.",    
            "Operators.numbers..area.",    "Literals.strings..area.", "Operators.strings..area.", "Operators.literals..area..1")
considering = setdiff(colnames(G), exclude)

mer.formula = reformulate(terms, response="Understood")
bestAIC = 999999
bestChoice = ""

for(j in 1:10) {
  for(i in 1:length(considering)) {
    choice = considering[i];
    oldTerms = terms
    
    print(choice)
    
    terms = c(terms, choice)
    formula = reformulate(terms, response="Understood")
    
    print(formula)
    
    mer.fit = glmer(formula, data=G_scaled, family = binomial(link="logit"))
    newAIC = AIC(mer.fit)
    
    print(newAIC)
    
    if(newAIC < bestAIC) {
      bestAIC = newAIC
      bestChoice = choice;
    }
    
    terms = oldTerms
  }
  
  bestTerms = c(bestChoice, bestTerms)
  terms = c(terms, bestChoice)
  print("The best formula so far is:")
  print(bestTerms)
}

bestFormula = reformulate(bestTerms, response="Understood")

mer.best.fit = glmer(bestFormula, data=G_scaled, family=binomial(link="logit"), glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e5)))

summary(mer.best.fit)
r.squaredGLMM(mer.best.fit)
vif.mer(mer.best.fit)

infl = influence(mer.best.fit, obs=TRUE, count=TRUE)
thresh = 4 / 324
cooks = cooks.distance(infl, sort=TRUE)
removed = c(27, 190, 275, 124, 219, 161, 199, 207)
G[removed, ]$TAU # Tend to have a TAU of 0
G[removed, ]$JavaProfessional # Removed a comparatively large number of Java Professionals
G[removed, ]$Professional

# Did not remove all 16 points over the threshold
# selected eight most influential
mer.best.fit.fixed = glmer(bestFormula, data=G_scaled[-removed, ], family=binomial(link="logit"), glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=1e5)))
summary(mer.best.fit.fixed)
r.squaredGLMM(mer.best.fit.fixed)

# Bootstrapping
require(boot)

b_par = bootMer(x=mer.best.fit, FUN=fixef, nsim=2000, parallel=c("multicore"), ncpus=8)
b_par

format.perc <- function(probs, digits) {
  paste(format(100 * probs, trim = TRUE,
               scientific = FALSE, digits = digits),
        "%")
}

bCI.tab <- function(b,ind=length(b$t0), type="perc", conf=0.95) {
  btab0 <- t(sapply(as.list(seq(ind)),
                    function(i)
                      boot.ci(b,index=i,conf=conf, type=type)$percent))
  btab <- btab0[,4:5]
  rownames(btab) <- names(b$t0)
  a <- (1 - conf)/2
  a <- c(a, 1 - a)
  pct <- format.perc(a, 3)
  colnames(btab) <- pct
  return(btab)
}


# all of the coefs such that p < 0.05 have CIs that do not include 0
boot.ci(b_par,type="perc",index=1)
bCI.tab(b_par)

library(pbkrtest)
library(parallel)

(nc <- detectCores())

cl <- makeCluster(rep("localhost", nc))

mer_new <-update(mer.best.fit, . ~ (1|participant_id))
PBmodcomp(mer.best.fit,mer_new,nsim=200,cl=cl)

stopCluster(cl)

################################

# Ignore these in the LASSO regression
ignore2 = colnames(G)[colSums(is.na(G)) > 0]
ignore2 = c(ignore2, c("participant_id", "system_name", "snippet_signature", "developer_position", "AU", "Understood", "JavaProfessional", "Professional", "PBU", "TAU", "TNPU", "PE.gen", "PE.spec..java."))

# Will say that they understood it if they passed (got a 0.666 or more)
y = factor(G$AU > 0.6)
x = scale(G[, !(names(G) %in% ignore2)])

# Dummy variables for glmnet
d_prof = as.matrix(as.numeric(G$PE.gen >= 5))
d_java = as.matrix(as.numeric(G$PE.spec..java. >= 5))
colnames(d_prof) = c("Professional")
colnames(d_java) = c("JavaProfessional")

x = cbind(x, d_prof)
x = cbind(x, d_java)


# Dummy variables required by glmnet
#x = cbind(x, model.matrix(~ Professional - 1, data=G))
#x = cbind(x, model.matrix(~ JavaProfessional - 1, data=G))

colnames(x)

# Random CV
res = NA
aucs = c()

kept = c()

# Change to 5000, but 1000 should be sufficient for replication of our result
for(i in 1:5000) {
  train = sample(1:nrow(G), nrow(G) * 0.75)
  test = (-train)
  
  grid = 10^seq(10,-2,length=100)
  
  cv.out = cv.glmnet(as.matrix(x[train,]), y[train], alpha=1, nfold = 10, family = "binomial")
  bestlam = cv.out$lambda.min
  
  out = glmnet(as.matrix(x)[train,], y[train], alpha=1, lambda=grid, family="binomial")
  lasso.pred = predict(out, s=bestlam, newx=x[test,], type="response")
  coef = as.matrix(predict(out, type="coefficients", s=bestlam))
  names_c = names(subset(coef, coef[,1] != 0)[,1])
  kept = c(kept, names_c)
  
  roc.test = roc(y[test] ~ lasso.pred, auc=T)
  aucs = c(aucs, pROC::auc(roc.test))
  print(i)
  
  ret = matrix(sapply(seq(0, 1, 0.025), function(x) coords(roc.test, x, input="specificity", output="sensitivity"))["sensitivity",], nrow=1, ncol=41)
  if(is.matrix(res)) {
    res = rbind(res, ret)
  } else {
    res = ret
  }
}


sps = seq(0, 1, 0.025)
quantiles = apply(res, 2, function(x) quantile(x, probs=c(0.025, 0.975)))
means = apply(res, 2, function(x) mean(x))

lower = quantiles[1,]
upper = quantiles[2,]

points(sps, lower)
points(sps, upper)
points(sps, means)

require(ggplot2)
require(plotROC)
library(RColorBrewer)

# Plot the ROC curve

cbbPalette <- brewer.pal(9, "Set1")

df = data.frame(sps=sps, sens=means, lower=lower, upper=upper)
ggplot(df, aes(x=(1-sps), y=sens)) + geom_ribbon(aes(ymin=lower, ymax=upper), fill="grey90") +
  theme_bw() +
  theme(legend.position = c(0.4, 0.2), 
        legend.direction="horizontal",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  ggtitle("Understood Average ROC") +
  geom_line(color="steelblue", size=1.2) +
  geom_abline(slope=1, intercept=0, size=1) +
  labs(y="True Positive Rate", x="False Positive Rate") + 
  theme(legend.position = "none") +
  theme_bw() +
  theme(legend.position = c(0.7, 0.9),
        legend.direction="horizontal",
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        plot.title = element_text(size=16,hjust=1))

# The variables retained in the model can be found in the variable `kept`
# We split this at "(Intercept)" to find the percentage occurrence
# of each parameter in the model





