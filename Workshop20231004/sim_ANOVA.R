#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#			Simulation for ANOVA
#		
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Author: Zach Siders
#Incept: 10/3/2023

###-----------------------------------------------------
#		Initialization
###-----------------------------------------------------
	library(glmmTMB)
	library(emmeans)
###-----------------------------------------------------
#		Data Generation (One-way ANOVA)
###-----------------------------------------------------
	# Kéry (2010) 9.2.1 pg. 119
	npops <- 5				# Number of populations
	nsample <- 10				# Number of snakes in each
	pop_means <- c(50, 40, 45, 55, 60) 	# Population mean SVL
	sigma <- 3				# Residual sd

	#Generate covariate (categorical for ANOVA)
	n <- npops * nsample 			# Total number of data points
	pop <- gl(n = npops, k = nsample, length = n)# Indicator for population
	means <- rep(pop_means, rep(nsample, npops))

	# Create design matrix
	X <- as.matrix(model.matrix(~ pop-1)) 
	head(X)

	#individual mean given effects
	ind_means <- X %*% as.matrix(pop_means)
	#random draw of y (DON'T FORGET THE ERROR!)
	y <- rnorm(n, 
			   mean = ind_means, 
			   sd = sigma)
	#alternative generation of y using matrix math
	eps <- rnorm(n, 
				mean = 0, 
				sd = sigma) # Residuals 

	y_alt <- ind_means + eps #additive error structure (normal)

	df1way <- data.frame(pop = pop,
						 len = y)
###-----------------------------------------------------
#		Data Viz (One-way ANOVA)
###-----------------------------------------------------
	boxplot(len ~ pop, data = df1way, 
			col="grey", las = 1, 
			xlab="Population", ylab="SVL", main="")
###-----------------------------------------------------
#		Model (One-way ANOVA)
###-----------------------------------------------------
	aov1way <- glmmTMB(len ~ pop,
					   data = df1way)
	summary(aov1way) # Model summary (p-value is Wald's test)
	summary(aov(aov1way)) #F-test for population effect
	#grab fixed effects
	summary(aov1way)$coef$cond #Fixed effect table
	#grab sigma
	summary(aov1way)$sigma
###-----------------------------------------------------
#		Expected Marginal Means (One-way ANOVA)
###-----------------------------------------------------
	(em1way <- emmeans(aov1way, ~pop)) #easy way to get back pop-level effect
	plot(em1way)

	(em1way_pairs <- pairs(em1way)) #pairwise comparison between populations
###-----------------------------------------------------
#		One-way ANOVA, multiplicative error
###-----------------------------------------------------
	sigma_me <- 0.1
	y_me <- rlnorm(n, 
				   mean = log(ind_means),
				   sd = sigma_me)
	df1way$len_me <- y_me
	
	boxplot(len_me ~ pop, data = df1way, 
			col="grey", las = 1, 
			xlab="Population", ylab="SVL", main="")

	aov1way_me_bad <- glmmTMB(len_me ~ pop,
							data = df1way,
							family = gaussian(link='log'))

	summary(aov1way_me_bad) # Model summary (p-value is Wald's test)
	summary(aov(aov1way_me_bad)) #F-test for population effect
	summary(aov1way_me_bad)$sigma
	(em1way_me_bad <- emmeans(aov1way_me_bad, ~pop,
							type = 'response')) #easy way to get back pop-level effect
	(em1way_pairs_me_bad <- pairs(em1way_me_bad)) #pairwise comparison between population

	aov1way_me <- glmmTMB(log(len_me) ~ pop,
							data = df1way)

	

	summary(aov1way_me) # Model summary (p-value is Wald's test)
	summary(aov(aov1way_me)) #F-test for population effect
	summary(aov1way_me)$sigma

	(em1way_me <- emmeans(aov1way_me, ~pop, 
							type='response')) #easy way to get back pop-level effect
	(em1way_pairs_me_ratio <- pairs(em1way_me)) #pairwise comparison between population
	(em1way_pairs_me <- pairs(em1way_me, 
								ratios=FALSE)) #pairwise comparison between population

	#compare bad and good
	cbind(bad = summary(aov1way_me_bad)$coef$cond[,1],
	      good = summary(aov1way_me)$coef$cond[,1])
###-----------------------------------------------------
#		Data Generation (Two-way ANOVA)
###-----------------------------------------------------
	#see Kéry (2010) 10.2 pg. 131 for another example
	# we are going to riff on this but keep our snake example
	nelev <- 3
	n <- npops * nelev * nsample 			# Total number of data points
	pop_2way <- gl(n = npops, k = nsample, length = n)
	elev <- gl(n = nelev, k = nsample / nelev, length = n) #indicator for elevation group
	elev_means <- c(5, 0, -3) # adjustmean elevation has on pop_mean
	pop_means_adj <- pop_means #need to make pop_means have intercept
	pop_means_adj[-1] <- pop_means[-1] - pop_means[1]
	int_effect <- rnorm((npops-1) * (nelev-1),
	                    mean = 0,
	                    sd = 1) #std normal draw for interaction effects
	all_effects <- c(pop_means_adj, 
	                 elev_means[-1],
	                 int_effect)
	all_effects[1] <- all_effects[1] + elev_means[1] #adjust intercept
	all_effects #check them out
	X_2way <- as.matrix(model.matrix(~pop_2way * elev))
 	ind_means_2way <- X_2way %*% as.matrix(all_effects)
	head(ind_means_2way,50)
	#visualize the matrix math
	head(X_2way %*% diag(all_effects), 50)

	#random draw of y (DON'T FORGET THE ERROR!)
	y_pe <- rnorm(n, 
	              mean = ind_means_2way, 
	              sd = sigma)
	df2way <- data.frame(pop = pop,
	                     elev = elev,
	                     len = y_pe)
###-----------------------------------------------------
#		Data Viz (Two-way ANOVA)
###-----------------------------------------------------
	boxplot(len ~ pop+elev, data = df2way, 
			col=rep(c('gray30','gray70','gray90'), each=npops), las = 1, 
			xlab="Population*Elevation", ylab="SVL", main="")
	legend('topleft', legend=unique(elev),
	       fill = c('gray30','gray70','gray90'))
###-----------------------------------------------------
#		Model (Two-way ANOVA)
###-----------------------------------------------------
	aov2way <- glmmTMB(len ~ pop*elev,
					   data = df2way)
	summary(aov2way) # Model summary (p-value is Wald's test)
	summary(aov(aov2way)) #F-test for population effect
	#grab fixed effects
	summary(aov2way)$coef$cond #Fixed effect table
	#grab sigma
	summary(aov2way)$sigma
###-----------------------------------------------------
#		Expected Marginal Means (Two-way ANOVA)
###-----------------------------------------------------
	(em2way_pop <- emmeans(aov2way, ~pop)) #easy way to get back pop-level effect
	(em2way_elev <- emmeans(aov2way, ~elev)) #easy way to get back elev-level effect
	(em2way <- emmeans(aov2way, ~pop+elev)) #easy way to get back pop*elev-level effect
	plot(em2way)

	(em2way_pop_pairs <- pairs(em2way_pop)) #pairwise comparison between populations
	(em2way_elev_pairs <- pairs(em2way_elev)) #pairwise comparison between elevations
	(em2way_pairs <- pairs(em2way)) #pairwise comparison between both


















	