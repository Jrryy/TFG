library(parallel)

if (!'CORElearn' %in% installed.packages()){
	install.packages('CORElearn')
}

library(CORElearn)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!'mixOmics' %in% installed.packages()){
	BiocManager::install("mixOmics", ask=FALSE)
}

library(mixOmics)

####### FISHER #######

# to calculate Fisher score of each feature: 
fisher = function(gene, pos, neg) {
	score = (mean(gene[pos]) - mean(gene[neg]))**2 / (length(pos)*var(gene[pos]) + length(neg)*var(gene[neg]))
    	return(score)
}

apply_fisher = function(dataset, positive, negative, features_to_keep = 5000, debug = TRUE){

  output = NULL
  
	scores = apply(dataset, 2, fisher, positive, negative)

	print(length(scores))

	stopifnot(length(scores) == dim(dataset)[2])

	sorted_scores = sort(scores, decreasing=TRUE, index.return=TRUE)
	if (debug){
	  plot(sorted_scores$x)
	  
	  dev.new()
	} else {
	  output$to_plot = sorted_scores$x
	}
	# We keep the best x features according to fda, 5000 by default but the user should be able to input it.
	output$data = dataset[, sort(sorted_scores$ix[0:features_to_keep])]

	return(output)
}


####### RELIEFF #######

set.seed(1337)

# Applies ReliefF. For now, only ReliefFequalK, but I plan on using exp too and adding a parameter to specify it
apply_relieff = function(dataset, classes, features_to_keep = 500, iterations = 0, estimator = 'ReliefFexpRank', debug = TRUE){
	
	stopifnot(length(classes) == dim(dataset)[1])

	df_ = as.data.frame(dataset)
	df = cbind(df_, class=classes)

	reliefF_attrs = attrEval('class', df, estimator=estimator, ReliefIterations=iterations)

	sorted_attrs = sort(reliefF_attrs, decreasing=TRUE, index.return=TRUE)

	plot(sorted_attrs$x)
	dev.new()
	
	after_relieff = dataset[, sort(sorted_attrs$ix[0:features_to_keep])]

	heatmap(after_relieff)
	dev.new()

	return(after_relieff)
}

####### PCA #######

apply_pca = function(dataset, classes, components = 10, debug = TRUE){
	pca_plot = pca(dataset, ncomp=components)
	plot(pca_plot)
	dev.new()

	plotIndiv(pca_plot, ind.names=FALSE, group=classes, legend=TRUE, ellipse=TRUE)
	dev.new()
}

####### PLS-DA #######

apply_plsda = function(dataset, classes, components = 10, cv_folds = 5, cv_repeats = 10, debug = TRUE){
	plsda_data = plsda(dataset, classes, ncomp=components, scale=FALSE)

  perf_plsda_data = perf(plsda_data, validation = 'Mfold', folds = cv_folds, progressBar = TRUE, nrepeat = cv_repeats)
  matplot(perf_plsda_data$error.rate$BER, type = 'l', lty = 1, col = color.mixo(1:3), main = 'Balanced Error Rate')
  legend('topright', c('max.dist', 'centroids.dist', 'mahalanobis.dist'), lty = 1, col = color.mixo(1:3))
  dev.new()
  
  list_keepX = c(5:10, seq(20, 100, 10))
  tune_splsda_data = tune.splsda(dataset, classes, ncomp = components, validation = 'Mfold', folds = cv_folds, dist = 'max.dist', progressBar = TRUE, measure = 'BER', test.keepX = list_keepX, nrepeat = cv_repeats)
  
  error = tune_splsda_data$error.rate
  final_ncomp = tune_splsda_data$choice.ncomp$ncomp
  select_keepX = tune_splsda_data$choice.keepX[1:final_ncomp]
  plot(tune_splsda_data, col = color.jet(components))
  dev.new()
  
  final_splsda = splsda(dataset, classes, ncomp = final_ncomp, keepX = select_keepX)
  plotIndiv(final_splsda, comp = c(1, 2), ind.names = FALSE, legend = TRUE, ellipse = TRUE, title = 'SPLS-DA, final result, components 1 and 2')
  dev.new()
  if (final_ncomp > 2){
    plotIndiv(final_splsda, comp = c(1, 3), ind.names = FALSE, legend = TRUE, ellipse = TRUE, title = 'SPLS-DA, final result, components 1 and 3')
    dev.new()
    
    plotIndiv(final_splsda, comp = c(2, 3), ind.names = FALSE, legend = TRUE, ellipse = TRUE, title = 'SPLS-DA, final result, components 2 and 3')
    dev.new()
  }
  auroc(final_splsda, roc.comp = 1)
  dev.new()
  
  auroc(final_splsda, roc.comp = 2)
  dev.new()
  
  if (final_ncomp > 2){
    auroc(final_splsda, roc.comp = 3)
    dev.new()
  }
  
  final_perf = perf(final_splsda, validation = 'Mfold', folds = cv_folds, dist = 'max.dist', nrepeat = cv_repeats)
  matplot(final_perf$error.rate$BER, type = 'l', lty = 1, col = color.mixo(1:3), main = 'Balanced Error Rate of the final model')
  legend('topright', c('max.dist', 'centroids.dist', 'mahalanobis.dist'), lty = 1, col = color.mixo(1:3))
  dev.new()
  
  results = NULL
  results$splsda = final_splsda
  results$perf = final_perf
  return(results)
}

####### MAIN FUNCTION TO EXECUTE #######
apply_workflow = function(dataset, classes, fisher_variables = 5000, relieff_variables = 500, pca_components = 10, cv_folds = 5, cv_repeats = 10, debug=TRUE){
	# Process initial data
	data_matrix = as.matrix.data.frame(dataset)
	positives = which(classes == classes[1])
	negatives = which(classes != classes[1])

	after_fisher = apply_fisher(data_matrix, positives, negatives, debug = debug)

	after_relieff = apply_relieff(after_fisher, classes, debug = debug)

	apply_pca(after_relieff, classes, debug = debug)
	
	final_results = apply_plsda(after_relieff, classes, debug = debug)
	
	return(final_results)
}

