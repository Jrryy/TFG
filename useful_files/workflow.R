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

load("dataframes.RData")

####### PREPARATION #######

mas5_matrix = as.matrix.data.frame(mas5)

ALS = which(classes == "ALS")
noALS = which(classes != "ALS")

####### FISHER #######

# to calculate Fisher score of each feature: 
fisher = function(gene, pos, neg) {
	score = (mean(gene[pos]) - mean(gene[neg]))**2 / (length(pos)*var(gene[pos]) + length(neg)*var(gene[neg]))
    	return(score)
}

apply_fisher = function(dataset, positive, negative, features_to_keep = 5000){

	scores = apply(dataset, 2, fisher, positive, negative)

    	print(length(scores))

    	stopifnot(length(scores) == dim(dataset)[2])

    	sorted_scores = sort(scores, decreasing=TRUE, index.return=TRUE)
    	plot(sorted_scores$x)

	dev.new()

    	# We keep the best x features according to fda, 5000 by default but the user should be able to input it.
    	after_fisher = dataset[, sort(sorted_scores$ix[0:features_to_keep])]

	return(after_fisher)
}


####### RELIEFF #######

set.seed(1337)

# Applies ReliefF. For now, only ReliefFequalK, but I plan on using exp too and adding a parameter to specify it
apply_relieff = function(dataset, classes, features_to_keep = 500, iterations = 0){
	
	stopifnot(length(classes) == dim(dataset)[1])

	df_ = as.data.frame(dataset)
	df = cbind(df_, class=classes)

	reliefF_attrs = attrEval('class', df, estimator='ReliefFexpRank', ReliefIterations=iterations)

	sorted_attrs = sort(reliefF_attrs, decreasing=TRUE, index.return=TRUE)

	plot(sorted_attrs$x)
	dev.new()
	
	after_relieff = dataset[, sort(sorted_attrs$ix[0:features_to_keep])]

	heatmap(after_relieff)
	dev.new()

	return(after_relieff)
}

####### PCA #######

apply_pca = function(dataset, classes, components = 10){
	pca_plot = pca(dataset, ncomp=components, center=TRUE, scale=TRUE)
	plot(pca_plot)
	dev.new()

	plotIndiv(pca_plot, comp=c(1, 2), ind.names=FALSE, group=classes, legend=TRUE, ellipse=TRUE)
}

####### MAIN FUNCTION TO EXECUTE #######
apply_workflow = function(dataset, classes, fisher_variables = 5000, relieff_variables = 500, pca_components = 10){
	# Process initial data
	data_matrix = as.matrix.data.frame(dataset)
	positives = which(classes == classes[1])
	negatives = which(classes != classes[1])

	after_fisher = apply_fisher(data_matrix, positives, negatives)

	after_relieff = apply_relieff(after_fisher, classes)

	apply_pca(after_relieff, classes)
}

