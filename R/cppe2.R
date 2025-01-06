library( ape ) 
library( Matrix )
library( ggtree ) 
library( ggplot2 )
library( glue )

.maketreedata <- function(tre, brquantile=0.01)
{
	n <- length( tre$tip.label)
	nnode = ape::Nnode(tre)

	# <based on treedater>
	tipEdges <- which( tre$edge[,2] <= n)
	i_tip_edge2label <- tre$tip.label[ tre$edge[tipEdges,2] ]
	daughters <- lapply( 1:(n+nnode), function(i) c() )
	parent <- rep(NA, n + nnode)
	for (k in 1:(n+nnode)){
		x <- tre$edge[which(tre$edge[,1] == k),2]
		if (length(x) > 0){
			daughters[[k]] <- c(daughters[[k]], x )
			for (u in x){
				if (!is.na(u)) parent[u] <- k
			}
		}
	}

	ndesc = sapply( daughters, length )
	nodetimes  <- ape::node.depth.edgelength( tre )
	maxHeight= rh <- max( nodetimes[1:n] )
	rhs = rh # for compatibility with multitree version 
	sts <- nodetimes[1:n]
	shs <- rh - sts 
	nhs = nodeheights = rh - nodetimes
	inhs <- sort( rh - nodetimes[ (n+1):(n + nnode) ] )

	# parent node heights 
	pnhs = nhs[ parent ] 

	coalescentcohorts = lapply( 1:nnode, function(j) {
		h = nhs[ n + j ]
		which( (nhs<h) & (pnhs>=h) ) 
	})

	whno = which( !is.na( parent ) )
	nr = length(whno)

	# precomp dep variable 
	nodeys = lapply( 1:ape::Nnode(tre), function(j){
		cohort= coalescentcohorts[[j]]
		dgtrs = daughters[[ j+n  ]]
		A = length(cohort) 
		if ( A == 1 ) return(0)
		y = rep( -(A/(2*(A-1)))*log(A-1), A)
		y[cohort %in% dgtrs] <- (A-1)*(A/(2*(A-1)))*log(A-1)   / 2 # TODO check /2 ; TODO check extra factor (A-1)
		y
	})
	nodey <- do.call( c, nodeys )
	
	# adjusted breanch lengths 
	brlens <- (nodetimes[whno]-nodetimes[parent[whno]])
	lbbrlen <- quantile( brlens[ brlens > 0 ], brquantile )
	brlens[ brlens <= 0 ] <- lbbrlen 

	rv = modifyList( tre 
		, list( n = ape::Ntip(tre), nedge = ape::Nedge( tre ), nnode = ape::Nnode( tre )
	 	 , daughters = daughters
	 	 , parent = parent
	 	 , nodetimes = nodetimes 
	 	 , nodeheights = nodeheights
	 	 , internalnodeheights = inhs 
	 	 , internalnodetimes = nodetimes[ (n+1):(n + nnode) ]
	 	 , parentnodeheights = pnhs
	 	 , ndesc = ndesc 
	 	 , coalescentcohorts = coalescentcohorts
	 	 , nr = nr 
	 	 , whno = whno 
	 	 , nodey = nodey 
	 	 , brlens = brlens 
	 	)
	)
	class(rv) <- c('phylo', 'cppephylo' ) 
	
	rv
}

.computenodew <- function( tr1, ipw  )
{
	stopifnot( inherits( tr1, 'cppephylo' ))
	nodey <- tr1$nodey 
	nodew <- rep( 1, length( nodey ))
	if( !is.null( ipw )){
		stopifnot( length(ipw )==ape::Ntip(tr1))
		stopifnot( is.numeric(ipw ) )
		if (!is.null(names(ipw ))){
			stopifnot( all( tr1$tip.label %in% names(ipw )))
			ipw <- ipw[ tr1$tip.label ]
		}
		ancestralw = ape::ace( ipw , tr1, type = 'continuous', method = 'pic' )$ace # pic for fast LS 
		# ancestralw <- rep( 1, tr1$nnode )
		ancestralw <- c( ipw, ancestralw )
		nodews <- lapply( 1:tr1$nnode, function(j){
			cohort= tr1$coalescentcohorts[[j]]
			dgtrs = tr1$daughters[[ j+tr1$n  ]]
			wj <- mean( ancestralw[ dgtrs ] )
			rep( wj, length( cohort ))
		})
		nodew <- do.call( c, nodews )
		nodew <- nodew / mean( nodew )
	}
	nodew 
}


#' Default arguments for tauprofile 
TPARGS <- list( logtaulb = -4, logtauub = 37, res = 21, startpc = 50, endpc = 100, nobj = 100 )

#' Fit a COD GMRF model using weighted least squares 
#' 
#' @param tr1 Phylogenetic tree in ape::phylo format 
#' @param tau Precision parameter. If NULL, will invoke `tauprofile` to find best value. 
#' @param profcontrol Optional list of arguments passed to `tauprofile`
#' @param inverseprobabilityweights An optional vector (named or unnamed) of sample weights for each tip in the input tree
#' @value A GPGMRF model fit 
#' @export 
codls <- function(tr1, logtau = c(0, NULL ), profcontrol = list(), inverseprobabilityweights=NULL )
{
	if (!inherits( tr1, 'cggephylo' ) & inherits(tr1,'phylo'))
	{
		tr1 <- .maketreedata(tr1)
	}
	logtau = logtau[1] 
	tpdf <- NULL 
	if ( is.null( logtau )){
		tpargs <- modifyList( TPARGS, profcontrol )
		tpargs$tr = tr1 
		tpargs$ipw = inverseprobabilityweights
		tpdf <- do.call( tauprofile, tpargs ) 
		logtau <- tpdf$logtau[ which.min( tpdf$loss ) ] 
		print( tpdf )
	}
	st1 <- Sys.time() 

	whno = tr1$whno 
	nr = tr1$nr 

	i <- which( !is.na( tr1$parent ) )
	ary = rep(0, nr)

	y = c( ary, tr1$nodey )

	arw = sqrt( exp(logtau)/tr1$brlens ) # sqrt(1/.) so variance prop to brlen 
	nodew <- .computenodew( tr1, inverseprobabilityweights )
	w = c( arw, nodew )
	st2b <- Sys.time() 

	# indices for X
	## AR component 
	ai = 1:nr
	aj = whno 
	ax = rep(1, nr )
	ai = c( ai, 1:nr )
	aj = c( aj, tr1$parent[whno] )
	ax = c( ax, rep(-1, nr) )

	# coindices is used by later functions, but is slow, maybe move this computation...
	k <- nr + 1 
	coindices <- list() 
	coii <- 1 
	for (co in tr1$coalescentcohorts)
	{
		coindices[[ coii ]] <- k:(k+length(co)-1) 
		coii <- coii + 1 
		k <- k + length(co) 
	}

	ncoi <- sum( sapply( tr1$coalescentcohorts, length ))
	coi <- (nr+1):(nr+ncoi)
	coj <- do.call( c, tr1$coalescentcohorts )
	cox <- rep(1, ncoi )
	ai <- c( ai, coi )
	aj <- c( aj, coj )
	ax <- c( ax, cox )

	X <- Matrix::sparseMatrix( i=ai, j=aj, x= ax  )
	W <- Matrix::Diagonal( x = w )
	QQ <- t(X) %*% W %*% X 
	b  <- t(X) %*% W %*% y 
	f2beta =  solve( QQ, b ) |> as.vector()

	st3 <- Sys.time() 

	structure( list( coef = f2beta 
		, logtau = logtau 
		, data = tr1 
		, X = X 
		, W = W 
		, y = y 
		, nr = nr 
		, arindices = 1:nr
		, logoddsindices = (nr+1):nrow(X)
		, istartnodeterms  = nr+1
		, coindices = coindices # indices in X corresponding to coalescentcohorts 
		, tauprofile = tpdf
		, runtime = st3  - st1 
	   )
	   , class = 'gpgmrf' )
}

#' @export 
print.gpgmrf <- function(f)
{
	stopifnot( inherits( f, 'gpgmrf' ))
	cat(' Genealogical placement GMRF model fit \n')
	print( f$data ) 
	cat('Range of coefficients: \n')
	print(range( coef(f)))
	cat(glue::glue( 'Precision parameter (log tau): {f$logtau} \n') )
	cat('\n')
	invisible(f) 
}

#' @export 
coef.gpgmrf <- function(f) 
{
	stopifnot( inherits( f, 'gpgmrf' ))
	f$coef 
}

#' @export 
plot.gpgmrf <- function( f )
{
	stopifnot( inherits( f, 'gpgmrf' ))
	f2beta = f$coef 
	tr1 = f$data 
	class(tr1) <- 'phylo' 
	fdf <- data.frame( node = 1:length(tr1$nodetimes), theta = f2beta) #coef(f1)  )
	gtr1 = ggtree::ggtree( tr1 ) %<+% fdf 
	gtr1 + aes(color = theta) + 
		scale_color_gradient2( low='blue'
			, mid = 'lightblue'
			, high = 'red'
			, midpoint = 0 # median( fdf$theta )
			, limits = range(fdf$theta)
			, name = "ψ" ) + 
		ggtree::theme_tree2() +
		ggtree::geom_tiplab()
}





#' Evaluate the loss function in the GPEGMRF model across a range of tau (precision parameter) values
#' 
#' 
#' 
#' @param tr A phylogenetic tree in ape::phylo format 
#' @param logtaulb Lower bound of precision parameteters 
#' @param logtauub Upper bound of precision parameteters 
#' @param res Number of tau values to evaluate 
#' @param startpc The initial per cent of nodes in the tree counting from root to tips where the loss function will be evaluated 
#' @param endpc The final per cent of nodes in the tree counting from root to tips where the loss function will be evaluated 
#' @param nobj The integer number of points along the tree where the loss function will be evaluated. If Inf, will use all points between startpc and endpc, but may be slow. 
#' @param ipw Optional inverse probability weights for each sample 
#' @value A data frame containing the loss function evaluated over a range of tau values 
#' @export 
tauprofile <- function(tr , logtaulb = -4, logtauub = 35, res = 11, startpc = 75, endpc = 100, nobj = 100, ipw = NULL ) 
{

	if (!inherits(tr, 'cggephylo' ) & inherits(tr,'phylo'))
	{
		tr <- .maketreedata(tr)
	}

	stopifnot( logtaulb < logtauub ) 
	logtaus = seq( (logtaulb), (logtauub), length = res )  

	# cohorts in order 
	nodeorder <- order( tr$internalnodetimes, decreasing=FALSE) 

	# indices for forecasting 
	starticohorts = floor( startpc*.01*length(nodeorder) )
	endicohorts = floor( endpc*.01*length(nodeorder) )-1
	icohorts <- starticohorts:endicohorts
	if ( nobj < length( icohorts )){
		icohorts <- seq( starticohorts, endicohorts, by = floor((endicohorts-starticohorts+1)/nobj))
	}

	logtau0 <- (logtaulb + logtauub)/2
	f = codls(tr, logtau0 ) # 

	nodew <- .computenodew( tr, ipw )

	loss <- c() 
	for (logtau in logtaus )
	{
		losses <- sapply( icohorts, function(i){
			keepinds <- c( 1:tr$nr, do.call(c, f$coindices[nodeorder[1:i]] )  )
			X1 = f$X[ keepinds, ]

			arw = sqrt( exp(logtau)/tr$brlens ) # sqrt(1/.) so variance propto brlen
			# nodew <- rep( 1, length(f$logoddsindices ))  # f used here 
			w <- c( arw, nodew )
			W1 <- Matrix::Diagonal( x = w )[keepinds, keepinds ] 
			y1 = f$y[ keepinds ] # f used here 
			QQ1 <- t(X1) %*% W1 %*% X1 
			b1 <- t(X1) %*% W1 %*% y1 
			beta1 =  solve( QQ1, b1 ) |> as.vector()
			sum( (f$y[f$coindices[[nodeorder[i+1]]]] - f$X[f$coindices[[nodeorder[i+1]]],] %*% beta1)^2 )
		})
		loss <- c( loss, mean(losses))
	}
	
	odf <- data.frame( logtau = logtaus, loss = loss )
	# o= optimise( ofun, lower = logtaulb , upper = logtauub , maximum = FALSE ) # optimise?
	odf 
}

.optimcodgmrf <- function( f , ...)
{

	coy = ifelse( f$y[ f$logoddsindices] > 0 , 1 , 0 )
	ip = apply( f$X[ f$logoddsindices, ], MAR=1, FUN=function(x)which(x==1)[1] )

	A <- rep(1, length(coy))
	for (co in f$coindices){
		a <- length( co ) 
		A[ co-f$istart+1 ] <- a
	}
	psiintercept <- -log(A-1)
	lbbrlen <- quantile(f$data$brlens[f$data$brlens>0], .01)
	brlens <- pmax( f$data$brlens , lbbrlen )

	ofun <- function(psi, logtau = f$logtau)
	{
		arterms <- dnorm( as.vector( f$X[f$arindices,] %*% psi ), 0
			, sd =  sqrt( brlens/exp(logtau)) 
			, log=TRUE )

		psi1 = psi[ ip ]
		pp <- 1 / (1 + exp(-(psi1 + psiintercept)))
		coodterms <-  coy*log(pp) + (1-coy)*log(1-pp)  

		ll1 = sum( arterms ) 
		ll2 = sum( coodterms )
		# print(c( ll1, ll2) )

		-(ll1 + ll2)
	}
	# ofun( f$coef )
	# , control = list(fnscale=-1)
	o = optim( par = f$coef, fn = ofun, method = 'BFGS', ...)
	# print( summary( o$par ))
	# print( summary( f$coef ))
	# plot( f$coef, o$par ); abline( a = 0, b = 1 )

	f$coef <- o$par 
	f$optim <- o 
	f
}

#' Fit a genealogical placement GMRF model using maximum likelihood 
#' 
#' This method optimises the sequential-bernoulli likelihood of a COD GMRF model using gradient descent and using `codls` to find initial conditions. 
#' If tau is not provided, `codls` is also used to optimise this parameter. 
#' This method is slower than `codls` and is not recommended for trees with more than several hundred samples. 
#' 
#' The ML COD GMRF model does not support inverse probability weighting of samples. Use `codls` if sample weighting is needed.
#' 
#' @param tr1 Phylogenetic tree in ape::phylo format 
#' @param tau Precision parameter. If NULL, will invoke `tauprofile` to find best value. 
#' @param profcontrol Optional list of arguments passed to `tauprofile`
#' @param ... Additional arguments are passed to `optim`
#' @value A GPGMRF model fit 
#' @export 
codml <- function(tr1, logtau = c(0, NULL ), profcontrol = list(), ... )
{
	f = codls( tr1, logtau, profcontrol )
	.optimcodgmrf(f, ... )
}


#' Compute phylogenetic clusters by cutting tree at branches with large changes in coalescent odds 
#' 
#' @param f A model fit from `codls`
#' @param clth Numeric threshold change in coalescent log odds 
#' @param rescale if TRUE (default), coalescent log odds are rescaled (mean zero, unit variance) prior to applying thresholds
#' @value A data frame with cluster asignment for each tip 
#' @export 
computeclusters <- function(f, clth, rescale=TRUE)
{
	clusters <- list() 
	tr <- f$data 
	if (rescale)
		scpsi <- scale( f$coef ) 
	tocut <- apply( tr$edge, MAR=1, FUN=function(uv){
		u = uv[1] 
		v = uv[2] 
		(scpsi[v] - scpsi[u]) > clth 
	}) |> which() 
	table( tocut )
	ntd <- lapply(  1:(Ntip(tr)+Nnode(tr)), function(x) c() )
	clusters <- list() 
	nadded <- rep(0, Ntip(tr) + Nnode(tr))
	for (i in 1:Ntip(tr)) ntd[[i]] <- i 
	for (i in (1+Ntip(tr)):(Ntip(tr)+Nnode(tr)) ) ntd[[i]] <- NA
	for ( ie in ape::postorder(tr))
	{
		u <- tr$edge[ie,1]
		v <- tr$edge[ie,2]
		tryCatch( ntd[[v]], error = function(e) browser() )
		ntd[[u]] <- c( ntd[[u]], ntd[[v]] ) 
		nadded[u] <- nadded[u]+1
			# if ( ie %in% tocut ) browser() 
		if ( (ie %in% tocut) & (nadded[u]==tr$ndesc[u]))  
		{
			clusters[[length(clusters)+1]] <-  na.omit( c( ntd[[u]], ntd[[v]] ) )
			ntd[[u]] <- NA 
		} 
	}
	cl1 <- setdiff( 1:Ntip(tr), na.omit( do.call(c,clusters)) )
	clusters[[length(clusters)+1]] <-  cl1
	clusterdf <- data.frame() 
	for (ic in seq_along( clusters )){
		cl <- clusters[[ic]]
		if( length( cl ) > 0 )
		{
			clusterdf <- rbind( clusterdf 
			, data.frame( tip = cl, tip.label = tr$tip.label[ cl ], clusterid = ic, psi = f$coef[cl] )
			)
		}
	}
	clusterdf$clusterid <- as.factor( clusterdf$clusterid )
	clusterdf
}

#' Plots a fit from `codls` and cluster assigment from `computeclusters`
#'
#' @export 
plotclusters <- function( f, clusterdf, ... )
{
	library( ggtree ) 
	library( ggplot2 )
	library( dplyr )
	cmat <- as.matrix( clusterdf[, c('clusterid') ]); colnames(cmat) <- 'cluster'
	rownames(cmat) <- clusterdf$tip.label 
	gheatmap( plot(f) , cmat, ... )
}


#' Computes the CH index for selecting clustering thresholds
#' 
#' @param clusterdf Output from `computeclusters`
#' @export 
chindex <- function( clusterdf )
{
	cldfs <- split( clusterdf, clusterdf$clusterid )
	if ( length( cldfs ) == 1 ) return(0) 
	sapply( cldfs, function(x) mean(x$psi) ) -> clpsi 
	sapply( cldfs, nrow ) -> cln
	opsi <- mean( clusterdf$psi ) 
	bcss <- sum( cln * (clpsi-opsi)^2 ) 
	wcss <- sapply( cldfs, function(x) mean((x$psi - mean(x$psi))^2)) |> sum() 
	n <- nrow(clusterdf ) 
	k <- length( cldfs ) 
	(bcss/(k-1)) / (wcss/(n-k)) 
}

#' Compute CH index over a range of clustering thresholds
#' 
#' @param f A fit from `codls`
#' @param clths Numeric vector of clustering thresholds
#' @param rescale Passed to `computeclusters`
#' @export 
chindices <- function(f, clths = seq( .1, 1.5, length = 20 ), rescale=TRUE)
{
	chs <- sapply( clths, function(clth){
		cldf <- computeclusters( f, clth , rescale=rescale)
		chindex(  cldf )
	})
	data.frame( threshold = clths, CH = chs ) 
}
