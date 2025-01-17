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
		y[cohort %in% dgtrs] <- (A-1)*(A/(2*(A-1)))*log(A-1)   / 2 # check /2 ; check extra factor (A-1)
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


# Default arguments for tauprofile 
TPARGS <- list( logtaulb = -4, logtauub = 37, res = 11, startpc = 50, endpc = 100, nobj = 100 )

#' Fit a COD GMRF model using weighted least squares 
#' 
#' @param tr1 Phylogenetic tree in ape::phylo format 
#' @param logtau Precision parameter. If NULL, will invoke `tauprofile` to find best value. 
#' @param profcontrol Optional list of arguments passed to `tauprofile`
#' @param weights An optional vector (named or unnamed) of sample weights for each tip in the input tree
#' @return A GPGMRF model fit 
#' 
#' @examples 
#' # A simple example that does not have population structure 
#' set.seed( 1111 )
#' tr <- ape::rtree( 100 )
#' f <- codls(tr)
#' coef(f) |> head() 
#' summary(f)
#' \dontrun{
#' plot(f)
#' }
#' 
#' # This example has population structure 
#' tr0 = rcoal(20); tr0$edge.length <- .01*tr0$edge.length 
#' tr1 = rcoal(80); 
#' dx <- (max(node.depth.edgelength( tr1 ))-max(node.depth.edgelength( tr0 )))
#' tr0$root.edge <- dx
#' tr <- bind.tree(tr0,tr1, position = dx)
#' f <- codls(tr)
#' summary(f) 
#' \dontrun{
#' plot(f)
#' }
#' 
#' @export 
codls <- function(tr1, logtau = NULL , profcontrol = list(), weights=NULL )
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
		tpargs$ipw = weights
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
	nodew <- .computenodew( tr1, weights )
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
	QQ <- Matrix::t(X) %*% W %*% X 
	b  <- Matrix::t(X) %*% W %*% y 
	f2beta =  Matrix::solve( QQ, b ) |> as.vector()

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

#' Root mean square coalescent log odds 
#'
#' This is a summary statistic that describes the amount of variation in coalescent rates across lineages in a phylogenetic tree  
#' It is defined as \deqn{ \sqrt{ \sum_i l_i \psi_i^2 / L }  }
#' where the sum is over all branches i in the tree and weighted by branch length \eqn{l_i} and where  \eqn{ L = \sum_i l_i  }
#'
#' @param f Fit from `codls`
#' @return Numeric RMSCLO 
#' @export 
rmsclo = phylopredictance <- function(f) 
{
	L = sum( f$data$brlens ) 
	l = f$data$brlens[match(1:length(coef(f)), f$data$edge[,2])]
	l[ is.na(l)] <- 0
	L = sum(l)
	sqrt(  sum( l*(coef(f)-mean(coef(f)))^2 ) / L  )
}

#' Effective number of extant lineages 
#' 
#' @param f Fit from `codls` 
#' @return Numeric effective number extant lineages 
#' @export 
neffextant <- function(f) 
{
	p <- 1 / (1 + exp(-coef(f)[1:ape::Ntip(f$data)]))
	p <- p / sum(p) 
	1 / sum(p^2)
}

#' @export 
summary.gpgmrf <- function(f)
{
	print(f)
	cat('\n' )
	odf <- data.frame( logprecision = f$logtau
		 , RMSCLO = rmsclo(f) 
		 , Neff = neffextant(f) 
		)
	print( odf )
	cat('\n' )
	invisible(odf)
}

#' @export 
print.gpgmrf <- function(f)
{
	stopifnot( inherits( f, 'gpgmrf' ))
	cat(' Genealogical placement GMRF model fit \n')
	print( f$data ) 
	cat('Range of coefficients: \n')
	cat( glue::glue('{range(coef(f))}') )
	cat( '\n' )
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
	fdf <- data.frame( node = 1:length(tr1$nodetimes), theta = f2beta)
	gtr1 = ggtree::`%<+%`( ggtree::ggtree( tr1 ), fdf )
	gtr1 + ggtree::aes(color = theta) + 
		ggplot2::scale_color_gradient2( low='blue'
			, mid = 'lightblue'
			, high = 'red'
			, midpoint = 0 # median( fdf$theta )
			, limits = range(fdf$theta)
			, name = "psi" ) + 
		ggtree::theme_tree2() +
		ggtree::geom_tiplab()
}





#' Evaluate the loss function of the cod model across a range of tau (precision parameter) values
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
#' @param weights Optional inverse probability weights for each sample 
#' @return A data frame containing the loss function evaluated over a range of tau values 
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
			QQ1 <- Matrix::t(X1) %*% W1 %*% X1 
			b1 <- Matrix::t(X1) %*% W1 %*% y1 
			beta1 =  Matrix::solve( QQ1, b1 ) |> as.vector()
			sum( (f$y[f$coindices[[nodeorder[i+1]]]] - f$X[f$coindices[[nodeorder[i+1]]],] %*% beta1)^2 )
		})
		loss <- c( loss, mean(losses))
	}
	
	odf <- data.frame( logtau = logtaus, loss = loss )
	odf$optimal <- ''
	odf$optimal[ which.min( odf$loss ) ] <- '***' 
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
	o = optim( par = f$coef, fn = ofun, method = 'BFGS', ...)

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
#' @return A GPGMRF model fit 
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
#' @return A data frame with cluster asignment for each tip 
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
	ntd <- lapply(  1:(ape::Ntip(tr)+ape::Nnode(tr)), function(x) c() )
	clusters <- list() 
	nadded <- rep(0, ape::Ntip(tr) + ape::Nnode(tr))
	for (i in 1:ape::Ntip(tr)) ntd[[i]] <- i 
	for (i in (1+ape::Ntip(tr)):(ape::Ntip(tr)+ape::Nnode(tr)) ) ntd[[i]] <- NA
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
	cl1 <- setdiff( 1:ape::Ntip(tr), na.omit( do.call(c,clusters)) )
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
#' @param A `codls` fit 
#' @param Output data frame from `computeclusters`
#' @export 
plotclusters <- function( f, clusterdf, ... )
{
	stopifnot( 'ggtree' %in% installed.packages()[,1] )
	cl <- sort( unique( clusterdf$clusterid ))
	cmat <- sapply( cl, function(x) (clusterdf$clusterid == x ) )
	colnames(cmat ) <- cl 
	rownames(cmat) <- clusterdf$tip.label 
	cmat <- cmat[ f$data$tip.label, ]
	ggtree::gheatmap( plot(f), cmat, ... )
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
	chdf <- data.frame( threshold = clths, CH = chs, optimal='' ) 
	chdf$optimal[ which.max(chdf$CH) ] <- '***'
	chdf
}


#' Automatically reweight sample units that may be oversampled 
#' 
#' Given a fit from `codls` and a set of samples which are suspected of over-sampling, this function will 
#' re-compute `codls` over a range of reweighted samples. This will identify the sample weight at which 
#' an association is lost between coalescent odds (psi) and the given set of samples. This is an appropriate
#' weight to use if there is an association between coalescent odds and psi that is due to sampling effects and 
#' not due to evolutionary effects, but note that this method may mask evolutionary effects if any are present.
#' 
#' @param f A `codls` fit 
#' @param rwtips Vector of samples (type character) which are suspected of over-sampling
#' @param wlb Numeric lower bound of sample weights to examine 
#' @param wub Numeric upper bound of sample weights to examine 
#' @param res Integer number of weights to examine 
#' @param alpha The p value threshold used for selecting the optimal weight 
#' @return A list with components `fit`: the reweigthed `codls` fit; `weights`: a new vector of sample weights; `optimalweight` the scalar weight applied to oversampled units; and `summary`: a data frame showing regression p values over a range of weights 
#' @export 
autoreweight <- function(f, rwtips, wlb = .01, wub = .5, res = 10, alpha = .05 ) {
	tr <- f$data 
	V <- rep(FALSE, ape::Ntip(tr)) |> setNames(tr$tip.label )
	V[ reweighttips ] <- TRUE 
	w <- rep(1, ape::Ntip(tr)) |> setNames( tr$tip.label)
	ws <- seq( wlb, wub, length=res )
	slms <- lapply( ws, function(ww) {
		w[ reweighttips ] <- ww 
		ff = codls( tr, logtau = f$logtau , weights = w)
		slm <- lm( coef(ff)[1:ape::Ntip(tr)] ~ V) |> summary() 
		slm
	} )
	ps <- sapply( slms, function(x) 
		x$coefficients[2,4]
	)
	sdf <- data.frame( sampleweight=ws, p = ps )
	if( max(sdf$p)>alpha ){
		optweight <- max( sdf$sampleweight[ sdf$p > alpha ] )
		w[ reweighttips ] <- optweight
		ff <- codls( tr, logtau = f$logtau , weights = w)
	} else{
		optweight <- NA 
		ff <- f 
		w<- rep(1, ape::Ntip(tr)) |> setNames( tr$tip.label)
	}
	
	list(
		fit = ff 
		, weights = w 
		, optimalweight = optweight 
		, summary = sdf 
	)
}
