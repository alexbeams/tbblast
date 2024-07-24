# Taken directly from Rasigade's github:
#	Use the functions defined herein to calcualte
#	the Timed Haplotype Density



#
# Timescaled haplotypic density - estimates epidemic
# success from genotypic data
# Copyright (C) 2015-2016 Jean-Philippe Rasigade
# Contact: jean-philippe.rasigade@univ-lyon1.fr
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/agpl.html>
#

# NOTE
# Code documentation adheres to roxygen2 format
# Example script in Section 4

#
# SECTION 1 SUPPORTING FUNCTIONS
#
# 

#' @title Hamming distance
#' @param x matrix of microsatellites or other categorical data
#' @export
hamming <- function(x) {
	n  <- nrow(x)
	gd <- matrix(0, n, n)
	for(i in 1:(n-1)) for(j in (i+1):n)
		gd[i,j] <- gd[j,i] <- sum(x[i,] != x[j,])
	gd
}

#' @title Truncated geometric distribution
#' @param x vector of quantiles
#' @param alpha geometric parameter
#' @param p truncation limit (maximum value of x)
#' @export
dtgeom <- function(x, alpha, p) {
	k <- (1 - alpha) / (1 - alpha^(p+1)) * alpha^x
	k[!is.finite(k)] <- 1/p
	return(k)
}

#' @title Cumulative truncated geometric function
#' @param x value (quantile)
#' @param a geomatric parameter
#' @param m number of markers
#' @return cumulative probability P(X < x)
#' @export
ptgeom <- function(x, a, m) (1 - a^x) / (1 - a^(m))

#
# SECTION 2 TIMESCALE / BANDWIDTH
#
# 

#' @title TMRCA maximum-likelihood estimate under Infinite Alleles Model
#' @param h Hamming distance (number of different sites)
#' @param m total number of sites
#' @param mu per-site mutation rate
#' @details After Walsh B. Genetics 2001, Eq. (5) p. 898.
#' @return vector of tMRCAs; the scale (years, generations, etc) depends on the scale of mu
#' @export
tmrca_iam <- function(h, m, mu) 1/(2*mu) * log(m/(m-h))

#' @title Inverse TMRCA relation
#' @param t TMRCA
#' @param m number of markers
#' @param mu mutation rate per marker
#' @return expected number of differences (rational)
#' @export
tmrca_iam_inv <- function(t, m, mu) m * exp(-2*t*mu) * (exp(2*t*mu) - 1)

#' @title Adjust THD bandwidth to a specific TMRCA quantile
#' @param t tMRCA
#' @param m number of markers
#' @param mu mutation rate
#' @param q target cumulative value
#' @return geometric bandwidth value such that the given TMRCA is the q'th quantile
#' of the geometric kernel function
#' @export
tmrca2bandwidth <- function(t, m, mu, q = .5) {
	# Expected number of differences
	h <- tmrca_iam_inv(t, m, mu)
	# Target function for quantile q
	f <- function(b) (ptgeom(h, b, m) - q)^2
	# Optimize
	b <- optimize(f, c(1e-6,1-1e-6))$minimum
	b
}

#
# SECTION 3 THD COMPUTATION
#
# 

#' @title Timescaled haplotypic density
#' @param h Hamming distance matrix
#' @param alpha geometric parameter
#' @param p truncation limit (maximum value of x)
#' @param from subset vector (boolean)
#' @param to subset vector (boolean)
#' @param skipself if TRUE (default), individuals do not contribute to their own density 
#' @return Vector of densities
#' @export
thd <- function(h, alpha = 0.5, p = max(h), from, to, skipself = T) {
	if(missing(from) & missing(to)) {
		# Pairwise densities
		n <- nrow(h)
		k <- NULL 
		if(skipself) 
			k <- sapply(1:n, function(i) {
					mean(dtgeom(h[,i][-i], alpha, p))
				 })
		else
			k <- sapply(1:n, function(i) {
					mean(dtgeom(h[,i], alpha, p))
				 })
		return(k)		
	} else {
		# Densities from one group to another
		# Check whether from = to and skip diagonal
		if(missing(from)) stop("Argument 'to' cannot be provided alone")
		if(missing(to)) to <- from
		from[is.na(from)] <- FALSE
		to[is.na(to)]     <- FALSE
		skipdiag <- all(from == to)
		hh <- h[to, from] # Scan columns rather than rows for efficiency
		if(sum(to) == 1)   hh <- matrix(hh, 1, sum(from))
		if(sum(from) == 1) hh <- matrix(hh, 1, sum(to))
		n <- ncol(hh)
		kdefunc <- if(skipdiag) {
			function(i, hh, alpha, p) {mean(dtgeom(hh[,i][-i], alpha, p))}
		} else {
			function(i, hh, alpha, p) {mean(dtgeom(hh[,i], alpha, p))}
		}
		k <- sapply(1:n, kdefunc, hh = hh, alpha = alpha, p = p)
		return(k)
	}
}

#
# SECTION 4 EXAMPLE
#
# 

## Simulate matrix of 15-loci minisatelite data for 200 individuals
#n <- 200
#m <- 15
#divergence <- 0.1 # Poisson parameter
#d <- matrix(rpois(n*m, divergence), n, m)
#
## plot(hclust(dist(d)))
#
## Specify evolutionary rate (e.g. per locus per year) and timescale (years)
#mu <- 5e-4
#timescale <- 50
#
## Compute bandwidth
#bandwidth <- tmrca2bandwidth(timescale, m, mu)
#
## Extract Hamming distance matrix
#H <- hamming(d)
#
## Compute THD
#THD <- thd(H, bandwidth, m)
#
## Display tree and THDs
#hc <- hclust(as.dist(H))
#
#par(mfrow = c(2,1))
#par(mar = c(1, 5, 1, 5))
#plot(hc, labels = FALSE, xlab = "", main = "", sub = "")
#barplot(THD[hc$order], ylim = rev(range(THD)), ylab = "THD")

