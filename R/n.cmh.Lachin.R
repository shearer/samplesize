n.cmh.Lachin<-function(power=0.8,alpha=0.05,eta1=0.5,pi1,pi2) {
	# Calculates sample size for the Cochran-Mantel-Haenszel test of a 2xJ dataset (2 groups, with J ordered categories) using the method of Lachin, Statistics in Medicine 30(2011):3057-66.
	# Sample sizes are calculated for two cases: (1) Assuming the scores are 1:n_categories and comparing their medians using the Cochran-Mantel_Haenszel test, and (2) Assuming the score ranks will be compared using the Wilcoxon (Mann-Whitney) test.
	# First a few checks
	if (power <= 0 | power >= 1) {
		stop("Power must be a numeric value between 0 and 1")
	}
	if (alpha <= 0 | alpha >= 1) {
		stop("alpha must be anumeric value between 0 and 1")
	}
	if (eta1 <= 0 | eta1 >= 1) {
		stop("t must be a numeric value between 0 and 1")
	}
	npi1 <- length(pi1)
	npi2 <- length(pi2)
	if (npi1 != npi2) {
		stop("pi1 and pi2 must be vectors of equal length")
	}
	if (npi1 < 2) {
		stop("pi1 and pi2 must be vectors of with at least two elements")
	}
	
	# Set the parameters
	eta2<-1-eta1
	pi<-eta1*pi1+eta2*pi2
	sc<-c(1:npi1)  #Scoring version assumes simple numerical score per category
	rnk<-c(1:npi1) #Initialise rnk
	sumPi=0
	for (i in rnk) {
		rnk[i]<-sumPi+pi[i]/2         #Eq 8, assumes scoring by ranks
		sumPi<-sumPi+pi[i]
	}
	
	sig_pi<-diag(pi)-pi%o%pi         #Covariance of pi based on Eq 2
	tau2_sc<-(sc%*%(pi1-pi2))^2/(sc%*%sig_pi%*%sc*(1/eta1+1/eta2))  #tau2 based on Eq 7
	tau2_rnk<-(rnk%*%(pi1-pi2))^2/(rnk%*%sig_pi%*%rnk*(1/eta1+1/eta2))  #tau2 based on Eq 9
	
	nSample_sc<-ceiling(gtx::chi2ncp(alpha,power)/tau2_sc)        #Package gtx has a function to determine noncentrality parameter for Chi2
	n1_sc <- round(nSample_sc * eta1, 0); n2_sc <- round(nSample_sc * eta2, 0)
	
	nSample_rnk<-ceiling(gtx::chi2ncp(alpha,power)/tau2_rnk)
	n1_rnk <- round(nSample_rnk * eta1, 0); n2_rnk <- round(nSample_rnk * eta2, 0)
	return(data.frame('Method'=c('Score','Rank'),'Total sample size'=c(nSample_sc,nSample_rnk),'Group 1'=c(n1_sc,n1_rnk),'Group 2'=c(n2_sc,n2_rnk)))
}
