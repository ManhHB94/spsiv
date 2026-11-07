********************************************************************
*	PROGRAM GENERATES SYNTHETIC INSTRUMENTAL VARIBALES
*				IN SPATIAL REGRESSION
*					-----------
*	According to Le Gallo & Paez (2013) and Fingleton (2022)
*					-----------
*	Version 1.3		Manh Hoang Ba (hbmanh9492@gmail.com)
*					
********************************************************************
*	Version history
*		02 Sep 25: 1.0 - first time contribution
*		22 Oct 25: 1.1 - upspeed and modify results
*		06 Oct 26: 1.2 - disable "mata clear"
*					   - correct tempname `p_r' in _eigvecvars
*					   - allow spmatrix object
*		07 Nov 26: 1.3 - upspeed:
*					     . no reshape/save/merge: panel is arranged
*					       as N x T matrix directly in Mata
*					     . eigenvectors kept in Mata (no svmat, no
*					       matsize limit, no E* variables)
*					     . correlations with all E_i in one matrix
*					       product; critical |r| replaces p-values
*					     . double-centering instead of QWQ products
*					   - if/in are marked together with varlist (marksample)
*					     and with the id/time variables
*					   - free Mata objects after use
********************************************************************

cap program drop spsiv
program define spsiv
	version 11.1	// According to spmat.ado
	
	syntax varlist(min=1 numeric) [if] [in] , Mmat(name) [Alpha(real 0.05)]

	if `alpha'<=0 | `alpha'>=1 {
		di in red "{bf:alpha()} must be between 0 and 1"
		exit 198
	}

*	Remove duplicate variables from varlist
	local varlist : list uniq varlist

*	Check spset
	qui cap confirm v _ID
	if _rc>0 {
		di in red "data must be {bf:spset}"
		exit 198
	}

	else {
*	Check xtset/tsset	
		qui cap xtset
		if _rc>0 {			// Cross-sectional data
			qui cap spset
			if _rc>0 {
				di in red "data must be {bf:spset}"
				exit 198
			}
		}
		else {				// Panel data
			if "`r(timevar)'"=="" {
				di in red "time variable not set; use {bf:xtset} {it:panelvar} {it:timevar}"
				exit 198
			}
			local ivar `r(panelvar)'
			local tvar `r(timevar)'
			qui xtset `ivar' `tvar'
			qui spset, modify nocoord
		}
	}
*
	local sp_id_var `r(sp_id_var)'

*	Mark sample: if/in, non-missing varlist, non-missing id and time
	marksample touse
	markout `touse' `sp_id_var' `tvar'
	qui count if `touse'
	if r(N)==0 error 2000

*	Temporary variables to hold the synthetic instruments
	local k : word count `varlist'
	local siv_tmp
	forvalues j = 1/`k' {
		tempvar s`j'
		local siv_tmp `siv_tmp' `s`j''
	}

*	Calculate eigenvectors E (kept in MATA)
	tempname Emat
	cap noi _eigvecvars, w(`mmat') emat(`Emat')
	if _rc {
		local rc = _rc
		cap mata: mata drop `Emat'
		exit `rc'
	}

*	Calculate Synthetic instrument for varlist (MATA)
	cap noi mata: _spsiv_core(`Emat', "`varlist'", "`siv_tmp'", ///
		"`touse'", "`sp_id_var'", "`tvar'", `alpha')
	local rc = _rc
	cap mata: mata drop `Emat'
	if `rc' exit `rc'

*	Name the new variables and calculate correlations
	local j = 0
	foreach name of local varlist {
		local ++j
		local s : word `j' of `siv_tmp'
		qui cap confirm v siv_`name'
		if _rc==0 qui drop siv_`name'
		rename `s' siv_`name'
		qui corr `name' siv_`name' if `touse'
		scalar corr_`name' = r(rho)
	}

	
*	Display results
	di _n in ye "Correlation between X and synthetic intrumental variable"

	local k : word count `varlist'
	local hline = 12*(`k'+1)

	// Loop through first 6 variables
	local groups = ceil(`k'/6)
	forvalues g = 1/`groups' {
		local start = (`g'-1)*6 + 1
		local end = min(`start'+5, `k')
		
		// First line: title "Variables" and varnames
		di as text "{hline `hline'}"
		di _column(1) "Variable" _continue
		
		local col = 15		
		forvalues j = `start'/`end' {
			local var : word `j' of `varlist'
			di _column(`col') %9s "`= abbrev("`var'", 9)'" _continue
			local col = `col' + 12
		}		
		
		di ""
		di as text "{hline `hline'}"
		
		// Second line: "Corr(X,SIV)" and correlation coef.
		di _column(1) "Correlation" _continue
		local col = 15
		forvalues j = `start'/`end' {
			local var : word `j' of `varlist'
			di _column(`col') %9.4f in ye corr_`var' _continue
			local col = `col' + 12
		}
		di ""
		di as text "{hline `hline'}"
		
		// If not the last group, print a blank line to space the next table.
		if `g' < `groups' {
			di ""
		}
	}
end


*********************************************************************
*	Program calculate Eigenvectors variables
*		from a Symmetric Adjacency matrix
*		(result is a MATA matrix named in emat())
*********************************************************************
cap program drop _eigvecvars	
program define _eigvecvars
	version 11.1
	syntax , Wmatrix(name) Emat(name)

	tempname w_m idv
	
*	Check spmat object
	qui cap spmat summarize `wmatrix'
	if _rc > 0 {
		qui cap spmatrix summarize `wmatrix'
		if _rc > 0 {
			di as err "{it:`wmatrix'} not found. It must be a {bf:spmat} or {bf:spmatrix} object"
			exit 498
		}
	}
	
*	Copy weighting matrix to MATA	
	cap spmat getmatrix `wmatrix' `w_m'
	if _rc > 0 {
		cap spmatrix matafromsp `w_m' `idv' = `wmatrix' 
		if _rc > 0 {
			di as err "cannot copy {it:`wmatrix'} to Mata"
			exit 498
		}
	}
	
*	Eigenvectors of P = (I-11'/N)W(I-11'/N)
	cap noi mata: `emat' = _spsiv_E(`w_m')
	local rc = _rc
	cap mata: mata drop `w_m'
	cap mata: mata drop `idv'
	if `rc' exit `rc'
end


*********************************************************************
*	MATA functions
*********************************************************************

mata:

//	Eigenvectors of the double-centered symmetric matrix W
	real matrix _spsiv_E(real matrix W)
	{
		cm = mean(W)				// 1 x r column means
		if (rows(W) != cols(W)) {
			errprintf("Weighting matrix must be square\n")
			exit(498)
		}
		if (!issymmetric(W)) {
			errprintf("Weighting matrix must be symmetrical\n")
			exit(498)
		}
		gm = mean(cm')				// grand mean

		// (I-11'/r) W (I-11'/r) = W - cm - cm' + gm  when W is symmetric
		N = rows(W)
		P = W :- J(N, 1, 1)*cm :- cm'*J(1, N, 1) :+ gm
		
		symeigensystem(P, V=., L=.)
		return(V)
	}

//	Rank of each value among the distinct sorted values (1,2,...)
	real colvector _spsiv_rank(real colvector v)
	{
		o = order(v, 1)
		s = v[o]
		q = (. \ s)
		q = q[|1 \ rows(s)|]
		g = runningsum(s :!= q)
		r = J(rows(v), 1, .)
		r[o] = g
		return(r)
	}

//	Synthetic instruments for balanced panel data
	void _spsiv_core( /*
*/		real matrix E, /*
*/		string scalar vars_str, /*
*/		string scalar outs_str, /*
*/		string scalar touse, /*
*/		string scalar idvar, /*
*/		string scalar tvar, /*
*/		real scalar alpha /*
*/	) {
		vars = tokens(vars_str)
		outs = tokens(outs_str)
		k = cols(vars)
		
		// Position of each selected observation in the N x T layout
		id = st_data(., idvar, touse)
		if (tvar == "") tm = J(rows(id), 1, 1)
		else            tm = st_data(., tvar, touse)
		
		M  = rows(id)
		iu = _spsiv_rank(id)
		it = _spsiv_rank(tm)
		N  = max(iu)
		T  = max(it)
		
		if (M != N*T) {
			errprintf("panel must be balanced in the estimation sample\n")
			exit(498)
		}
		p = (iu :- 1) :* T :+ it
		chk = J(N*T, 1, 0)
		chk[p] = J(M, 1, 1)
		if (sum(chk) != N*T) {
			errprintf("panel must be balanced in the estimation sample\n")
			exit(498)
		}
		if (rows(E) != N) {
			errprintf("number of spatial units in the sample (%g) differs from the dimension of the weighting matrix (%g)\n", N, rows(E))
			exit(498)
		}
		if (N < 3) {
			errprintf("at least 3 spatial units are required\n")
			exit(498)
		}
		
		// Critical value of |r| equivalent to p-value < alpha of the slope
		df    = N - 2
		tc    = invttail(df, alpha/2)
		rcrit = tc / sqrt(df + tc^2)
		
		// Centered and standardized eigenvectors (computed once)
		Ec = E :- mean(E)
		nE = sqrt(colsum(Ec:^2))
		ok = nE :> 1e-8
		Zc = Ec :* (ok :/ (nE :+ !ok))
		
		Xall = st_data(., vars, touse)
		
		for (v = 1; v <= k; v++) {
			
			// Arrange variable v as N x T
			xv = J(N*T, 1, .)
			xv[p] = Xall[., v]
			X  = colshape(xv, T)
			
			xm  = mean(X)
			Xc  = X :- xm
			nx  = sqrt(colsum(Xc:^2))
			okx = (colmax(X) :> colmin(X)) :& (nx :> 0)
			Zx  = Xc :* (okx :/ (nx :+ !okx))
			
			// Correlations of every x_t with every E_i (T x N)
			R = cross(Zx, Zc)
			
			S = J(N, T, .)
			for (t = 1; t <= T; t++) {
				idx = selectindex(abs(R[t, .]) :> rcrit)
				if (length(idx) == 0) {
					S[., t] = J(N, 1, xm[t])
				}
				else {
					Zs = Ec[., idx]
					S[., t] = xm[t] :+ Zs * (invsym(cross(Zs, Zs)) * cross(Zs, Xc[., t]))
				}
			}
			
			// Store to Stata (same observation order as st_data)
			sv = rowshape(S, 1)'
			j  = st_addvar("double", outs[v])
			st_store(., j, touse, sv[p])
		}
	}
end
