********************************************************************************
********************************************************************************
********************************************************************************
*
* Demo Stata 16 svy: proportion shortcomings
* and the programmer's command, svypd, that was written 
* to overcome them.
*
* May 25, 2020
*
* Dale Rhoda
* Biostat Global Consulting
* Dale.Rhoda@biostatglobal.com
*
*
********************************************************************************
* Revision history
*
*
* 2020-05-25    Dale Rhoda      Original version
* 2020-05-26	Dale Rhoda		Update comments re: ICC for y50b, y50c, y75b 
*                               and y75c.  Where
*                               coverage is precisely the same in every
*                               cluster; ICC there is -1 / (n0 - 1)
*                               where n0 is the weighted average number of
*                               respondents per cluster
*
*                               Note that the organ pipe plot PowerPoint
*                               from the 2018 Stata conference has this 
*                               same small error; it says that the ICC is
*                               zero when coverage is equal across clusters
*                               and in a sense we treat it as zero because we
*                               usually do not allow negative values when
*                               we analyze vaccination coverage surveys, so
*                               if the observed value is < 0, we set it to 0.
*
*                               But technically, the value is -1 / (n0-1).
*                               I will update the Powerpoint in the organ-pipe 
*                               GitHub repositories in the coming days to 
*                               correct this error.
* 2020-05-28	Dale Rhoda      Added this revision history.
*
********************************************************************************
********************************************************************************

* Generate a small svy dataset with several challenging, but plausible outcomes
*
* The dataset has fifty PSUs with four respondents each

clear
set seed 8675309

set obs 200
gen psu = ceil(_n / 4)
tab psu ,m

bysort psu : gen respid = _n

* For the outcome y0, no one has the outcome of interest
gen y0 = 0


* For the 3 outcomes named y50a y50b and y50c, half of the respondents
* have the outcome of interest

* For y50a, everyone in half the clusters has it
gen y50a = _n <= 100

* For y50b, half the respondents in every cluster has it
gen y50b = respid < 3

* For y50c, a random half of respondents have the outcome; the 
* outcome is independent of location

gen rand = runiform()
sort rand
gen y50c = _n <= 100
drop rand

* For y50d, half of the respondents have the outcome of interest, but
* there is some positive intracluster correlation 
gen rand = runiform()*50 + psu
sort rand
gen y50d = _n <= 100
drop rand

sort psu respid

gen y75a = _n <= 150

gen y75b = respid < 4

gen rand = runiform()
sort rand
gen y75c = _n <= 150
drop rand

gen rand = runiform()*50 + psu
sort rand
gen y75d = _n <= 150
drop rand

gen y99 = !inlist(_n,1,200)

gen y01 = !y99

* For the outcome y100, everyone has it
gen y100 = 1

* Re-sort dataset into psu and respondent id order
sort psu respid

* Organ pipe plots show positive outcomes as a colored portion of a vertical
* bar: one bar per cluster.  See references below.

********************************************************************************
* For the opplot syntax to run, you'll need to put the .ado files that 
* accompany this file somewhere in your Stata adopath.  This can 
* be accompished in several ways:
* 
* Type the command: adopath + "<full path to folder where you save the .ado files>"
*
* Or save the .ado files in what Stata considers to be your PERSONAL folder.
* (Type the command "adopath" to learn where your PERSONAL folder is.)
********************************************************************************

* The distribution of respondents with the outcome across clusters varies 
* markedly here
opplot y0  , clustvar(psu) xtitle(Primary Sampling Units (PSUs)) title(y0)   name(y0  , replace) export(y0.png)  exportwidth(2000) //   0%
opplot y01 , clustvar(psu) xtitle(Primary Sampling Units (PSUs)) title(y01)  name(y01 , replace) export(y01.png) exportwidth(2000) //   1%
                                                                                                 
opplot y50a, clustvar(psu) xtitle(Primary Sampling Units (PSUs)) title(y50a) name(y50a, replace) export(y50a.png) exportwidth(2000) // 50%; icc is 1; psu COMPLETELY determines outcome
opplot y50b, clustvar(psu) xtitle(Primary Sampling Units (PSUs)) title(y50b) name(y50b, replace) export(y50b.png) exportwidth(2000) // 50%; icc is -1/3; psu is independent of outcome
opplot y50c, clustvar(psu) xtitle(Primary Sampling Units (PSUs)) title(y50c) name(y50c, replace) export(y50c.png) exportwidth(2000) // 50%; icc is very near zero
opplot y50d, clustvar(psu) xtitle(Primary Sampling Units (PSUs)) title(y50d) name(y50d, replace) export(y50d.png) exportwidth(2000) // 50%; icc is ~0.30, so outcome is spatially heterogeneous
                                                                                                 
opplot y75a, clustvar(psu) xtitle(Primary Sampling Units (PSUs)) title(y75a) name(y75a, replace) export(y75a.png) exportwidth(2000) // 75%; icc is nearly 1; psu NEARLY determines outcome
opplot y75b, clustvar(psu) xtitle(Primary Sampling Units (PSUs)) title(y75b) name(y75b, replace) export(y75b.png) exportwidth(2000) // 75%; icc is -1/3; psu is independent of outcome
opplot y75c, clustvar(psu) xtitle(Primary Sampling Units (PSUs)) title(y75c) name(y75c, replace) export(y75c.png) exportwidth(2000) // 75%; icc is very near zero
opplot y75d, clustvar(psu) xtitle(Primary Sampling Units (PSUs)) title(y75d) name(y75d, replace) export(y75d.png) exportwidth(2000) // 75%; icc is ~-.26, so outcome is spatially heterogeneous
                                                                                                
opplot y99,  clustvar(psu) xtitle(Primary Sampling Units (PSUs)) title(y99)  name(y99 , replace) export(y99.png) exportwidth(2000) // 99%
                                                                                                 
opplot y100, clustvar(psu) xtitle(Primary Sampling Units (PSUs)) title(y100) name(y100, replace) export(y100.png) exportwidth(2000) // 100%

label variable y0     "No one has the outcome"
label variable y01    "1% have the outcome"
label variable y50a   "50% - everyone in half the PSUs"
label variable y50b   "50% - half of everyone in all PSUs"
label variable y50c   "50% - everyone in half the PSUs: ICC = -1/3"
label variable y50d   "50% - 0, 1, 2, 3 or 4 per PSU: ICC = 0.30"
label variable y75a   "75% - all in 37 PSUs and half in one PSU"
label variable y75b   "75% - 3/4 in all PSUs"
label variable y75c   "75% - everyone in 3/4 the PSUs: ICC = -1/3"
label variable y75d   "75% - 0, 1, 2, 3, or 4 per PSU: ICC = 0.26"
label variable y99    "99% have the outcome"
label variable y100   "Everyone has the outcome"

label variable psu    "PSU ID"
label variable respid "Respondent ID within the PSU"

* This is a self-weighted cluster sample survey; specify the svyset
svyset psu

compress
aorder

* Save the dataset
save survey_with_plausible_yet_perverse_proportions, replace

********************************************************************************
********************************************************************************
********************************************************************************

use survey_with_plausible_yet_perverse_proportions, clear

describe

* Use Stata's proportion command and svy: proportion command to 
* summarize coverage for each outcome.


foreach outcome in y0 y01 y50a y50b y50c y75a y75b y75c y99 y100 {
	proportion `outcome'
	svy: proportion `outcome'
}

* Note that when the Std. Err. is 0 because every PSU has
* the same sample proportion, svy: proportion does not give a meaningful
* 95% confidence interval.

* Note also that where the Std. Err. is non-zero, the 2-sided
* confidence interval (CI) for the svy: proportion command is appropriately
* wider than that for the simple proportion command, because the svy: 
* CI calculation is based on 49 degrees of freedom (df) (number of clusters
* minus the number of strata) whereas the simple
* random sample calculation from the proportion command is based on 199 df 
* (or N-1).

* But the primary thing I want you to see here is that svy: proportion 
* fails to give a meaningful confidence interval for y0, y50b, y75b and y100.

* Some survey literature (Korn & Graubard, 1998; Dean & Pagano, 2015)
* recommends ignoring the clustering when all or none of the respondents 
* have the outcome of interest or when all the clusters have the same
* observed proportion.
*
* In these situations, there is no cluster effect, so we might do the 
* calculation ignoring the sample design, and for example, a meaningful 
* Clopper-Pearson interval is calculable.  It is probably better to 
* calculate such an interval than to return a null or nonsensical interval.  
* It would also be possible to calculate a Wilson, Agresti, Jeffreys
* or Fleiss interval.
*

ci proportions y0
ci proportions y50b
ci proportions y75b
ci proportions y100

* In every case, that command yields a meaningful confidence interval!  You 
* can also give it an option for agresti or jeffreys or wilson.

* So what we want is to be able to 'trap' the condition where coverage in 
* every PSU is exactly equal and run the ci proportions command instead
* of the svy: proportion command.  Conceptually, that is what the 
* Biostat Global Consulting command svypd.ado does.

* It turns out that in addition to the Clopper-Pearson, there are several 
* confidence interval formulae that provide meaningful results when the cluster  
* level coverage does not vary in the sample.  The Wilson, and Jeffreys and 
* Agresti-Coull intervals, and what I call the Fleiss, which is a Wilson 
* interval with a continuity correction all yield meaningful intervals even if 
* the standard error is zero.

* As far as I can see, the Stata svy: proportion syntax always returns either 
* no interval or a degenerate interval of 0 width when the outcome is
* constant across clusters.

svy: proportion y50a, citype(jeffreys)
svy: proportion y50b, citype(jeffreys)
svy: proportion y50c, citype(jeffreys)


svy: proportion y50a, citype(agresti)
svy: proportion y50b, citype(agresti)
svy: proportion y50c, citype(agresti)

* (By the way, note how wide the CIs are for y50a versus y50c and for y75a 
*  versus y75c.  That is because the design effect takes its largest value
*  when ALL respondents have the outcome in some PSUs and NO respondents
*  have the outcome in the other clusters.  The intracluster correlation 
*  coefficient (ICC) takes its max value of 1 in that situation, which 
*  yields a design effect equal to the average number of respondents per
*  cluster.  In the case of this toy dataset, that is 4.  
*
* Conceptually, in the situation where every cluster has the same sample
* coverage, the design effect takes its minimum value because the ICC is a 
* small negative number (-1 / (n0 - 1)).  But the Stata syntax is not going 
* to show us that.  Stata is tripped up by the fact that the standard error 
* of cluster-sample coverage is 0.

* In our work we often write programs that estimate hundreds or thousands of
* survey proportions
* in a single run, so we can't stop to do special calculations by hand in the
* circumstance where sample coverage is 0% or 100% (or some other value that
* is equal across clusters) so we've written a wrapper program named svypd.ado 
* that does the following:

* 1. Detect right away if the standard error is 0 or missing, and 
*    in that case, calculate a Clopper-Pearson CI if the user has 
*    asked for a Wald or logit interval, that would not
*    be defined because of the 0 or missing std err.
*
* 2. If the user has requested a calculable interval (Wilson, Jeffreys,
*     Agresti-Coull, Clopper-Pearson, or Fleiss) and the standard error
*     is zero, reset the degrees of freedom to be N-1 instead of the
*     survey value of Nclusters - Nstrata.
*
* 3. If StdErr > 0, adjust the estimation to use the survey design degrees of 
*    freedom, if requested.
*
* 4. Truncate the estimate of the design effect to be no smaller than 1.0
*    if requested
*
* 5. Calculate a CI with a LEVEL anywhere between 0.01 and 99.99%.  (Stata
*    imposes an arbitrary miniimum of 10.0 for this parameter without
*    stating or defending a reason why.)
*
* 6. Allow simultaneous calculation of numerous confidence intervals (i.e., 
*    the 99% CI, 98% CI, 97% CI, ... 3% CI, 2% CI, 1% CI) for plotting what
*    we call inchworm plots.  Those are described in the 2018 WHO reference
*    manual cited in the references below and in the 2016 Stata Conference
*    presentation, too, but are outside the scope of the points we are 
*    making in this program.
*
* The svypd program may only be run on a dataset that has been svyset and it
* may only be run on an outcome that takes the values 0, 1 or missing (.).
*
* It bases most of its calculations on formulae from the very clear Dean & 
* Pagano 2015 paper.  It adds a Fleiss interval from Fleiss Levin and Park, 
* 2003.
*
* It returns the survey-weighted proportion estimate, its standard error, and 
* (by default) bounds of a 2-sided 95% CI and a 2-sided 90% CI (which may be 
* taken as 1-sided upper- and lower- 95% bounds).  The user may specify 
* other confidence levels.
*
* It also returns the design effect, the number of degrees of freedom,
* the number of clusters, the weighted and unweighted sample size and the
* name of the type of confidence interval that was calculated.
*
* Again, if the user asks for Wald or logit interval that cannot be 
* calculated because of lack of variation in the cluster-level outcomes, 
* the program substitutes a Clopper-Pearson interval, and returns the 
* string "Clopper-Pearson" in r(method) instead of the method that the 
* user requested.

* If the user requests one of the intervals that is calculable, 
* svypd will return that type of interval.

********************************************************************************
********************************************************************************
********************************************************************************

* Demonstrate several points

********************************************************************************

* 1. The Wald interval gives nonsensical interval bounds (outside the 0-100%)
*    range when the sample proportion is near 0% or 100%.

*    The Wald lower bound can be smaller than 0%.

foreach t in wald logit wilson jeffreys exact agresti {
	svy: proportion y01, citype(`t')
}

*    Or the Wald upper bound can be higher than 100%.
foreach t in wald logit wilson jeffreys exact agresti {
	svy: proportion y99, citype(`t')
}

****************************************

* 2. svy: proportion fails to give a helpful interval when 
*    cluster level coverage does not vary

*    Or the upper bound can be higher than 100%.
foreach t in wald logit wilson jeffreys exact agresti {
	svy: proportion y50b, citype(`t')
}

foreach t in wald logit wilson jeffreys exact agresti {
	svy: proportion y75b, citype(`t')
}


* 3. svypd returns meaningful intervals in that same situation,
*    substituting Clopper-Pearson (exact) if the user requests
*    Wald or Wilson or logit, and calculating what the user
*    requests if they ask for Jeffreys or Fleiss or Agresti-Coull or 
*    Clopper-Pearson.
*
*    Note that svypd is a programmer's command that returns its
*    results in scalars and macros & matrices.  It does not send any output
*    to the output window so we use Stata's "return list" commmand to
*    see the output.

foreach t in wald logit wilson jeffreys exact agresti fleiss {
	di _n "Showing `t' intervals for the problematic proportions:"
	foreach y in y0 y50b y75b y100 {
		di _n "Summarize `y'"
		di "Summarize `y'"
		svypd `y', adjust truncate method(`t') level(95)
		return list
	}
}

* Note that in every case demostrated above, svypd recognized that 
* it was okay to use the full degrees of freedom (N-1) to calculate 
* a narrow interval, ignoring the survey sample design, because of
* a lack of cluster effect.

* Note also that when we use svypd to estimate proportions where the
* cluster-level outcomes vary, it uses the appropriate survey degrees
* of freedom (Nclusters - Nstrata)

foreach t in wald logit wilson jeffreys exact agresti fleiss {
	di _n "Showing `t' intervals for the problematic proportions:"
	foreach y in y01 y50a y50c y50d y75a y75c y75d y99  {
		di _n "Summarize `y'"
		svypd `y', adjust truncate method(`t') level(95)
		return list
	}
}

* Note that for each of those calculations, the df = 49 instead of 199.
*
*
*
*
* Note that svypd does not prevent the user from asking for a Wald interval,
* even in situations where the Wald returns a nonsensical upper or lower-bound.
*
********************************************************************************
********************************************************************************
********************************************************************************

* Recommendation:
*
* In any situation where you are estimating a proportion from a complex 
* survey sample, you can use our svypd command instead of svy: proportion
* and you can have peace of mind knowing that it will return a) the 
* same weighted proportion that svy: proportion would return, and b) a 
* meaningful confidence interval, no matter what the sample proportion
* and no matter how the outcomes are distributed across the primary 
* sampling units.
*
* A number of useful parameters are returned in scalars, macros and
* in the matrix of confidence interval bounds.
*
* This svypd command lies at the heart of the suite of Stata programs we have 
* written for the World Health Organization and it is also the computational
* kernel for our iwplot.ado program that makes inchworm plots to help 
* visualize estimated proportions.
*
* We would be happy to answer questions or receive feedback.
*
* We will continue to rely on svypd because of needing values of level smaller
* than 10 and because of needing vectors of intervals for inchworm plots, 
* but we recommend that Stata update the svy: proportion command to return
* meaningful intervals (Wilson, Jeffreys, Agresti, and Clopper-Pearson) when 
* the StdErr == 0 (or is infinitesimal).  This will help a variety of users 
* obtain a meaningful interval in the not uncommon situation where the 
* sample proportion is 0% or 100%.
* 
* -Dale Rhoda
*  Dale.Rhoda@biostatglobal.com
*
*
********************************************************************************
********************************************************************************
********************************************************************************

* Dean, Natalie, and Marcello Pagano. 2015. “Evaluating Confidence Interval 
*    Methods for Binomial Proportions in Clustered Surveys.” Journal of 
*    Survey Statistics and Methodology 3 (4): 484–503. 
*    https://doi.org/10.1093/jssam/smv024.

* Fleiss, Joseph, Bruce Levin, and Myunghee Cho Paik. 2003. Statistical 
*    Methods for Rates and Proportions. 3rd Ed. 3rd ed. New York: Wiley.

* Franco, Carolina, Roderick JA Little, Thomas A Louis, and Eric V Slud. 
*    2019. “Comparative Study of Confidence Intervals for Proportions in 
*    Complex Sample Surveys.” Journal of Survey Statistics and Methodology 
*    7 (3): 334–64.

* Korn, Edward L., and Barry I. Graubard. 1998. “Confidence Intervals for 
*    Proportions with Small Expected Number of Positive Counts Estimated 
*    from Survey Data.” Survey Methodology 24: 193–201.

* Prier, Mary, and Dale A. Rhoda. 2018. “Organ Pipe Plots for Clustered 
*    Datasets - Visualize Disparities in Cluster-Level Coverage.” 
*    Presented at the Stata Conference 2018, Columbus, Ohio, July 19. 
*    https://github.com/BiostatGlobalConsulting/organ-pipe-plots/blob/master/opplot_presentation.pptx.

* Rhoda, Dale A. 2022 svypd-stata-pubic GitHub repository holding the program
*    that was developed as part of the World Health Organization's Vaccination
*    Coverage Quality Indicators (VCQI) suite of programs.
*    https://github.com/BiostatGlobalConsulting/svypd-stata-public

* Rhoda, Dale A. 2016. “Inchworm Plots: Visual Representation of Inferential 
*    Uncertainty.” Presented at the Stata Conference 2016, Chicago, July. 
*    https://www.stata.com/meeting/chicago16/slides/chicago16_rhoda.pptx.

* Rhoda, Dale A., and Yulia Fungard. 2019. Stata Code to Produce Inchworm 
*    Plots of Estimated Survey Proportions (version 1.25). Stata. VCQI. 
*    Biostat Global Consulting. 
*    https://github.com/BiostatGlobalConsulting/inchworm-plots-stata.

* Rhoda, Dale A., and Mary L. Prier. 2019. OPPLOT: Stata Module to Generate 
*    a Vertical Bar Chart to Summarize a Binary Outcome in Cluster Survey 
*    Data (version 1.13). Stata. Stata. VCQI. Biostat Global Consulting. 
*    https://ideas.repec.org/c/boc/bocode/s458627.html.

* World Health Organization. 2018. “Vaccination Coverage Cluster Surveys: 
*    Reference Manual.” World Health Organization, Geneva. 
*    https://apps.who.int/iris/handle/10665/272820
