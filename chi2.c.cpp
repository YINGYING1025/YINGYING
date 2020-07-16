 ////////////////////////////////////////////////////////////////////////////////
 /// \f$ \chi^{2} \f$ test for comparing weighted and unweighted histograms
 ///
 /// Function: Returns p-value. Other return values are specified by the 3rd parameter
 ///
 /// \param[in] h2 the second histogram
 /// \param[in] option
 ///   - "UU" = experiment experiment comparison (unweighted-unweighted)
 ///   - "UW" = experiment MC comparison (unweighted-weighted). Note that
 ///      the first histogram should be unweighted
 ///   - "WW" = MC MC comparison (weighted-weighted)
 ///   - "NORM" = to be used when one or both of the histograms is scaled
 ///              but the histogram originally was unweighted
 ///   - by default underflows and overflows are not included:
 ///      * "OF" = overflows included
 ///      * "UF" = underflows included
 ///   - "P" = print chi2, ndf, p_value, igood
 ///   - "CHI2" = returns chi2 instead of p-value
 ///   - "CHI2/NDF" = returns \f$ \chi^{2} \f$/ndf
 /// \param[in] res not empty - computes normalized residuals and returns them in this array
 ///
 /// The current implementation is based on the papers \f$ \chi^{2} \f$ test for comparison
 /// of weighted and unweighted histograms" in Proceedings of PHYSTAT05 and
 /// "Comparison weighted and unweighted histograms", arXiv:physics/0605123
 /// by N.Gagunashvili. This function has been implemented by Daniel Haertl in August 2006.
 ///
 /// #### Introduction:
 ///
 /// A frequently used technique in data analysis is the comparison of
 /// histograms. First suggested by Pearson [1] the \f$ \chi^{2} \f$  test of
 /// homogeneity is used widely for comparing usual (unweighted) histograms.
 /// This paper describes the implementation modified \f$ \chi^{2} \f$ tests
 /// for comparison of weighted and unweighted  histograms and two weighted
 /// histograms [2] as well as usual Pearson's \f$ \chi^{2} \f$ test for
 /// comparison two usual (unweighted) histograms.
 ///
 /// #### Overview:
 ///
 /// Comparison of two histograms expect hypotheses that two histograms
 /// represent identical distributions. To make a decision p-value should
 /// be calculated. The hypotheses of identity is rejected if the p-value is
 /// lower then some significance level. Traditionally significance levels
 /// 0.1, 0.05 and 0.01 are used. The comparison procedure should include an
 /// analysis of the residuals which is often helpful in identifying the
 /// bins of histograms responsible for a significant overall \f$ \chi^{2} \f$ value.
 /// Residuals are the difference between bin contents and expected bin
 /// contents. Most convenient for analysis are the normalized residuals. If
 /// hypotheses of identity are valid then normalized residuals are
 /// approximately independent and identically distributed random variables
 /// having N(0,1) distribution. Analysis of residuals expect test of above
 /// mentioned properties of residuals. Notice that indirectly the analysis
 /// of residuals increase the power of \f$ \chi^{2} \f$ test.
 ///
 /// #### Methods of comparison:
 ///
 /// \f$ \chi^{2} \f$ test for comparison two (unweighted) histograms:
 /// Let us consider two  histograms with the  same binning and the  number
 /// of bins equal to r. Let us denote the number of events in the ith bin
 /// in the first histogram as ni and as mi in the second one. The total
 /// number of events in the first histogram is equal to:
 /// \f[
 ///  N = \sum_{i=1}^{r} n_{i}
 /// \f]
 /// and
 /// \f[
 ///  M = \sum_{i=1}^{r} m_{i}
 /// \f]
 /// in the second histogram. The hypothesis of identity (homogeneity) [3]
 /// is that the two histograms represent random values with identical
 /// distributions. It is equivalent that there exist r constants p1,...,pr,
 /// such that
 /// \f[
 ///\sum_{i=1}^{r} p_{i}=1
 /// \f]
 /// and the probability of belonging to the ith bin for some measured value
 /// in both experiments is equal to pi. The number of events in the ith
 /// bin is a random variable with a distribution approximated by a Poisson
 /// probability distribution
 /// \f[
 ///\frac{e^{-Np_{i}}(Np_{i})^{n_{i}}}{n_{i}!}
 /// \f]
 ///for the first histogram and with distribution
 /// \f[
 ///\frac{e^{-Mp_{i}}(Mp_{i})^{m_{i}}}{m_{i}!}
 /// \f]
 /// for the second histogram. If the hypothesis of homogeneity is valid,
 /// then the  maximum likelihood estimator of pi, i=1,...,r, is
 /// \f[
 ///\hat{p}_{i}= \frac{n_{i}+m_{i}}{N+M}
 /// \f]
 /// and then
 /// \f[
 ///  X^{2} = \sum_{i=1}^{r}\frac{(n_{i}-N\hat{p}_{i})^{2}}{N\hat{p}_{i}} + \sum_{i=1}^{r}\frac{(m_{i}-M\hat{p}_{i})^{2}}{M\hat{p}_{i}} =\frac{1}{MN} \sum_{i=1}^{r}\frac{(Mn_{i}-Nm_{i})^{2}}{n_{i}+m_{i}}
 /// \f]
 /// has approximately a \f$ \chi^{2}_{(r-1)} \f$ distribution [3].
 /// The comparison procedure can include an analysis of the residuals which
 /// is often helpful in identifying the bins of histograms responsible for
 /// a significant overall \f$ \chi^{2} \f$ value. Most convenient for
 /// analysis are the adjusted (normalized) residuals [4]
 /// \f[
 ///  r_{i} = \frac{n_{i}-N\hat{p}_{i}}{\sqrt{N\hat{p}_{i}}\sqrt{(1-N/(N+M))(1-(n_{i}+m_{i})/(N+M))}}
 /// \f]
 /// If hypotheses of  homogeneity are valid then residuals ri are
 /// approximately independent and identically distributed random variables
 /// having N(0,1) distribution. The application of the \f$ \chi^{2} \f$ test has
 /// restrictions related to the value of the expected frequencies Npi,
 /// Mpi, i=1,...,r. A conservative rule formulated in [5] is that all the
 /// expectations must be 1 or greater for both histograms. In practical
 /// cases when expected frequencies are not known the estimated expected
 /// frequencies \f$ M\hat{p}_{i}, N\hat{p}_{i}, i=1,...,r \f$ can be used.
 ///
 /// #### Unweighted and weighted histograms comparison:
 ///
 /// A simple modification of the ideas described above can be used for the
 /// comparison of the usual (unweighted) and weighted histograms. Let us
 /// denote the number of events in the ith bin in the unweighted
 /// histogram as ni and the common weight of events in the ith bin of the
 /// weighted histogram as wi. The total number of events in the
 /// unweighted histogram is equal to
 ///\f[
 ///  N = \sum_{i=1}^{r} n_{i}
 ///\f]
 /// and the total weight of events in the weighted histogram is equal to
 ///\f[
 ///  W = \sum_{i=1}^{r} w_{i}
 ///\f]
 /// Let us formulate the hypothesis of identity of an unweighted histogram
 /// to a weighted histogram so that there exist r constants p1,...,pr, such
 /// that
 ///\f[
 ///  \sum_{i=1}^{r} p_{i} = 1
 ///\f]
 /// for the unweighted histogram. The weight wi is a random variable with a
 /// distribution approximated by the normal probability distribution
 /// \f$ N(Wp_{i},\sigma_{i}^{2}) \f$ where \f$ \sigma_{i}^{2} \f$ is the variance of the weight wi.
 /// If we replace the variance \f$ \sigma_{i}^{2} \f$
 /// with estimate \f$ s_{i}^{2} \f$ (sum of squares of weights of
 /// events in the ith bin) and the hypothesis of identity is valid, then the
 /// maximum likelihood estimator of  pi,i=1,...,r, is
 ///\f[
 ///  \hat{p}_{i} = \frac{Ww_{i}-Ns_{i}^{2}+\sqrt{(Ww_{i}-Ns_{i}^{2})^{2}+4W^{2}s_{i}^{2}n_{i}}}{2W^{2}}
 ///\f]
 /// We may then use the test statistic
 ///\f[
 ///  X^{2} = \sum_{i=1}^{r} \frac{(n_{i}-N\hat{p}_{i})^{2}}{N\hat{p}_{i}} + \sum_{i=1}^{r} \frac{(w_{i}-W\hat{p}_{i})^{2}}{s_{i}^{2}}
 ///\f]
 /// and it has approximately a \f$ \sigma^{2}_{(r-1)} \f$ distribution [2]. This test, as well
 /// as the original one [3], has a restriction on the expected frequencies. The
 /// expected frequencies recommended for the weighted histogram is more than 25.
 /// The value of the minimal expected frequency can be decreased down to 10 for
 /// the case when the weights of the events are close to constant. In the case
 /// of a weighted histogram if the number of events is unknown, then we can
 /// apply this recommendation for the equivalent number of events as
 ///\f[
 ///  n_{i}^{equiv} = \frac{ w_{i}^{2} }{ s_{i}^{2} }
 ///\f]
 /// The minimal expected frequency for an unweighted histogram must be 1. Notice
 /// that any usual (unweighted) histogram can be considered as a weighted
 /// histogram with events that have constant weights equal to 1.
 /// The variance \f$ z_{i}^{2} \f$ of the difference between the weight wi
 /// and the estimated expectation value of the weight is approximately equal to:
 ///\f[
 ///  z_{i}^{2} = Var(w_{i}-W\hat{p}_{i}) = N\hat{p}_{i}(1-N\hat{p}_{i})\left(\frac{Ws_{i}^{2}}{\sqrt{(Ns_{i}^{2}-w_{i}W)^{2}+4W^{2}s_{i}^{2}n_{i}}}\right)^{2}+\frac{s_{i}^{2}}{4}\left(1+\frac{Ns_{i}^{2}-w_{i}W}{\sqrt{(Ns_{i}^{2}-w_{i}W)^{2}+4W^{2}s_{i}^{2}n_{i}}}\right)^{2}
 ///\f]
 /// The  residuals
 ///\f[
 ///  r_{i} = \frac{w_{i}-W\hat{p}_{i}}{z_{i}}
 ///\f]
 /// have approximately a normal distribution with mean equal to 0 and standard
 /// deviation  equal to 1.
 ///
 /// #### Two weighted histograms comparison:
 ///
 /// Let us denote the common  weight of events of the ith bin in the first
 /// histogram as w1i and as w2i in the second one. The total weight of events
 /// in the first histogram is equal to
 ///\f[
 ///  W_{1} = \sum_{i=1}^{r} w_{1i}
 ///\f]
 /// and
 ///\f[
 ///  W_{2} = \sum_{i=1}^{r} w_{2i}
 ///\f]
 /// in the second histogram. Let us formulate the hypothesis of identity of
 /// weighted histograms so that there exist r constants p1,...,pr, such that
 ///\f[
 ///  \sum_{i=1}^{r} p_{i} = 1
 ///\f]
 /// and also expectation value of weight w1i equal to W1pi and expectation value
 /// of weight w2i equal to W2pi. Weights in both the histograms are random
 /// variables with distributions which can be approximated by a normal
 /// probability distribution \f$ N(W_{1}p_{i},\sigma_{1i}^{2}) \f$ for the first histogram
 /// and by a distribution \f$ N(W_{2}p_{i},\sigma_{2i}^{2}) \f$ for the second.
 /// Here \f$ \sigma_{1i}^{2} \f$ and \f$ \sigma_{2i}^{2} \f$ are the variances
 /// of w1i and w2i with estimators \f$ s_{1i}^{2} \f$ and \f$ s_{2i}^{2} \f$ respectively.
 /// If the hypothesis of identity is valid, then the maximum likelihood and
 /// Least Square Method estimator of pi,i=1,...,r, is
 ///\f[
 ///  \hat{p}_{i} = \frac{w_{1i}W_{1}/s_{1i}^{2}+w_{2i}W_{2} /s_{2i}^{2}}{W_{1}^{2}/s_{1i}^{2}+W_{2}^{2}/s_{2i}^{2}}
 ///\f]
 /// We may then use the test statistic
 ///\f[
 /// X^{2} = \sum_{i=1}^{r} \frac{(w_{1i}-W_{1}\hat{p}_{i})^{2}}{s_{1i}^{2}} + \sum_{i=1}^{r} \frac{(w_{2i}-W_{2}\hat{p}_{i})^{2}}{s_{2i}^{2}} = \sum_{i=1}^{r} \frac{(W_{1}w_{2i}-W_{2}w_{1i})^{2}}{W_{1}^{2}s_{2i}^{2}+W_{2}^{2}s_{1i}^{2}}
 ///\f]
 /// and it has approximately a \f$ \chi^{2}_{(r-1)} \f$ distribution [2].
 /// The normalized or studentised residuals [6]
 ///\f[
 ///  r_{i} = \frac{w_{1i}-W_{1}\hat{p}_{i}}{s_{1i}\sqrt{1 - \frac{1}{(1+W_{2}^{2}s_{1i}^{2}/W_{1}^{2}s_{2i}^{2})}}}
 ///\f]
 /// have approximately a normal distribution with mean equal to 0 and standard
 /// deviation 1. A recommended minimal expected frequency is equal to 10 for
 /// the proposed test.
 ///
 /// #### Numerical examples:
 ///
 /// The method described herein is now illustrated with an example.
 /// We take a distribution
 ///\f[
 /// \phi(x) = \frac{2}{(x-10)^{2}+1} + \frac{1}{(x-14)^{2}+1}       (1)
 ///\f]
 /// defined on the interval [4,16]. Events distributed according to the formula
 /// (1) are simulated to create the unweighted histogram. Uniformly distributed
 /// events are simulated for the weighted histogram with weights calculated by
 /// formula (1). Each histogram has the same number of bins: 20. Fig.1 shows
 /// the result of comparison of the unweighted histogram with 200 events
 /// (minimal expected frequency equal to one) and the weighted histogram with
 /// 500 events (minimal expected frequency equal to 25)
 /// Begin_Macro
 /// ../../../tutorials/math/chi2test.C
 /// End_Macro
 /// Fig 1. An example of comparison of the unweighted histogram with 200 events
 /// and the weighted histogram with 500 events:
 ///   1. unweighted histogram;
 ///   2. weighted histogram;
 ///   3. normalized residuals plot;
 ///   4. normal Q-Q plot of residuals.
 ///
 /// The value of the test statistic \f$ \chi^{2} \f$ is equal to
 /// 21.09 with p-value equal to 0.33, therefore the hypothesis of identity of
 /// the two histograms can be accepted for 0.05 significant level. The behavior
 /// of the normalized residuals plot (see Fig. 1c) and the normal Q-Q plot
 /// (see Fig. 1d) of residuals are regular and we cannot identify the outliers
 /// or bins with a big influence on \f$ \chi^{2} \f$.
 ///
 /// The second example presents the same two histograms but 17 events was added
 /// to content of bin number 15 in unweighted histogram. Fig.2 shows the result
 /// of comparison of the unweighted histogram with 217 events (minimal expected
 /// frequency equal to one) and the weighted histogram with 500 events (minimal
 /// expected frequency equal to 25)
 /// Begin_Macro
 /// ../../../tutorials/math/chi2test.C(17)
 /// End_Macro
 /// Fig 2. An example of comparison of the unweighted histogram with 217 events
 /// and the weighted histogram with 500 events:
 ///   1. unweighted histogram;
 ///   2. weighted histogram;
 ///   3. normalized residuals plot;
 ///   4. normal Q-Q plot of residuals.
 ///
 /// The value of the test statistic \f$ \chi^{2} \f$ is equal to
 /// 32.33 with p-value equal to 0.029, therefore the hypothesis of identity of
 /// the two histograms is rejected for 0.05 significant level. The behavior of
 /// the normalized residuals plot (see Fig. 2c) and the normal Q-Q plot (see
 /// Fig. 2d) of residuals are not regular and we can identify the outlier or
 /// bin with a big influence on \f$ \chi^{2} \f$.
 ///
 /// #### References:
 ///
 ///  - [1] Pearson, K., 1904. On the Theory of Contingency and Its Relation to
 ///    Association and Normal Correlation. Drapers' Co. Memoirs, Biometric
 ///    Series No. 1, London.
 ///  - [2] Gagunashvili, N., 2006. \f$ \sigma^{2} \f$ test for comparison
 ///    of weighted and unweighted histograms. Statistical Problems in Particle
 ///    Physics, Astrophysics and Cosmology, Proceedings of PHYSTAT05,
 ///    Oxford, UK, 12-15 September 2005, Imperial College Press, London, 43-44.
 ///    Gagunashvili,N., Comparison of weighted and unweighted histograms,
 ///    arXiv:physics/0605123, 2006.
 ///  - [3] Cramer, H., 1946. Mathematical methods of statistics.
 ///    Princeton University Press, Princeton.
 ///  - [4] Haberman, S.J., 1973. The analysis of residuals in cross-classified tables.
 ///    Biometrics 29, 205-220.
 ///  - [5] Lewontin, R.C. and Felsenstein, J., 1965. The robustness of homogeneity
 ///    test in 2xN tables. Biometrics 21, 19-33.
 ///  - [6] Seber, G.A.F., Lee, A.J., 2003, Linear Regression Analysis.
 ///    John Wiley & Sons Inc., New York.
 
 Double_t TH1::Chi2Test(const TH1* h2, Option_t *option, Double_t *res) const
 {
    Double_t chi2 = 0;
    Int_t ndf = 0, igood = 0;
 
    TString opt = option;
    opt.ToUpper();
 
    Double_t prob = Chi2TestX(h2,chi2,ndf,igood,option,res);
 
    if(opt.Contains("P")) {
       printf("Chi2 = %f, Prob = %g, NDF = %d, igood = %d\n", chi2,prob,ndf,igood);
    }
    if(opt.Contains("CHI2/NDF")) {
       if (ndf == 0) return 0;
       return chi2/ndf;
    }
    if(opt.Contains("CHI2")) {
       return chi2;
    }
 
    return prob;
 }
 
 ////////////////////////////////////////////////////////////////////////////////
 /// The computation routine of the Chisquare test. For the method description,
 /// see Chi2Test() function.
 ///
 /// \return p-value
 /// \param[in] h2 the second histogram
 /// \param[in] option
 ///  - "UU" = experiment experiment comparison (unweighted-unweighted)
 ///  - "UW" = experiment MC comparison (unweighted-weighted). Note that the first
 ///        histogram should be unweighted
 ///  - "WW" = MC MC comparison (weighted-weighted)
 ///  - "NORM" = if one or both histograms is scaled
 ///  - "OF" = overflows included
 ///  - "UF" = underflows included
 ///      by default underflows and overflows are not included
 /// \param[out] igood test output
 ///    - igood=0 - no problems
 ///    - For unweighted unweighted  comparison
 ///      - igood=1'There is a bin in the 1st histogram with less than 1 event'
 ///      - igood=2'There is a bin in the 2nd histogram with less than 1 event'
 ///      - igood=3'when the conditions for igood=1 and igood=2 are satisfied'
 ///    - For  unweighted weighted  comparison
 ///      - igood=1'There is a bin in the 1st histogram with less then 1 event'
 ///      - igood=2'There is a bin in the 2nd histogram with less then 10 effective number of events'
 ///      - igood=3'when the conditions for igood=1 and igood=2 are satisfied'
 ///    - For  weighted weighted  comparison
 ///      - igood=1'There is a bin in the 1st  histogram with less then 10 effective
 ///        number of events'
 ///      - igood=2'There is a bin in the 2nd  histogram with less then 10 effective
 ///        number of events'
 ///      - igood=3'when the conditions for igood=1 and igood=2 are satisfied'
 /// \param[out] chi2 chisquare of the test
 /// \param[out] ndf number of degrees of freedom (important, when both histograms have the same empty bins)
 /// \param[out] res normalized residuals for further analysis
 
 Double_t TH1::Chi2TestX(const TH1* h2,  Double_t &chi2, Int_t &ndf, Int_t &igood, Option_t *option,  Double_t *res) const
 {
 
    Int_t i_start, i_end;
    Int_t j_start, j_end;
    Int_t k_start, k_end;
 
    Double_t sum1 = 0.0, sumw1 = 0.0;
    Double_t sum2 = 0.0, sumw2 = 0.0;
 
    chi2 = 0.0;
    ndf = 0;
 
    TString opt = option;
    opt.ToUpper();
 
    if (fBuffer) const_cast<TH1*>(this)->BufferEmpty();
 
    const TAxis *xaxis1 = GetXaxis();
    const TAxis *xaxis2 = h2->GetXaxis();
    const TAxis *yaxis1 = GetYaxis();
    const TAxis *yaxis2 = h2->GetYaxis();
    const TAxis *zaxis1 = GetZaxis();
    const TAxis *zaxis2 = h2->GetZaxis();
 
    Int_t nbinx1 = xaxis1->GetNbins();
    Int_t nbinx2 = xaxis2->GetNbins();
    Int_t nbiny1 = yaxis1->GetNbins();
    Int_t nbiny2 = yaxis2->GetNbins();
    Int_t nbinz1 = zaxis1->GetNbins();
    Int_t nbinz2 = zaxis2->GetNbins();
 
    //check dimensions
    if (this->GetDimension() != h2->GetDimension() ){
       Error("Chi2TestX","Histograms have different dimensions.");
       return 0.0;
    }
 
    //check number of channels
    if (nbinx1 != nbinx2) {
       Error("Chi2TestX","different number of x channels");
    }
    if (nbiny1 != nbiny2) {
       Error("Chi2TestX","different number of y channels");
    }
    if (nbinz1 != nbinz2) {
       Error("Chi2TestX","different number of z channels");
    }
 
    //check for ranges
    i_start = j_start = k_start = 1;
    i_end = nbinx1;
    j_end = nbiny1;
    k_end = nbinz1;
 
    if (xaxis1->TestBit(TAxis::kAxisRange)) {
       i_start = xaxis1->GetFirst();
       i_end   = xaxis1->GetLast();
    }
    if (yaxis1->TestBit(TAxis::kAxisRange)) {
       j_start = yaxis1->GetFirst();
       j_end   = yaxis1->GetLast();
    }
    if (zaxis1->TestBit(TAxis::kAxisRange)) {
       k_start = zaxis1->GetFirst();
       k_end   = zaxis1->GetLast();
    }
 
 
    if (opt.Contains("OF")) {
       if (GetDimension() == 3) k_end = ++nbinz1;
       if (GetDimension() >= 2) j_end = ++nbiny1;
       if (GetDimension() >= 1) i_end = ++nbinx1;
    }
 
    if (opt.Contains("UF")) {
       if (GetDimension() == 3) k_start = 0;
       if (GetDimension() >= 2) j_start = 0;
       if (GetDimension() >= 1) i_start = 0;
    }
 
    ndf = (i_end - i_start + 1) * (j_end - j_start + 1) * (k_end - k_start + 1) - 1;
 
    Bool_t comparisonUU = opt.Contains("UU");
    Bool_t comparisonUW = opt.Contains("UW");
    Bool_t comparisonWW = opt.Contains("WW");
    Bool_t scaledHistogram  = opt.Contains("NORM");
 
    if (scaledHistogram && !comparisonUU) {
       Info("Chi2TestX", "NORM option should be used together with UU option. It is ignored");
    }
 
    // look at histo global bin content and effective entries
    Stat_t s[kNstat];
    GetStats(s);// s[1] sum of squares of weights, s[0] sum of weights
    Double_t sumBinContent1 = s[0];
    Double_t effEntries1 = (s[1] ? s[0] * s[0] / s[1] : 0.0);
 
    h2->GetStats(s);// s[1] sum of squares of weights, s[0] sum of weights
    Double_t sumBinContent2 = s[0];
    Double_t effEntries2 = (s[1] ? s[0] * s[0] / s[1] : 0.0);
 
    if (!comparisonUU && !comparisonUW && !comparisonWW ) {
       // deduce automatically from type of histogram
       if (TMath::Abs(sumBinContent1 - effEntries1) < 1) {
          if ( TMath::Abs(sumBinContent2 - effEntries2) < 1) comparisonUU = true;
          else comparisonUW = true;
       }
       else comparisonWW = true;
    }
    // check unweighted histogram
    if (comparisonUW) {
       if (TMath::Abs(sumBinContent1 - effEntries1) >= 1) {
          Warning("Chi2TestX","First histogram is not unweighted and option UW has been requested");
       }
    }
    if ( (!scaledHistogram && comparisonUU)   ) {
       if ( ( TMath::Abs(sumBinContent1 - effEntries1) >= 1) || (TMath::Abs(sumBinContent2 - effEntries2) >= 1) ) {
          Warning("Chi2TestX","Both histograms are not unweighted and option UU has been requested");
       }
    }
 
 
    //get number of events in histogram
    if (comparisonUU && scaledHistogram) {
       for (Int_t i = i_start; i <= i_end; ++i) {
          for (Int_t j = j_start; j <= j_end; ++j) {
             for (Int_t k = k_start; k <= k_end; ++k) {
 
                Int_t bin = GetBin(i, j, k);
 
                Double_t cnt1 = RetrieveBinContent(bin);
                Double_t cnt2 = h2->RetrieveBinContent(bin);
                Double_t e1sq = GetBinErrorSqUnchecked(bin);
                Double_t e2sq = h2->GetBinErrorSqUnchecked(bin);
 
                if (e1sq > 0.0) cnt1 = TMath::Floor(cnt1 * cnt1 / e1sq + 0.5); // avoid rounding errors
                else cnt1 = 0.0;
 
                if (e2sq > 0.0) cnt2 = TMath::Floor(cnt2 * cnt2 / e2sq + 0.5); // avoid rounding errors
                else cnt2 = 0.0;
 
                // sum contents
                sum1 += cnt1;
                sum2 += cnt2;
                sumw1 += e1sq;
                sumw2 += e2sq;
             }
          }
       }
       if (sumw1 <= 0.0 || sumw2 <= 0.0) {
          Error("Chi2TestX", "Cannot use option NORM when one histogram has all zero errors");
          return 0.0;
       }
 
    } else {
       for (Int_t i = i_start; i <= i_end; ++i) {
          for (Int_t j = j_start; j <= j_end; ++j) {
             for (Int_t k = k_start; k <= k_end; ++k) {
 
                Int_t bin = GetBin(i, j, k);
 
                sum1 += RetrieveBinContent(bin);
                sum2 += h2->RetrieveBinContent(bin);
 
                if ( comparisonWW ) sumw1 += GetBinErrorSqUnchecked(bin);
                if ( comparisonUW || comparisonWW ) sumw2 += h2->GetBinErrorSqUnchecked(bin);
             }
          }
       }
    }
    //checks that the histograms are not empty
    if (sum1 == 0.0 || sum2 == 0.0) {
       Error("Chi2TestX","one histogram is empty");
       return 0.0;
    }
 
    if ( comparisonWW  && ( sumw1 <= 0.0 && sumw2 <= 0.0 ) ){
       Error("Chi2TestX","Hist1 and Hist2 have both all zero errors\n");
       return 0.0;
    }
 
    //THE TEST
    Int_t m = 0, n = 0;
 
    //Experiment - experiment comparison
    if (comparisonUU) {
       Double_t sum = sum1 + sum2;
       for (Int_t i = i_start; i <= i_end; ++i) {
          for (Int_t j = j_start; j <= j_end; ++j) {
             for (Int_t k = k_start; k <= k_end; ++k) {
 
                Int_t bin = GetBin(i, j, k);
 
                Double_t cnt1 = RetrieveBinContent(bin);
                Double_t cnt2 = h2->RetrieveBinContent(bin);
 
                if (scaledHistogram) {
                   // scale bin value to effective bin entries
                   Double_t e1sq = GetBinErrorSqUnchecked(bin);
                   Double_t e2sq = h2->GetBinErrorSqUnchecked(bin);
 
                   if (e1sq > 0) cnt1 = TMath::Floor(cnt1 * cnt1 / e1sq + 0.5); // avoid rounding errors
                   else cnt1 = 0;
 
                   if (e2sq > 0) cnt2 = TMath::Floor(cnt2 * cnt2 / e2sq + 0.5); // avoid rounding errors
                   else cnt2 = 0;
                }
 
                if (Int_t(cnt1) == 0 && Int_t(cnt2) == 0) --ndf;  // no data means one degree of freedom less
                else {
 
                   Double_t cntsum = cnt1 + cnt2;
                   Double_t nexp1 = cntsum * sum1 / sum;
                   //Double_t nexp2 = binsum*sum2/sum;
 
                   if (res) res[i - i_start] = (cnt1 - nexp1) / TMath::Sqrt(nexp1);
 
                   if (cnt1 < 1) ++m;
                   if (cnt2 < 1) ++n;
 
                   //Habermann correction for residuals
                   Double_t correc = (1. - sum1 / sum) * (1. - cntsum / sum);
                   if (res) res[i - i_start] /= TMath::Sqrt(correc);
 
                   Double_t delta = sum2 * cnt1 - sum1 * cnt2;
                   chi2 += delta * delta / cntsum;
                }
             }
          }
       }
       chi2 /= sum1 * sum2;
 
       // flag error only when of the two histogram is zero
       if (m) {
          igood += 1;
          Info("Chi2TestX","There is a bin in h1 with less than 1 event.\n");
       }
       if (n) {
          igood += 2;
          Info("Chi2TestX","There is a bin in h2 with less than 1 event.\n");
       }
 
       Double_t prob = TMath::Prob(chi2,ndf);
       return prob;
 
    }
 
    // unweighted - weighted  comparison
    // case of error = 0 and content not zero is treated without problems by excluding second chi2 sum
    // and can be considered as a data-theory comparison
    if ( comparisonUW ) {
       for (Int_t i = i_start; i <= i_end; ++i) {
          for (Int_t j = j_start; j <= j_end; ++j) {
             for (Int_t k = k_start; k <= k_end; ++k) {
 
                Int_t bin = GetBin(i, j, k);
 
                Double_t cnt1 = RetrieveBinContent(bin);
                Double_t cnt2 = h2->RetrieveBinContent(bin);
                Double_t e2sq = h2->GetBinErrorSqUnchecked(bin);
 
                // case both histogram have zero bin contents
                if (cnt1 * cnt1 == 0 && cnt2 * cnt2 == 0) {
                   --ndf;  //no data means one degree of freedom less
                   continue;
                }
 
                // case weighted histogram has zero bin content and error
                if (cnt2 * cnt2 == 0 && e2sq == 0) {
                   if (sumw2 > 0) {
                      // use as approximated  error as 1 scaled by a scaling ratio
                      // estimated from the total sum weight and sum weight squared
                      e2sq = sumw2 / sum2;
                   }
                   else {
                      // return error because infinite discrepancy here:
                      // bin1 != 0 and bin2 =0 in a histogram with all errors zero
                      Error("Chi2TestX","Hist2 has in bin (%d,%d,%d) zero content and zero errors\n", i, j, k);
                      chi2 = 0; return 0;
                   }
                }
 
                if (cnt1 < 1) m++;
                if (e2sq > 0 && cnt2 * cnt2 / e2sq < 10) n++;
 
                Double_t var1 = sum2 * cnt2 - sum1 * e2sq;
                Double_t var2 = var1 * var1 + 4. * sum2 * sum2 * cnt1 * e2sq;
 
                // if cnt1 is zero and cnt2 = 1 and sum1 = sum2 var1 = 0 && var2 == 0
                // approximate by incrementing cnt1
                // LM (this need to be fixed for numerical errors)
                while (var1 * var1 + cnt1 == 0 || var1 + var2 == 0) {
                   sum1++;
                   cnt1++;
                   var1 = sum2 * cnt2 - sum1 * e2sq;
                   var2 = var1 * var1 + 4. * sum2 * sum2 * cnt1 * e2sq;
                }
                var2 = TMath::Sqrt(var2);
 
                while (var1 + var2 == 0) {
                   sum1++;
                   cnt1++;
                   var1 = sum2 * cnt2 - sum1 * e2sq;
                   var2 = var1 * var1 + 4. * sum2 * sum2 * cnt1 * e2sq;
                   while (var1 * var1 + cnt1 == 0 || var1 + var2 == 0) {
                      sum1++;
                      cnt1++;
                      var1 = sum2 * cnt2 - sum1 * e2sq;
                      var2 = var1 * var1 + 4. * sum2 * sum2 * cnt1 * e2sq;
                   }
                   var2 = TMath::Sqrt(var2);
                }
 
                Double_t probb = (var1 + var2) / (2. * sum2 * sum2);
 
                Double_t nexp1 = probb * sum1;
                Double_t nexp2 = probb * sum2;
 
                Double_t delta1 = cnt1 - nexp1;
                Double_t delta2 = cnt2 - nexp2;
 
                chi2 += delta1 * delta1 / nexp1;
 
                if (e2sq > 0) {
                   chi2 += delta2 * delta2 / e2sq;
                }
 
                if (res) {
                   if (e2sq > 0) {
                      Double_t temp1 = sum2 * e2sq / var2;
                      Double_t temp2 = 1.0 + (sum1 * e2sq - sum2 * cnt2) / var2;
                      temp2 = temp1 * temp1 * sum1 * probb * (1.0 - probb) + temp2 * temp2 * e2sq / 4.0;
                      // invert sign here
                      res[i - i_start] = - delta2 / TMath::Sqrt(temp2);
                   }
                   else
                      res[i - i_start] = delta1 / TMath::Sqrt(nexp1);
                }
             }
          }
       }
 
       if (m) {
          igood += 1;
          Info("Chi2TestX","There is a bin in h1 with less than 1 event.\n");
       }
       if (n) {
          igood += 2;
          Info("Chi2TestX","There is a bin in h2 with less than 10 effective events.\n");
       }
 
       Double_t prob = TMath::Prob(chi2, ndf);
 
       return prob;
    }
 
    // weighted - weighted  comparison
    if (comparisonWW) {
       for (Int_t i = i_start; i <= i_end; ++i) {
          for (Int_t j = j_start; j <= j_end; ++j) {
             for (Int_t k = k_start; k <= k_end; ++k) {
 
                Int_t bin = GetBin(i, j, k);
                Double_t cnt1 = RetrieveBinContent(bin);
                Double_t cnt2 = h2->RetrieveBinContent(bin);
                Double_t e1sq = GetBinErrorSqUnchecked(bin);
                Double_t e2sq = h2->GetBinErrorSqUnchecked(bin);
 
                // case both histogram have zero bin contents
                // (use square of content to avoid numerical errors)
                 if (cnt1 * cnt1 == 0 && cnt2 * cnt2 == 0) {
                    --ndf;  //no data means one degree of freedom less
                    continue;
                 }
 
                 if (e1sq == 0 && e2sq == 0) {
                    // cannot treat case of booth histogram have zero zero errors
                   Error("Chi2TestX","h1 and h2 both have bin %d,%d,%d with all zero errors\n", i,j,k);
                   chi2 = 0; return 0;
                }
 
                Double_t sigma = sum1 * sum1 * e2sq + sum2 * sum2 * e1sq;
                Double_t delta = sum2 * cnt1 - sum1 * cnt2;
                chi2 += delta * delta / sigma;
 
                if (res) {
                   Double_t temp = cnt1 * sum1 * e2sq + cnt2 * sum2 * e1sq;
                   Double_t probb = temp / sigma;
                   Double_t z = 0;
                   if (e1sq > e2sq) {
                      Double_t d1 = cnt1 - sum1 * probb;
                      Double_t s1 = e1sq * ( 1. - e2sq * sum1 * sum1 / sigma );
                      z = d1 / TMath::Sqrt(s1);
                   }
                   else {
                      Double_t d2 = cnt2 - sum2 * probb;
                      Double_t s2 = e2sq * ( 1. - e1sq * sum2 * sum2 / sigma );
                      z = -d2 / TMath::Sqrt(s2);
                   }
                   res[i - i_start] = z;
                }
 
                if (e1sq > 0 && cnt1 * cnt1 / e1sq < 10) m++;
                if (e2sq > 0 && cnt2 * cnt2 / e2sq < 10) n++;
             }
          }
       }
       if (m) {
          igood += 1;
          Info("Chi2TestX","There is a bin in h1 with less than 10 effective events.\n");
       }
       if (n) {
          igood += 2;
          Info("Chi2TestX","There is a bin in h2 with less than 10 effective events.\n");
       }
       Double_t prob = TMath::Prob(chi2, ndf);
       return prob;
    }
    return 0;
 }
 ////////////////////////////////////////////////////////////////////////////////
 /// Compute and return the chisquare of this histogram with respect to a function
 /// The chisquare is computed by weighting each histogram point by the bin error
 /// By default the full range of the histogram is used.
 /// Use option "R" for restricting the chisquare calculation to the given range of the function
 /// Use option "L" for using the chisquare based on the poisson likelihood (Baker-Cousins Chisquare)
 
 Double_t TH1::Chisquare(TF1 * func, Option_t *option) const
 {
    if (!func) {
       Error("Chisquare","Function pointer is Null - return -1");
       return -1;
    }
 
    TString opt(option); opt.ToUpper();
    bool useRange = opt.Contains("R");
    bool usePL = opt.Contains("L");
 
    return ROOT::Fit::Chisquare(*this, *func, useRange, usePL);
 }
