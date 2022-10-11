# -*- coding: utf-8 -*-
"""
*** Distribution Fitting Toolbox ***

DESCRIPTION:
    Takes plain-text input file containing a single integer on each line and
    fits the dicrete uniform, beta binomial and zipfian distributions to data
    distribution. 
    
    A summary of the goodness-of-fit metrics (chi-square, R-square, RMSE, 
    Kolmogorov-Smirnov) are displayed for each fit, along with general
    statistics of the input integer data and plot of the distribution and fits,
    along with teh fitting parameters for each distribution.

INPUTS:
    filename ==> string of filename (and path) containing integer data
    printResults ==> Boolean value (True/False), True prints results in console
    plotResults ==> Boolean vale (True/False), True plots results
    
    
OUTPUTS:
    [stats,dists,fit] = distributionFitting(filename,printResults,plotResults)
    stats ==> 1x5 pandas data frame containing mean, standard deviation, 
              range, variance and inter-quartile range of integer data
    dists ==> 4x3 pandas data frame containing chi-square, R-square, RMSE 
              and Kolmogorov-Smirnov metrics, for fittings integer
              distributions to the Discrete Uniform Distribution (DUD),
              Beta Binomial Distribution (BBD) and Zipfian Distribution (ZD)
    fit ==>  1x3 fitting parameters a,b,c from BBD and ZD, defined in METHODS
              
METHODS:
    * Find frequency distribution of integer data [fr]
    * Convert to probability distribution [p]
    * Fit DUD to probability distribution
        - see https://mathworld.wolfram.com/DiscreteUniformDistribution.html
        - Each integer value has uniform probability of 1/n
        - n is number of integers in range of input data
    * Fit BBD to probability distribution
        - see https://mathworld.wolfram.com/BetaBinomialDistribution.html
        - PDF: P(x) = nCx(n,x)*B(x+a,n-x-b)/B(a,b)
        - B ==> beta function
        - nCx(n,x) ==> binomial coefficient of n choose x
        - n ==> len(x)
        - a,b ==> fitting parameters
    * Fit ZD to probability distribution
        - see https://mathworld.wolfram.com/ZipfDistribution.html
        - PDF: P(x) = x^(-c) / z(c)
        - z ==> Riemann zeta function
        - c ==> fitting parameter
    
GOODNESS-OF-FIT METRICS:
    Given observed data, O, and expected data (fit), E
    * chi-square
        - Used to test if a sample of data came from a specific ditribution
        - X^2 = sum_{i=1}^N (O_i - E_i)^2 / E_i
        - small chi = more likely to come from specific distribution
    * R-square
        - Statistical measure of fit, indicate how much of variation of a
          dependent variable is explined by the indeendent variable
        - R^2 = 1 - (sum_{i=1}^N (O_i - E_i)^2)/(sum_{i=1}^N (O_i - mean(O))^2)
        - R^2 = 1 ==> all variation of y-axis is account for by x-axis
    *  Root mean square error (RMSE)
        - Measure of distance between observed and predicted model
        - RMSE = sqrt(1/N * sum_{i = 1}^N (O_i - E_i)^2)
        - small RMSE = smaller distance between observed and fitted data
    * Kolmogorov-Smirnov (K-S)
        - Used with a sample from a population of an unknown distribution
        - FO ==> cumulative distribution of O
        - FE ==> cumulative distribution of E
        - Test statistic, D, gives greatest vertical distance between FO and FE
        - D = max_{i=1}^N | FO_i - FE_i|

LIBRARIES (+ version used to write script):
    numpy       (version 1.21.5)    
    scipy       (version 1.7.3)
    matplotlib  (version 3.5.1)
    pandas      (version 1.4.2)
    tabulate    (version 0.8.9)
    [python      (version 3.9.12)]
    
Author: Stella L. Harrison
"""
def distributionFitting(filename,printResults,plotResults):
    ## Imports
    import numpy as np
    from scipy.optimize import curve_fit
    from scipy.special import beta,binom,zetac
    from scipy.stats import iqr
    import matplotlib.pyplot as plt
    import pandas as pd
    from tabulate import tabulate
    
    #--------------------------------------------------------------------------
    ## Open plain-text file containing integers (1 number per line)
    f = open(filename,'r')
    f2 = f.read()
    RR = f2.splitlines()
    R = list(map(float,RR));
    
    #--------------------------------------------------------------------------
    ## Extract frequency and probability of each integer occuring
    xs = list(np.arange(np.min(R),np.max(R)+1)) # unique integer vaulues in data range
    fr = list([0]*len(xs)) # frequency (count) of each integer
    
    for i in range(0,len(R)):
        fr[xs.index(R[i])] += 1
        
    N = len(R) # total number of integers in file
    n = len(xs) # total number of unique integers in data range
    p = [x/N for x in fr] # convert fr (counts) to p (probability)
    
    #--------------------------------------------------------------------------
    ## Summary of integer statistics
    
    Rmean = np.mean(R)
    Rstd = np.std(R)
    Rrange = np.max(R) - np.min(R)
    Rvar = np.var(R)
    Riqr = iqr(R)
    
    basicStats = [Rmean, Rstd, Rrange, Rvar, Riqr]
    #--------------------------------------------------------------------------
    ## Distribution functions
    
    ## Discrete Uniform Distribution
    #  PDF: P(x) = 1/n
    #  n ==> len(x)
    #  All x values should have probability of 1/n, where n = numel(x)
    def uniDist(x, n):
        y = [1/n]*len(x)
        return y
    
    
    ## Beta Binomial Distribution (BBD)
    #  PDF: P(x) = nCx(n,x)*B(x+a,n-x-b)/B(a,b)
    #  B ==> beta function
    #  nCx(n,x) ==> binomial coefficient of n choose x
    #  n ==> len(x)
    #  a,b ==> fitting parameters
    def betaBinom(x,a,b):
        x = np.array(x)
        n = len(x)
        y = binom(n,x)*beta(a+x,b+n-x)/beta(a,b)
        return y
        
    def betaBinomN(x,n,a,b):
        x = np.array(x)
        y = binom(n,x)*beta(a+x,b+n-x)/beta(a,b)
        return y
    
    ## Zipfian Distributon (ZD)
    #  PDF: P(x) = x^(-c) / z(c)
    #  z ==> Riemann zeta function
    #  c ==> fitting parameter
    def zipfDist(x,c):
        y = [xx**(-1*(c)) for xx in x]/zetac(c)
        return y
     
    #--------------------------------------------------------------------------
    ## Distribution Fitting
    
    # DUD
    # no fit needed as y = 1/n
    
    # BBD
    try:
        b1,b2 = curve_fit(betaBinom,xs,p,bounds=((0,0),(1000,1000)))
    except: # if unable to optimise fit to Beta Binomial Distribution
        b1 = [np.nan]*2
        print("Could not optimize Beta Binomial fit")
        pass
        
    
    # ZD
    try:
        z1,z2 = curve_fit(zipfDist,xs,p)
    except: # if unable to optimise fit to Zipfian Distribution
        z1 = [np.nan]
        print("Could not optimize Zipfian fit")
        pass
    
    fittingParams = [b1[0],b1[1],z1[0]] # returned by function (see end)
    #--------------------------------------------------------------------------
    ## Goodness-of-fit metrics
    
    def chiSquare(O,E):
        counts = len(O)
        chi_sq = 0
        for i in range(0,counts):
            x = (O[i] - E[i]) ** 2
            x = x / E[i]
            chi_sq += x
        return chi_sq 
    
    
    def RSquare(O,E):
        SSres = [(O[i] - E[i])**2 for i in range(0,len(O))]
        SStot = [(O[i] - np.mean(O))**2 for i in range(0,len(O))]
        RR = 1 - (np.sum(SSres)/np.sum(SStot))
        return RR
        
    def RMSE(O,E):
        SE = [(O[i] - E[i])**2 for i in range(0,len(O))]
        RM = np.sqrt(np.mean(SE))
        return RM
    
    def KS(O,E):
        cumO = np.cumsum(O)
        cumE = np.cumsum(E)
        d = [abs(cumO[i] - cumE[i]) for i in range(0,len(O))]
        D = np.max(d)
        return D
                             
    
    # Append each goodness-of-fit parameter for each distribution
    chis = []
    chis.append(chiSquare(p, uniDist(xs, n)))
    chis.append(chiSquare(p, betaBinom(xs, b1[0],b1[1])))
    chis.append(chiSquare(p, zipfDist(xs, z1[0])))
    
    Rs = []
    Rs.append(RSquare(p, uniDist(xs, n)))
    Rs.append(RSquare(p, betaBinom(xs, b1[0],b1[1])))
    Rs.append(RSquare(p, zipfDist(xs, z1[0])))
    
    RMSEs = []
    RMSEs.append(RMSE(p, uniDist(xs, n)))
    RMSEs.append(RMSE(p, betaBinom(xs, b1[0],b1[1])))
    RMSEs.append(RMSE(p, zipfDist(xs, z1[0])))
    
    KSs = []
    KSs.append(KS(p, uniDist(xs, n)))
    KSs.append(KS(p, betaBinom(xs, b1[0],b1[1])))
    KSs.append(KS(p, zipfDist(xs, z1[0])))
    
    # unrounded
    cell_data = [chis, Rs, RMSEs, KSs]
    
    # Round to 3dp for table
    chi_text = ['%1.3f' %(x) for x in chis]
    R_text = ['%1.3f' %(x) for x in Rs]
    RMSE_text = ['%1.3f' %(x) for x in RMSEs]
    KS_text = ['%1.3f' %(x) for x in KSs]

    # round to 3dp for table
    cell_text = [['%1.3f' %(x) for x in cell_data[j]] for j in range(0,np.shape(cell_data)[0])]
    
    #--------------------------------------------------------------------------
    # Display Data as pandas data frame
    
    
    # Make pandas data frame of basic stats data, to 3dp    
    df1 = pd.DataFrame({"Mean":['%1.3f' %(Rmean)], 
                        "Std":['%1.3f' %(Rstd)], 
                        "Range":['%1.3f' %(Rrange)],
                        "Variance":['%1.3f' %(Rvar)], 
                        "IQR":['%1.3f' %(Riqr)]})
    
  
    # Make pandas data frame of goodness-of-fit metrics, to 3dp  
    df2 = pd.DataFrame({"Fit Metric":["Chi-Square","R-Square","RMSE","K-S"],
                       "Discrete Uniform":[chi_text[0],R_text[0],RMSE_text[0],KS_text[0]],
                       "Beta Binomial":[chi_text[1],R_text[1],RMSE_text[1],KS_text[1]],
                       "Zipfian":[chi_text[2],R_text[2],RMSE_text[2],KS_text[2]]})
    
    # print if printResults == True
    if printResults == True: 
        print("\nSummary statistics of integers:")
        print(tabulate(df1, headers = 'keys', tablefmt = 'psql',showindex=False))
        print("\nGooness-of-fit metrics from fitting distributions:")
        print(tabulate(df2, headers = 'keys', tablefmt = 'psql',showindex=False))
        
        print("\n*Fitting parameters*")
        print("\nBeta Binomial Distribution:\na = %1.3f, b = %1.3f" %(b1[0],b1[1]))
        print("\nZipfian Distribution:\nc = %1.3f" %(z1[0]))
                                      
    #--------------------------------------------------------------------------
    # Plot results
    
    # plot if plotResults == True
    if plotResults == True:
        f, (a0, a1) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4, 1]})
        f.set_size_inches(9,9)
        
        
        a0.bar(xs, p, color = 'powderblue', width=0.8)
        a0.set_xlabel('Integer Values',fontsize=12)
        a0.set_ylabel('Probability',fontsize=12)
        a0.set_xticks(xs)
        
        xfine = np.arange(np.min(R), np.max(R), 0.01)
        p1 = a0.plot(xfine, uniDist(xfine,n),linewidth=2,color='lightcoral',label='Discrete Uniform')
        p2 = a0.plot(xfine, betaBinomN(xfine,n,b1[0],b1[1]),linewidth=2,color='cornflowerblue',label='Beta Binomial')
        p3 = a0.plot(xfine, zipfDist(xfine,z1[0]),linewidth=2,color='chartreuse',label='Zipfian')
        
        a0.legend(loc='upper right',fontsize=12)
        
        rows = ['Chi-Square','R-Square','RMSE','K-S']
        cols = ['Discrete Uniform', 'Beta Binomial\n(a = %1.3f, b = %1.3f)' %(b1[0],b1[1]), 'Zipfian\n(c = %1.3f)' %(z1[0])]
        the_table = a1.table(cell_text,rowLabels=rows,colLabels=cols,loc='center',colColours=['lightcoral','cornflowerblue','chartreuse'])
        the_table.set_fontsize(10)
        a1.get_xaxis().set_visible(False)
        a1.get_yaxis().set_visible(False)
        plt.box(on=None)
        
        for j in range(0,3):
            cell = the_table[0,j]
            cell.set_height(0.3)
            
        
    #--------------------------------------------------------------------------
    # Output Data
        
    df_stats = pd.DataFrame({"Mean":[Rmean], 
                        "Std":[Rstd], 
                        "Range":[Rrange],
                        "Variance":[Rvar], 
                        "IQR":[Riqr]})
    
    df_dists = pd.DataFrame({"Fit Metric":["Chi-Square","R-Square","RMSE","K-S"],
                       "Discrete Uniform":[chis[0],Rs[0],RMSEs[0],KSs[0]],
                       "Beta Binomial":[chis[1],Rs[1],RMSEs[1],KSs[1]],
                       "Zipfian":[chis[2],Rs[2],RMSEs[2],KSs[2]]})
    
    
    ## comment as required 
    # out = [df_stats, df_dists, fittingParams] # output as arrays
    out = [basicStats, cell_data, fittingParams] # output as pandas data frames

    return out
