# distributionFittingToolbox
Fit distributions to univariate integer data and use goodness-of-fit metric to compare the performance of each distribution

# Description
Takes plain-text input file containing a single integer on each line and fits the dicrete uniform, beta binomial and zipfian distributions to data distribution. 

A summary of the goodness-of-fit metrics (chi-square, R-square, RMSE, Kolmogorov-Smirnov) are displayed for each fit, along with general statistics of the input integer data and plot of the distribution and fits, along with teh fitting parameters for each distribution.

# How to use function (example)
[stats, dists, fit] = distributionFitting(filename, printResults, plotResults)

# Input arguments
filename: string of filename (and path) containing integer data

printResults: Boolean value (True/False), True prints results in console

plotResults: Boolean vale (True/False), True plots results

# Output
output[0] stats: 1x5 pandas data frame containing mean, standard deviation, range, variance and inter-quartile range of integer data

output[1] dists: 4x3 pandas data frame containing chi-square, R-square, RMSE and Kolmogorov-Smirnov metrics, for fittings integer distributions to the Discrete Uniform Distribution (DUD), Beta Binomial Distribution (BBD) and Zipfian Distribution (ZD)

output[2] fit:  1x3 fitting parameters a,b,c from BBD and ZD, defined in METHODS
 
# Structure/method of code
The structure of the code was chosen to split the different tasks into sections. This makes it easier for the user to understand what the scipt is doing and allows the user to easily add new distributions or goodness-of-fit parameters as they desire.

* Find frequency distribution of integer data [fr]

* Convert to probability distribution [p]

* Fit DUD to probability distribution
    - see https://mathworld.wolfram.com/DiscreteUniformDistribution.html
    - Each integer value has uniform probability of 1/n
    - n is number of integers in range of input data
    
* Fit BBD to probability distribution
    - see https://mathworld.wolfram.com/BetaBinomialDistribution.html
    - PDF: P(x) = nCx(n,x)*B(x+a,n-x-b)/B(a,b)
    - B: beta function
    - nCx(n,x): binomial coefficient of n choose x
    - n: len(x)
    - a,b: fitting parameters
    
* Fit ZD to probability distribution
    - see https://mathworld.wolfram.com/ZipfDistribution.html
    - PDF: P(x) = x^(-c) / z(c)
    - z: Riemann zeta function
    - c: fitting parameter

# Goodness-of-fit metrics
Given observed data, O, and expected data (fit), E

* chi-square
    - Used to test if a sample of data came from a specific ditribution
    - X^2 = sum_{i=1}^N (O_i - E_i)^2 / E_i
    - small chi = more likely to come from specific distribution
    
* R-square
    - Statistical measure of fit, indicate how much of variation of a dependent variable is explined by the indeendent variable
    - R^2 = 1 - (sum_{i=1}^N (O_i - E_i)^2)/(sum_{i=1}^N (O_i - mean(O))^2)
    - R^2 = 1 --> all variation of y-axis is account for by x-axis
    
*  Root mean square error (RMSE)
    - Measure of distance between observed and predicted model
    - RMSE = sqrt(1/N * sum_{i = 1}^N (O_i - E_i)^2)
    - small RMSE = smaller distance between observed and fitted data
    
* Kolmogorov-Smirnov (K-S)
    - see: 
    - Used with a sample from a population of an unknown distribution
    - FO: cumulative distribution of O
    - FE: cumulative distribution of E
    - Test statistic, D, gives greatest vertical distance between FO and FE
    - D = max_{i=1}^N | FO_i - FE_i|

# Libraries & versions
* numpy       (version 1.21.5)    
* scipy       (version 1.7.3)
* matplotlib  (version 3.5.1)
* pandas      (version 1.4.2)
* tabulate    (version 0.8.9)
* [python      (version 3.9.12)]

# Assumptions
1. The curve_fit function has optimised the fit of the distribution to the data and that increasing the number of iterations in this optimisation will non-negligibly affect the fit
2. The curve_fit optimisation function has found the global minimum and not a local minimum
3. The fit parameters a and b are both greater than zero
4. The user will include the path to the integer plain-text file in the input argument "filename", or include this path as a pythonpath
5. The plain-text file is delimited by "\n" (i.e. one integer per line ONLY)
6. The plain-text file is small enough to be loaded to the computer's RAM

# Limitations
* It can be slow and memory intensive to load all the data from the plain-text file if it is large
* If the curve fitting function has not converged after the 600 iterations, it will output NaN as the fitting coefficients (limitation for Beta Binomial and Zipfian distributions)
* Can be slow to append element to a python list, as the list is changing size after each element is added
