# import the function
import distributionFitting as dFit

# run function
# INPUTS
# * filename: string to file location
# * printResults: Booldean value, True to print results in console
# * plotResults: Boolean value, True to plot results
#
# OUTPUTS
# * stats: array containing [mean, standard deviation, range, variance, inter-quartie range] of input integr data
# * dists
[stats,dists,fit] = distributionFitting(filename = 'exampleInput.txt', printResults = True, plotResults = True)                                     
