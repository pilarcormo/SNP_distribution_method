Analysing the shape and variance of the distributions
===

After getting rid of the centromeres, I calculated the standard deviation, kurtosis and skewness of the whole chromosome distributions. 

--------

- B & C are [mob1/mob2](http://www.sciencedirect.com/science/article/pii/S1931312814003850)

- F - unpublished

- [OCF2](http://onlinelibrary.wiley.com/doi/10.1111/j.1365-313X.2012.04993.x/full#ss9)

- [BCF2](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3772335/#SM3)

- [sup#1 chr1 & chr4](http://pcp.oxfordjournals.org/content/52/4/716.long)

------

<table>
  <tr><th>Sample <th>Chromosome</th><th>Correlation with normal distribution (r<sup>2</sup>)</th><th>Standard deviation (Mb)</th><th>Kurtosis</th><th>Skewness</th>
  
  <tr><th>B </th> <th> 5</th><th>>0.9<th>7.20 </th><th>2.69</th><th>-0.240</th>
  
  
  
   <tr><th>C </th> <th> 5</th><th>>0.9<th>5.79 </th><th>4.48</th><th>-0.967</th>  
  
  <tr><th> F</th> <th>1</th><th>>0.9<th>4.93</th><th>3.15</th><th>0.021</th>
  
  
  
  <tr><th> OCF2</th> <th>2</th><th>>0.9<th>6.01</th><th>1.85</th><th>-0.392</th>
  
  <tr><th> BCF2</th> <th>3</th><th>>0.9<th>3.20</th><th>2.02</th><th>0.201</th>
  
  
  <tr><th> sup#1</th> <th>1</th><th>~0.87<th>7.36</th><th>3.43</th><th>-1.138</th>
  
  
  
  <tr><th> sup#1</th> <th>4</th><th>>0.9<th>3.66</th><th>3.50</th><th>0.370</th>
  
</table>


####Standard deviation
If we ignore the example that doesn't correlate so well with the normal distribution (lower r<sup>2</sup>), the standard deviation will oscilate between 3.2 and 7.2 Mb .

####Kurtosis

Kurtosis reflects the shape of a distribution apart from the variance ([DeCarlo, 1997](http://www.columbia.edu/~ld208/psymeth97.pdf)). The normal distribution has a kurtosis of 3, so if kurtosis > 3 (positive), distributions are _**leptokurtic**_. In terms of shape, a leptokurtic distribution has a **more acute peak around the mean and fatter tails**, while a  _**platykurtic**_ (negative) distribution shows a lower, **wider peak around the mean and thinner tails** ([wikipedia](http://en.wikipedia.org/wiki/Kurtosis)). 

4 out of 6 samples shown here are leptokurtic. The kurtosis oscilates between +0.15 and +1.48. Three examples (mob1, OCF2, BCF2) show a negative kurtosis, with a value oscilating between -0.31 and -1.15. 

As for the _**skewness**_, this value implies that the distribution of the data is skewed to the left (or negatively skewed) or to the right (positively skewed). The skewness will depend on the location of the causal mutation (peak) in the chromosome, so I don't expect to see any conserved tendency. The 3 examples with negative skewness present the mutation at the end of the chromosome, while the other 3 examples with positive skewness present the mutation at the beginning or middle of the chromosome. 
