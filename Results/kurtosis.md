Analysing the shape and variance of the distributions
===

After getting rid of the centromeres, I calculated the standard deviation, kurtosis and skewness of the whole chromosome distributions. 

--------

- B & C are [mob1/mob2](http://www.sciencedirect.com/science/article/pii/S1931312814003850)

- A, E, F - unpublished

- [OCF2](http://onlinelibrary.wiley.com/doi/10.1111/j.1365-313X.2012.04993.x/full#ss9)

- [BCF2](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3772335/#SM3)

- [sup#1 chr1 & chr4](http://pcp.oxfordjournals.org/content/52/4/716.long)

------

<table>
  <tr><th>Sample <th>Length</th><th>Chromosome</th><th>Correlation with normal distribution (r<sup>2</sup>)</th><th>Standard deviation (Mb)</th><th>Kurtosis</th><th>Skewness</th>
  
    
  <tr><th>A </th> <th>16</th> <th> 2</th><th>~0.9<th>6.99 </th><th>1.62</th><th>-0.559</th>
  
  <tr><th>B </th>  <th>25</th><th> 5</th><th>>0.9<th>7.20 </th><th>2.69</th><th>-0.240</th>
 
  
   <tr><th>C </th> <th>41</th> <th> 5</th><th>~0.9<th>3.71</th><th>2.37</th><th>0.445</th>  
   
   <tr><th> E</th>  <th>18</th><th> 2</th><th>>0.9<th>2.34 </th><th>2.17</th><th>-0.104</th>
  
  <tr><th> F</th>  <th>93</th><th>1</th><th>>0.95<th>4.93</th><th>3.15</th><th>0.021</th>
  
  
  
  <tr><th> OCF2</th> <th>151</th> <th>2</th><th>>0.9<th>6.01</th><th>1.85</th><th>-0.392</th>
  
  <tr><th> BCF2</th>  <th>15</th><th>3</th><th>>0.9<th>3.20</th><th>2.02</th><th>0.201</th>
  
  
  <tr><th> sup#1</th> <th>2165</th> <th>1</th><th>~0.87<th>7.36</th><th>3.43</th><th>-1.138</th>
  
  
  
  <tr><th> sup#1</th>  <th>4633</th><th>4</th><th>>0.95<th>3.66</th><th>3.50</th><th>0.370</th>
  
</table>


####Standard deviation
If we ignore the example that doesn't correlate so well with the normal distribution (lower r<sup>2</sup>), the standard deviation will oscilate between 2.34 and 7.2 Mb. It is necessary to take also into account the size of the sample (length) 

####Kurtosis

Kurtosis reflects the shape of a distribution apart from the variance ([DeCarlo, 1997](http://www.columbia.edu/~ld208/psymeth97.pdf)). The normal distribution has a kurtosis of 3, so if kurtosis > 3 (positive), distributions are _**leptokurtic**_. In terms of shape, a leptokurtic distribution has a **more acute peak around the mean and fatter tails**, while a  _**platykurtic**_ (negative) distribution shows a lower, **wider peak around the mean and thinner tails** ([wikipedia](http://en.wikipedia.org/wiki/Kurtosis)). 

3 out of 9 samples shown here are leptokurtic. The kurtosis oscilates between +0.15 and +0.50. The remaining six examples show a negative kurtosis, with a value oscilating between -0.31 and -1.38. 

As for the _**skewness**_, this value implies that the distribution of the data is skewed to the left (or negatively skewed) or to the right (positively skewed). The skewness will depend on the location of the causal mutation (peak) in the chromosome, so I don't expect to see any conserved tendency. 
