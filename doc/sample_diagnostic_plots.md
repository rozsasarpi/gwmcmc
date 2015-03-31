
Diagnostic plots
=================

### Table of Contents
* [Autocorrelation plot](#autocorrelation-plot)
* [Trace plot](#trace-plot)
* [Triangle plot](#triangle-plot)
* [Histogram plot](#histogram-plot)

These illustrative diagnostic plots are corresponding to a dummy example:
* trivariate normal distribution
* inferring the the mean from sample
* vague priors

Autocorrelation plot
--------------------
![autocorr_01]

Trace plot
----------
![trace_01]

Triangle plot
-------------
The _hd_ refers to highest density credible interval (_CI_). The diagonal elements are the posterior distribution functions of the parameters.
The red lines at each diagonal plot represent the credible intervals.
The lower triangle shows the pairwise 2D marginal, posterior distributions with 2D credible contours. These are also highest density contours, and 2D kernel density functions are used to construct them.
The upper triangle indicates the Spearman correlation between the corresponding two samples. The fontsize is proportional to the correlation.
![triangle_01]

Histogram plot
--------------
The _hd_ refers to highest density credible interval (_CI_). The red area contains the 90% of the probability distribution function.

![1Dhist_01]

[autocorr_01]: https://github.com/rozsasarpi/gwmcmc/blob/master/doc/autocorrelation_plot_sample_01.png "autocorrelation plot"
[trace_01]: https://github.com/rozsasarpi/gwmcmc/blob/master/doc/trace_plot_sample_01.png "trace plot"
[triangle_01]: https://github.com/rozsasarpi/gwmcmc/blob/master/doc/triangle_plot_sample_01.png "triangle plot"
[1Dhist_01]: https://github.com/rozsasarpi/gwmcmc/blob/master/doc/1D_marginal_plot_sample_01.png "1D posterior marginal plot"
