### Date:  2020-05-04
### Title: MetaboAnalyst on metabolite data

#### Exploring PLS-DA with auto-scaled data

I ran MetaboAnalyst on the web interface using the metabolite data. 
I used autoscaling. There were no zero values and metabolites were detected in all samples so I did not do any additional filtering at this step.

The PLS-DA gave the plots below (which should be the same that Krista sees):

[![](https://raw.githubusercontent.com/shellytrigg/pteropod_pHxDO_metabolomics/master/Multivariate_analysis/pls_pair_autoscale_all_data.png)](https://raw.githubusercontent.com/shellytrigg/pteropod_pHxDO_metabolomics/master/Multivariate_analysis/pls_pair_autoscale_all_data.png)

I identified outlier samples from the PCA selecting the show sample names option and found samples M11A4_013, M12A2_018, and M4D1_001 to be outliers. 

I excluded these samples from the analysis, redid the autoscaling normalization and reran PLS-DA. I got the following plots:

[![](https://raw.githubusercontent.com/shellytrigg/pteropod_pHxDO_metabolomics/master/Multivariate_analysis/pls_pair_autoscale_exclude_outliers.png)](https://raw.githubusercontent.com/shellytrigg/pteropod_pHxDO_metabolomics/master/Multivariate_analysis/pls_pair_autoscale_exclude_outliers.png)

There is not much difference when outliers are removed

#### PLS-DA class order

By default, for PLS-DA 'class order matters' is checked off. The description next to the check box says "the default approach is meaningful when the group labels correspond to time series, disease severity, or treatment dosages. However, when group labels do not reflect quantitative differences, the user should uncheck the option "class order matters." In the crab paper, I had just used the default settings, but it is seeming like that might not be the best thing to do here. A point for discussion in our next meeting. 

When this box is unchecked, the PLS-DA gives the following plots: 

[![](https://raw.githubusercontent.com/shellytrigg/pteropod_pHxDO_metabolomics/master/Multivariate_analysis/pls_pair_autoscale_all_data_noClassOrder.png)](https://raw.githubusercontent.com/shellytrigg/pteropod_pHxDO_metabolomics/master/Multivariate_analysis/pls_pair_autoscale_all_data_noClassOrder.png)

When this box is unchecked AND outliers are removed, the PLS-DA gives the following plots: 

[![](https://raw.githubusercontent.com/shellytrigg/pteropod_pHxDO_metabolomics/master/Multivariate_analysis/pls_pair_autoscale_exclude_outliers_noClassOrder.png)](https://raw.githubusercontent.com/shellytrigg/pteropod_pHxDO_metabolomics/master/Multivariate_analysis/pls_pair_autoscale_exclude_outliers_noClassOrder.png)

#### PLS-DA Permutation test

I ran permutation tests of the PLS-DA results selecting 'prediction accuracy during training' as the test statistic.

I generated the following plots:

- class order matters, no outliers removed:  [![](https://raw.githubusercontent.com/shellytrigg/pteropod_pHxDO_metabolomics/master/Multivariate_analysis/Screen%20Shot%202020-05-04%20at%202.23.33%20PM.png)](https://raw.githubusercontent.com/shellytrigg/pteropod_pHxDO_metabolomics/master/Multivariate_analysis/Screen%20Shot%202020-05-04%20at%202.23.33%20PM.png)
- class order matters, outliers removed:  [![](https://raw.githubusercontent.com/shellytrigg/pteropod_pHxDO_metabolomics/master/Multivariate_analysis/Screen%20Shot%202020-05-04%20at%202.29.09%20PM.png)](https://raw.githubusercontent.com/shellytrigg/pteropod_pHxDO_metabolomics/master/Multivariate_analysis/Screen%20Shot%202020-05-04%20at%202.29.09%20PM.png)
- class order DOESN'T matter, outliers removed:  [![](https://raw.githubusercontent.com/shellytrigg/pteropod_pHxDO_metabolomics/master/Multivariate_analysis/Screen%20Shot%202020-05-04%20at%202.30.29%20PM.png)](https://raw.githubusercontent.com/shellytrigg/pteropod_pHxDO_metabolomics/master/Multivariate_analysis/Screen%20Shot%202020-05-04%20at%202.30.29%20PM.png)

I think the reason Krista's permutation tests look different are because they were run with the 'separation distance' test statistic. When I selected this, I produced the following plot that looks like Krista's:

- class order matters, no outliers removed:  [![](https://raw.githubusercontent.com/shellytrigg/pteropod_pHxDO_metabolomics/master/Multivariate_analysis/Screen%20Shot%202020-05-04%20at%203.17.14%20PM.png)](https://raw.githubusercontent.com/shellytrigg/pteropod_pHxDO_metabolomics/master/Multivariate_analysis/Screen%20Shot%202020-05-04%20at%203.17.14%20PM.png)

In the code the difference in these two methods is denoted by "bw" or "accu":

- 'prediction accuracy during training' selected:
	- mSet<-PLSDA.Permut(mSet, 1000, "accu") 
- 'separation distance' selected: 
	- mSet<-PLSDA.Permut(mSet, 1000, "bw")

**R history is here:** [https://github.com/shellytrigg/pteropod_pHxDO_metabolomics/blob/master/Multivariate_analysis/Rhistory.R](https://github.com/shellytrigg/pteropod_pHxDO_metabolomics/blob/master/Multivariate_analysis/Rhistory.R)
