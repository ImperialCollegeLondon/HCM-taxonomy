# Association of common and rare HCM-associated variants on cardiac shape and motion

* [Variant Classification](https://github.com/ImperialCollegeLondon/HCM_expressivity/tree/master/variant_classification)  
Variant interpretation in HCM.  

* [2D Cardiac Image analysis](https://github.com/baiwenjia/ukbb_cardiac)   
Automated pipeline for image segmentation and motion analysis.  

* [3D Cardiac Image analysis](https://github.com/ImperialCollegeLondon/4DSegment2.0)  
Deep learning cardiac image analysis. 

* [3D Mass univariate regression analysis](https://github.com/ImperialCollegeLondon/HCM_expressivity/tree/master/3D_regression_analysis)  
3D modelling of genotype-phenotype relationship.  

* [Statistical analysis](https://github.com/ImperialCollegeLondon/HCM_expressivity/tree/master/statistical_analysis)  
Code for reproducing the tables and plots.  

## Background information:

### Mass Univariate regression
The basic principle of the mass univariate regression is to implement a linear regression at each vertex of the 3D atlas of the heart, extracting the regression coefficient associated with the variable of interest, in an attempt to establish – in the context of this pipeline - any potential associations between the genotype and the given phenotype under study (e.g. wall thickness of the left ventricle). The general linear model itself follows the form of:

 Y = βX + covariate_1 + covariate_2 + … + ε
 
Where Y is the dependent variable (e.g. phenotype) and X is the explanatory variable (e.g. genotype). This enables the approximation of the extent to which a given variable may be explained by another. In this notebook we are using the wall thickness of the left ventricle as the phenotype (Y), and whether or not a given individual possesses an HCM-associated variant (X) as the genotype. The example in this notebook is of the form:

Y = βX + age + sex + race + SBP + BSA + ε

Repeating this process at every vertex of the 3D atlas yields a beta map that contains all of the beta coefficients which represent any associations between the wall thickness of the left ventricle and the genotype of the individuals. 



### Threshold-free cluster enhancement (TFCE)
A few steps are taken to improve confidence in the robustness in the analysis. TFCE is a technique that calculates a score at each vertex of the atlas, the value of which depends on the extent of concordance between signals of nearby vertices. Regions of concordance (beta coefficients of the same sign) have their scores boosted. By contrast, regions of discordance (beta coefficients of opposite signs) have their scores reduced. 

The TFCE-derived p-values are calculated by performing permutations on the input data, applying the Freedman-Lane procedure. 

The derived p-values are subject to the multiplicity problem. Applying a technique such as the Bonferroni correction is often too harsh in many statistical analyses, but even moreso here. This is because of the underlying biology of the heart, where if there is hypertrophy in the LV it will be regional, as opposed to concentrated at one point. For example, if an association between the WT of the LV and the genotype at a given vertex is representative of a real biological effect, then it is more likely that surrounding vertices will cohere with this signal. Using the TFCE approach it is therefore possible to compute a statistic that accounts for regions of covariance at each vertex. The regular linear model yields a map of the t-statistics (t map), and from here, the TFCE approach can yield a TFCE map. Permutation testing enables the derivation of TFCE-derived p-values, and subsequent correction using the BH technique corrects for the multiplicity problem.  




### Data locations





