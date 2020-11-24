
# Mass Univariate Analysis 

A 3D modelling of genotype-phenotype relationships.  

## [functions](https://github.com/ImperialCollegeLondon/HCM_expressivity/tree/master/statistical_analysis/functions)

Consists of : 

* [multiply.cpp](https://github.com/ImperialCollegeLondon/HCM_expressivity/tree/master/statistical_analysis/functions/multiply.cpp) function uses "eigenMapMatMult" for faster matrices multiplication.

* [murq.R](https://github.com/ImperialCollegeLondon/HCM_expressivity/tree/master/statistical_analysis/functions/murq.R) for the mass univariate regression using the QR decomposition of the matrix into an orthogonal matrix and a triangular matrix.

* [permFL_fast](https://github.com/ImperialCollegeLondon/HCM_expressivity/blob/master/statistical_analysis/functions/permFL_fast.R) for the faster mass univariate regression analysis computing TFCE derived p-values using Freedman-Lane procedure.

This function first computes the nuisance matrix (Z) and the residual-forming matrix (Rz) using "eigenMapMatMult" for faster matrix multiplication. 

The procedure is then similar to permFL.R function from the package [mutools3D](https://github.com/UK-Digital-Heart-Project/mutools3D).

* TFCE and mur functions are also used from the package [mutools3D](https://github.com/UK-Digital-Heart-Project/mutools3D).

## [data](https://github.com/ImperialCollegeLondon/HCM_expressivity/tree/master/statistical_analysis/data)

Consists of the data used in the [Mass_univaraite_analysis_UKBB](https://github.com/ImperialCollegeLondon/HCM_expressivity/blob/master/statistical_analysis/Mass_univariate_analysis_UKBB.R).


The code has been fully documented and descriptions are available from within.