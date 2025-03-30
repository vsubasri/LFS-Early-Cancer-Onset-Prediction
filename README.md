# DNA methylation predicts early onset of primary tumor in patients with Li-Fraumeni Syndrome #

Li-Fraumeni syndrome (LFS) is an autosomal dominant cancer predisposition syndrome. Approximately 80% of LFS patients harbour a germline TP53 mutation rendering them susceptible to a wide spectrum of early onset malignancies. A comprehensive surveillance regimen termed the ‘Toronto Protocol’, has recently been adopted for early tumor detection, demonstrating significant improvement in survival among these patients. However, the protocol’s “one-size-fits-all” approach fails to consider an individual patient's risk of cancer. We have built a machine learning model that predicts early onset of primary tumors in LFS patients by estimating the probability of cancer onset before the age of six, leveraging a patient's peripheral blood leukocyte methylation profile. 

## Preprocessing ##

This directory contains a collection of scripts used to preprocess the DNA methylation data. This includes identifying an optimal preprocessing and normalization, outlier and confounder removal.

~~~
   /preprocessing
~~~

## Model Development ##

This directory contains a collection of scripts used to perform feature selection, cross validation and model selection.

~~~
    /model
~~~

## Feature Analysis ##

~~~
    /feature_analysis
~~~

## Resampling Analysis ##

~~~
    /sampling
~~~

## Immune Cell Proportions ##

~~~
    /idol
~~~


