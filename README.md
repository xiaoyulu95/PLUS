# PLUS
**P**ositive and unlabeled **L**earning from **U**nbalanced cases and **S**parse structures, or **PLUS**, represents the first one to use positive and unlabeled learning framework to specifically model the under-diagnosis issue in predicting cancer metastasis potential. PLUS is specifically tailored for studying metastasis that deals with the unbalanced instance allocation as well as unknown metastasis prevalence, which are not capable by any other methods. Its robustness grants the possibility to harness the power of big data by integrating large scale datasets from different cancer types. Insights gleaned from this research will prove useful to the diagnosis and treatment of clinical metastatic disease.

## The Motivation of PLUS
![image](https://github.com/xiaoyulu95/PLUS/blob/master/fig/Figure1.png)

(a) Among the patients who were diagnosed as non-metastatic, some were under-diagnosed. Traditional classification using diagnosis as response may underestimate the cancer metastasis potential. PLUS is designed to recognize the bias in under-diagnosis, so that patients with higher metastasis potential could be accurately classified. (b) In TCGA Pan-cancer study, for patients who are clinically diagnosed as non-metastatic (M0) at baseline in each cancer type (columns), the top three rows shows the proportions of patients with follow-up information who were found alive and with non- progressed disease (NP-Alive), alive and with progressed disease (P-Alive), and dead (Dead); and the bottom three rows show the same proportions for patients who were diagnosed as metastatic (M1) at baseline. (c) The median follow-up time for patients who were diagnosed as non-metastatic (blue) and as metastatic (yellow) at baseline for each cancer type.


## Installation

Download [**PLUS function**](https://github.com/xiaoyulu95/PLUS/blob/master/R/PLUS.R)

```
#install dependent pkg

install.packages("glmnet")

```

## Usage
```
Prediction=PLUS(train_data=X,Label.obs=Label,Sample_use_time=30,l.rate=1,qq=0.1)
```

## Arguments

* `train_data`: Gene expression matrix (N by M) which has N samples and M variables

* `Label.obs`: Positive Unlabeled observation label for each sample, a binary vector, 1 means true positive labels, 0 means unlabeled labels

* `Sample_use_time`: Used in stop criteria, how many times each samples to be used in training process

* `l.rate`: Control how much information from last iteration will be used in next

* `qq`: Quantile of the probability for positive samples, used to determine the cutoff between positive and negtive

## Value
Result list contains three elements: `pred.y` shows the probability for each same to be predicted as positive; `cutoff` is the reference cutoff to transfer continues probability to binary 0/1 label; `pred.coef1` take the variable coefficient used in prediction model. 

## Example

We provide an [exmaple data](https://github.com/xiaoyulu95/PLUS/blob/master/data/example_data.RData)
```
### The R packages involved in PLUS package
library(glmnet)

source('PLUS.R')
```

```
X=example_data$train_data
```

X is input data matrix each row is a sample and each column is a variable. 

<div align=center> <img src="https://github.com/xiaoyulu95/PLUS/blob/master/fig/example_data.png" width="700"  height="100"> </div>



```
Label=example_data$Label.obs
```

Observation label for each sample. In cancer metastasis prediction, 1 means true metastasis, 0 means unlabeled samples.



```
Prediction=PLUS(train_data=X,Label.obs=Label,Sample_use_time=30,l.rate=1,qq=0.1)
```

`Prediction$pred.y` is a vector which contians the probability for each patient to be predicted as metastasis. You can either directly use the probability or use `Prediction$cutoff` as a cutoff to classify probabilities into binary Label. `pred_Label` is the final label for patient metastasis prediction. Selected genes in prediction can be achieved using `Prediction$pred.coef1`.

```
pred_Label=ifelse(Prediction$pred.y<Prediction$cutoff,0,1)
```

<!--![image](https://github.com/xiaoyulu95/PLUS/blob/master/fig/density.png)-->



<div align=center> <img src="https://github.com/xiaoyulu95/PLUS/blob/master/fig/density.png" width="500"  height="500"> </div>

<div align=center> Density plot for predicted probabilities and corresponding cutoff. </div>



## PU data simulation
PLUS package also provide the capability to simulate positive unlabeled data in different setting. Detail: https://github.com/xiaoyulu95/PLUS/tree/master/PU%20data%20simulation

## Validation on Single Cell RNA-seq 
Data and code for single cell RNA-seq analysis are available in: https://github.com/xiaoyulu95/PLUS/tree/master/SC%20data


## Contact Information

- [Xiaoyu Lu](https://zcslab.github.io/people/xiaoyu/)
(lu14@iu.edu)

Ph.D. candidate, Indiana University School of Medicine

- [Junyi Zhou](https://fsph.iupui.edu/about/directory/zhou-junyi.html)
(junyzhou@iu.edu)

Ph.D. candidate, Department of Biostatistics, Indiana University


## Reference
Zhou, J., Lu, X., Chang, W., Wan, C., Zhang, C. and Cao, S., 2020. PLUS: predicting pan-cancer metastasis potential based on positive and unlabeled learning

