# Positive Unlabeled data simulation
To comprehensively assess the proposed method, we designed a series of simulation studies. Four different scenarios for the noise distribution are simulated that correspond to: population balanced with two classes well separated, or clear balanced scenario; population balanced with two classes not well separated, or noisy balance scenario; population unbalanced with two classes well separated, or clear unbalanced scenario; population unbalanced with two classes not well separated, or the noisy unbalanced scenario. The alteration of population unbalancedness and separation are achieved by designing different propensity score function as follows:

<div align=center> <img src="https://github.com/xiaoyulu95/PLUS/blob/master/fig/distribution.png" width="700"  height="150"> </div>

The function looks similar, but their non-linear properties are distinguished. Here we show the distribution of 𝑃𝑟(𝑌 = 1) for each of the scenarios. In addition, direction of the population unbalancedness is also altered, which either makes the positive instance the larger or smaller class.

<div align=center> <img src="https://github.com/xiaoyulu95/PLUS/blob/master/fig/data_simulation.png" width="500"  height="500"> </div>



## Usage
```
simul_data=PU_data_simulation(p=100,N=200,confident_rate=0.5,scenario='noisy_balance',valid='01')
```

## Arguments

* `p`: Number of variables

* `N`: Number of samples

* `confident_rate`: Proportion of positive samples

* `scenario`: Have four options `noisy_balance` `clear_balanced` `clear_unbalanced` `noisy_unbalanced` represent four different scenarios

* `valid`: `valid='01` means positive samples are from the right distribution, `valid='10` means positive samples are from the left distribution

## Value
Result list contains three elements: `train_data` takes the N by p simulated data matrix; `Label.obs` is the PU label for each sample under choosen confident_rate; `Label.true` takes true label from the data distribution . 

## Reference
Zhou, J., Lu, X., Chang, W., Wan, C., Zhang, C. and Cao, S., 2020. PLUS: predicting pan-cancer metastasis potential based on positive and unlabeled learning

