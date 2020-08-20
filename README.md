## Project

__Inter-trial effects in priming of pop-out: Comparison of Bayesian updating models__

## Team members

* Fredrik Allenmark<sup>1</sup>
* Ahu Gokce<sup>2</sup>
* Thomas Geyer<sup>1</sup>
* Hermann J. Müller<sup>1</sup>
* Zhuanghua Shi<sup>1</sup>

1. Department of Psychology, Ludwig-Maximilians-Universität München, Munich, Germany
2. Department of Psychology, Kadir Has University, Istanbul, Turkey

## Abstract

In visual search tasks, repeating features or the position of the target results in faster response times. Such inter-trial ‘priming’ effects occur not just for repetitions from the immediately preceding trial but also from trials further back. A paradigm known to produce particularly long-lasting inter-trial effects – of the target-defining feature, target position, and response (feature) – is the ‘priming of pop-out’ (PoP) paradigm, which typically uses sparse search displays and random swapping across trials of target- and distractor-defining features. However, the mechanisms underlying these inter-trial effects are still not well understood. To address this, we applied a modeling framework combining an evidence accumulation (EA) model with different Bayesian updating rules of the model parameters (i.e., the drift rate and starting point of EA) for different aspects of stimulus history, to data from a (previously published) PoP study that had revealed significant inter-trial effects from several trials back for repetitions of the target color, the target position, and (response-critical) target feature. By performing a systematic model comparison, we aimed to determine which EA model parameter and which updating rule for that parameter best accounts for each inter-trial effect and the associated n-back temporal profile. We found that, in general, our modeling framework could accurately predict the n-back temporal profiles. Further, target color- and position-based inter-trial effects were best understood as arising from redistribution of a limited-capacity weight resource which determines the EA rate. In contrast, response-based inter-trial effects were best explained by a bias of the starting point towards the response associated with a previous target; this bias appeared largely tied to the position of the target. These findings elucidate how our cognitive systems continually tracks, and updates an internal predictive model of, a number of separable stimulus and response parameters in order to optimize task performance.

## Acknowledgements

This work was supported by German research foundation DFG MU773/16-1 to HJM and DFG SH166/3-1 to ZS. 

## Contents

The file data_analysis.R contains the code for plotting the main behavioral results figure, showing the inter-trial effects of target color, target position and the response defining target notch position.

The file factorial_Models.R runs the model comparison. 

The file models_analysis.R creates figures of relative AIC, temporal profiles of the inter-trial effects and model predictions as well as examples of the updating process for the updating rules of the best model. This uses the results saved by running factorial_Models.R.

The files bayesian_updates.R and figure_functions.R contains various functions that are required by the model comparison, e.g. implementations of the different updating rules, and/or in order to create the figures.