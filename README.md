# generative-model-parameter-fit
This repo includes key codes for "Parameter estimation for generative connectome models: Accuracy, reliability, and a fast parameter fitting method".
## Files included in this repo are:
1. generative_model_with_sequence.m: 
		matlab function modifying the work of Betzel et al (2016) on NeuroImage. Implementing the matching index model. 
		Modified code orders connections in synthetic network by their sequence of formation.
2. dependent_network_strategy_demo.m:
		a demo for the use of proposed dependent network strategy when fitting model parameters to individual connectomes
		(by calling the function generative_model_with_sequence.m).
3. sample_size_estimate.m:
		can be used to estimate the sample sizes required to detect expected between-group differences in generative model parameters.
4. dependent_network_strategy_data.mat:
		data used in the dependent_network_strategy_demo.m
5. sample_size_estimate_data.mat:
		data used in the sa,ple_size_estimate.m
* Note that 1 and 2 have dependency on Brain Connectivity Toolbox.

