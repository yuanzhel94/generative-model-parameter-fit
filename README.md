# generative-model-parameter-fit
This repo includes key codes for "Parameter estimation for connectome generative models: Accuracy, reliability, and a fast parameter fitting method".
## Files included in this repo are:
1. **flag.m:** 
		
    Matlab function modifying the work of Betzel et al (2016) on NeuroImage. Implementing the matching index model. Modified code orders connections in synthetic network by their sequence of formation.

2. **flag_demo.m:**

    A demo for the use of FLaG pipeline when fitting model parameters to individual connectomes (by calling the function flag.m).

3. **sample_size_estimate.m:**

    Can be used to estimate the sample sizes required to detect expected between-group differences in generative model parameters.

4. **flag_data.mat:**

    Data used in the flag_demo.m
		
5. **sample_size_estimate_data.mat:**

    Data used in the sample_size_estimate.m
		
* **Note that 1 and 2 have dependency on Brain Connectivity Toolbox.**

