# generative-model-parameter-fit
This repo includes codes to perform FLaG(Fast landscape generation)-implemented parameter estimation for connectome generative models as demonstrated in "Parameter estimation for connectome generative models: Accuracy, reliability, and a fast parameter fitting method".
## Files included in this repo are:
1. **flag_search.m:** 

    This is the function that will be called to perform parameter estimation with FLaG.

    Brief explaination of the FLaG:
        In the manuscript, we described the algorithm as assigning networks in growth to subjects with matched network density. However, in practice, we implement the algorithm in a reversed order. In this Matlab function, We first generated a network whose density matched to the subject with the most connections, and then gradually remove connections from the latest generated ends. For example, if a dataset has 3 subjects with 300, 305,310 connections respectively. We first generate a network with 310 connections for subject 3, next we remove the last 5 and 10 connections to generate FLaG networks for the first two subjects. 

2. **Peripheral function in peripheral folder**
        
    These are peripheral functions that are called by flag_search.m

    Including:
        generate_connections.m:     function modified from the original generative_model function by Betzel et al (2016). Modifications are made to return generated connections according to their generation orders.
        energy_from_flag.m:         function that compute the energy of a FLaG landscape networks against a group of population.
        network_topo.m:             function returns the four topological distributions of interest (i.e., k, c, b, e) for given networks.
        fcn_ks.m:                   function compute the ks statistics from two topological distributions (copy from Betzel et al (2016)).
        fcn_voronoi_select.m:       function that select parameters for Voronoi Tessellation search (copy from Betzel et al (2016)). Not used in the grid FLaG. Included for future implementation of FLaG to Voronoi tessellation.

3. **multilandscape_estimate.m:**

    Function that performs multilandscape estimation from landscape parameter points and their energy costs, as demonstrated in the manuscript.

4. **demo.m:**

    This is a demo code illustrate the use of flag_search.m and multilandscape_estimate.m functions to perform FLaG-implemented parameter estimation.
		
5. **demo_data.mat:**

    Data used in the demo.m for illustration purpose
		
* **Note that flag_serach.m have dependency on Brain Connectivity Toolbox.**

