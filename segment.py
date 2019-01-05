import math
import numpy as np
from scipy.stats import norm
from scipy.integrate import quadrature

from scipy.integrate import quad

import warnings
warnings.simplefilter('error', RuntimeWarning)

import generate_simulations
import timeit

np.random.seed(10)

def pelt(phis, read_depth, score_function):
    """
    Generate partitioning of data using PELT method

    Args:
        phis (array of floats): list of phi values
        read_depth (int): average read depth for all phis
        score_function (function): segment scoring function
    
    Return:
        List containing optimal segmentation of mutations

    TODO:
        Make calculate read depth at segment level
        Figure out proper penalty
    """
    n_bins = len(phis)

    # TODO: Figure out how to properly penalize segment.
    n_sigs = 6
    penalty = n_sigs * np.log(n_bins * read_depth)

    # Matrix which contains scores for all subpartitions
    sp_scores = np.zeros((n_bins, n_bins))
    sp_scores[:] = np.nan

    # Array containing maximum score for subpartition of length index
    max_sp_scores = np.zeros(n_bins)

    # Set of potential changepoints to disclude
    prune_set = []

    for sp_len in range(n_bins):
        for last_cp in range(sp_len + 1):
            r_phis = phis[last_cp:sp_len + 1]

            if last_cp in prune_set:
                sp_scores[sp_len][last_cp] = -math.inf
                continue

            # Segments of length 1 are not considered
            if (len(r_phis) <= 1):
                sp_scores[sp_len, last_cp] = -math.inf
                continue
            
            # Currently using left and right most phi in segment as integration bounds.
            # Changing to mean between bounds does not greatly affect result.
            r_seg_score = 2 * score_function(r_phis, r_phis[0], r_phis[-1], read_depth)

            # Do not penalize if not introducing changepoint
            if last_cp == 0:
                sp_scores[sp_len, last_cp] = r_seg_score
            else:
                l_score = max_sp_scores[last_cp-1]
                sp_scores[sp_len, last_cp] = l_score + r_seg_score - penalty
        
        max_sp_scores[sp_len] = np.nanmax(sp_scores[sp_len,:])

        # Check all changepoints for prune condition
        for cp in range(sp_len):
            if sp_scores[sp_len, cp] == -math.inf:
                continue
            if sp_scores[sp_len, cp] + penalty < max_sp_scores[sp_len]:
                prune_set.append(cp)

    np.core.arrayprint._line_width = 180
    np.set_printoptions(suppress=True,
   formatter={'float_kind':'{:f}'.format})

    return sp_scores

def recover_changepoints(score_matrix):
    """
    Recover the segmentation resulting in the greatest score under the BIC.
    
    Args:
        score_matrix (2D matrix of floats): Log likelihoods of subsegments
    Returns: 
        Changepoints corresponding to optimal segmentation
    """
    opt_cps = [np.argmax(score_matrix[-1])]
    
    while opt_cps[-1] != 0:
        optimal_cp = np.nanargmax((score_matrix[opt_cps[-1] - 1]))
        opt_cps.append(optimal_cp)
    return opt_cps

def score_integ(phis, left_bound, right_bound, read_depth):
    """
    Score a segment by taking the integral over all means in segment
    
    Args:
        phis (list of float): vaf of mutations in sub-segment
        read_depth (int): average read depth in sub-segment
    """

    integ = quad(norm_likelihood, left_bound, right_bound, args=(phis, read_depth, left_bound, right_bound))

    # TODO: Integration boundaries have issues. Current solution might be a hack.
    if (integ[0] == 0):
        return -math.inf
    
    if (integ[0] < 0):
        return 0

    score = np.log(integ[0])
    
    return score

def score_mle(phis, left_bound, right_bound, read_depth):
    """
    Score a segment using maximum likelihood
    """
    return sum(norm.logpdf(phis, np.mean(phis), np.sqrt(np.var(phis))))

def norm_likelihood(mean, phis, read_depth, left_bound, right_bound):
    """
    Compute likelihood under normal approximation

    Args:
        mean (float): Integration variable under quad
        phis (list of float): All ssm phi's in segment
        read_depth (int): Average read depth in segment
        left_bound (float): Integration left bound
        right_bound (float): Integration right bound

    Returns:
        Likelihood under normal approximation
    """
    n_phi = len(phis)
    variance = 2 * (((2 - mean) * mean )) / read_depth

    # Below: Actual variance, but does not greatly affect result
    #variance = 2 * ((2-mean)* mean / read_depth) + ((2 - mean) * mean / (read_depth + 1))

    const_numr = 0.5 * (np.log(2 * np.pi * variance) - np.log(n_phi))
    const_denom = 0.5 * n_phi * np.log(2 * np.pi * variance)
    const = const_numr - const_denom
    
    exp_term = (-n_phi / (2 * variance)) * (np.mean(np.square(phis)) - np.square(np.mean(phis)))
    
    pdf = norm.logpdf(np.mean(phis), mean, np.sqrt(variance / n_phi))

    # Make absolute to avoid sorting idiosyncrasies
    uniform_prior = -np.log(np.abs(right_bound - left_bound))
    
    score = uniform_prior + const + exp_term + pdf
    exp_score = np.exp(score)
    return exp_score

def calculate_mean_from_cp(phis, changepoints):
    """
    Determine subpartition mean from changepoint boundaries

    Args:
        phis (list of float): All ssms inputted into pelt
        changepoints (list of int): Changepoints detected by PELT

    Returns:
        List of means of each subclonal cluster
    """
    if changepoints[0] == 0:
        changepoints = changepoints[1:len(changepoints)]
    
    if len(changepoints) == 0:
        return [phis], [np.mean(phis)]

    
    assignments = []
    assignments.append(phis[0:changepoints[0]])
    for i in range(len(changepoints)-1):
        assignments.append(phis [changepoints[i] : changepoints[i+1]])
    assignments.append(phis[changepoints[len(changepoints)-1]:len(phis)])

    means = []

    for i in range(len(assignments)):
        means.append(np.mean(assignments[i]))

    return assignments, means

def run(ssm_per_pop, pop_freqs, read_depth):
    """
    Runs PELT on simulated data.

    Args:
        ssm_per_pop (int): Number of mutations per subclonal population
        pop_freqs (float): Population frequency of subclonal population
        read_depth (int): Average read depth over all samples
    """
    
    print(ssm_per_pop, pop_freqs, read_depth)

    # Generate simulated data
    counts, depths = generate_simulations.generate_mutations_binom(ssm_per_pop, pop_freqs, read_depth)
    phis = generate_simulations.generate_phis(counts, depths)

    phis = sorted(phis)
    start = timeit.default_timer()

    # Run PELT
    score_matrix = pelt(phis, read_depth, score_integ)
    recovered_cps = recover_changepoints(score_matrix)
    
    recovered_cps.reverse()

    stop = timeit.default_timer()

    print("Detected changepoints: ", recovered_cps)

    # Visualize
    split_points = []
    if len(recovered_cps) == 1:
        split_points.append(phis[:recovered_cps[0]])
        split_points.append(phis[recovered_cps[0]:])
    else:
        for i in range(len(recovered_cps)):
            if i is 0:
                split_points.append(phis[0:recovered_cps[i]])
            elif i is len(recovered_cps)-1:
                split_points.append(phis[recovered_cps[i-1]:recovered_cps[i]])
                split_points.append(phis[recovered_cps[i]:len(phis)])
            else:
                split_points.append(phis[recovered_cps[i-1]:recovered_cps[i]])

    import matplotlib.pyplot as plt
    cm = plt.cm.get_cmap('Set3')

    for points in range(len(split_points)):
        plt.scatter(split_points[points], split_points[points], 5,color=cm(points/len(split_points)))
    plt.show()
  
    assignments, means = calculate_mean_from_cp(phis, recovered_cps)

    print("Detected cluster means:" + str(means))
    
    flat_list = [item for sublist in assignments for item in sublist]

    assert len(flat_list) == len(phis)
    
    print("Number of segments: " + str(len(means)))
    print("=========================================")
    print("")

def batch_run(rd_array, ssm_array):
    """
    Aggregate run on simulated data generated accoring to param arrays.

    Args:
        rd_array (list of int): All desired read depths
        ssm_array (list of int): All desired ssm_per_populations
    """

    pop_freqs = [[0.71, 0.44, 0.25, 0.11, 0.03]]
    #pop_freqs = [[0.44, 0.11], [0.56, 0.25, 0.06], [0.64, 0.36, 0.16, 0.04], [0.71, 0.44, 0.25, 0.11, 0.03]]

    for pop_freq in pop_freqs:
        for rd in rd_array:
            for ssm in ssm_array:
                run(ssm, pop_freq, rd)
    
    '''
    DEPRECATED: Used to troubleshoot appropriate variances
    varis = [[], [], [] ]
    real = [[], [], []]
    for pop_freq in pop_freqs:
        for rd in rd_array:
            for ssm in ssm_array:
                v, r = run(ssm, pop_freq, rd)
                varis[0].append(v[0])
                real[0].append(r[0])
                varis[1].append(v[1])
                varis[2].append(v[2])
                real[1].append(r[1])
                real[2].append(r[2])
    
    import matplotlib.pyplot as plt

    plt.plot(rd_array, varis[0], color='red')
    plt.plot(rd_array, real[0], color='blue')
    plt.show()
    
    plt.plot(rd_array, varis[1], color='red')
    plt.plot(rd_array, real[1], color='blue')
    
    plt.show()
    
    plt.plot(rd_array, varis[2], color='red')
    plt.plot(rd_array, real[2], color='blue')
    
    plt.show()
    '''

if __name__ == "__main__":

    ssm_per_pop = 100
    pop_freqs = [.64,.36,.16]
    read_depth = 50

    run(100, [0.64, 0.36, 0.16], 50)
    #batch_run([300], [50])
