import numpy as np

np.random.seed(10)

def generate_mutations_normal(num_ssm, clonal_freqs):
    """
    Generate mutations by sampling from normal distribution
    DEPRECATED
    """
    sim_mut = []
    for clonal_freq in clonal_freqs:
        if clonal_freq < 0 or clonal_freq > 1:
            raise ValueError("clonal frequency must be between 0 and 1")
        
        for _ in range(num_ssm):
            ssm_freq = np.random.normal(clonal_freq, 0.1)
            sim_mut.append(ssm_freq)
     
    return sim_mut

def generate_mutations_binom(num_ssm, clonal_freqs, read_depth):
    """
    Generate num_ssm simulated reference read counts for each clonal population.
    
    Args:
        num_ssms (int): number of ssms per population
        clonal_freqs (array of float): the population frequency of each subpopulation
        read_depth (int): simulation read depth
    Returns: 
        Array of sampled reference reads, and sampled 
    """
    sim_ref_counts = []
    sim_read_depth = []
    for clonal_freq in clonal_freqs:
        if clonal_freq < 0 or clonal_freq > 1:
            raise ValueError("clonal frequency must be between 0 and 1")

        for _ in range(num_ssm):
            num_reads = np.random.poisson(read_depth)
            ref_reads = np.random.binomial(num_reads, 1 - clonal_freq + 0.5 * clonal_freq) 
            if (num_reads == ref_reads):
                ref_reads = num_reads * (1 - clonal_freq + 0.5 * clonal_freq)
            sim_ref_counts.append(ref_reads)
            sim_read_depth.append(num_reads)
            
    return (sim_ref_counts, sim_read_depth)

def generate_phis(ref_counts, read_depths):
    """
    Generate phi values from reference read counts sampled from binomial distribution
    
    Args:
        ref_counts (int): number of sampled reference reads
        read_depths (int): total number of sampled reads
    Returns:
        Sorted array of phi values
    """
    phis = []

    for i in range(len(ref_counts)):
        r = ref_counts[i]
        v = read_depths[i] - r

        phi = 2 * np.random.beta(v + 1, r + 1)
        phis.append(phi)
        
    return phis

def calculate_phi(ref_counts, read_depth):
    """
    DEPRECATED: Phis should be re-sampled from variant/reference reads, as opposed to taken from them directly.

    Generate pseudo-phi values for simulation purposes.
    
    Args:
        ref_counts (array of int): number of reference reads for each mutation locus
        read_depth (array of int): read depth at each mutation locus
    Return:
        Array of phis
    """
    phis = []
    for i in range(len(ref_counts)):
        phi = 2 * (read_depth[i] - ref_counts[i]) / read_depth[i]
        phis.append(phi)
    
    return phis
    
def trim_zeros(generated_mutations):
    """
    Removes all zero read counts from generated mutations
    DEPRECATED
    """
    new_mutations = [[],[]]
    
    for i in range(len(generated_mutations[0])):
        if generated_mutations[0][i] != generated_mutations[1][i]:
            new_mutations[0].append(generated_mutations[0][i])
            new_mutations[1].append(generated_mutations[1][i])
    return new_mutations

def visualize(phis, clonal_freqs, num_ssms, read_depth):
    """
    Used to visualize variance for troubleshooting
    DEPRECATED
    """
    import matplotlib.pyplot as plt
    cm = plt.cm.get_cmap('viridis')
    c = ['r', 'b', 'y', 'k', 'c']
    mus = [0.56, 0.25, 0.06]
    var = []
    real = []
    for i in range(len(clonal_freqs)):
        print("Cluster at mean: " + str(mus[i]))
        cluster = phis[num_ssms*i:num_ssms*(i+1)]
        mu = mus[i]
        #est_var = ((2*(2-mu)* mu) / read_depth)
        est_var = ((2-mu)* mu / read_depth) + ((2 - mu) * mu / (read_depth + 1))
        print("Estimated variance: " + str(est_var))
        print("Actual variance: " + str(np.var(cluster)))
        print("")
        plt.scatter(phis[num_ssms*i:num_ssms*(i+1)], phis[num_ssms*i:num_ssms*(i+1)], 5, color=c[i])
        var.append(est_var)
        real.append(np.var(cluster))
    plt.show()
    return var, real

if __name__ == "__main__":
    num_ssms = 100
    clonal_freqs = [0.12, 0.45, 0.78]
    read_depth = 50
    
    simulated_data = generate_mutations_binom(num_ssms, clonal_freqs, read_depth)
