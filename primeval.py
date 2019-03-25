import msprime
import math

###############################################################################
# Gutenkunst et al., 2009
###############################################################################
def gutenkunst_model(mu=1.5e-8, phi=0, length=1e4, n_afr=0, n_eas=0, n_eur=0, seed=100, debug=False):
    # First we set out the maximum likelihood values of the various parameters
    # given in Table 1.
    N_A = 7300
    N_B = 2100
    N_AF = 12300
    N_EU0 = 1000
    N_AS0 = 510
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_AF = 220e3 / generation_time
    T_B = 140e3 / generation_time
    T_EU_AS = 21.2e3 / generation_time
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU = 0.004
    r_AS = 0.0055
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    # Migration rates during the various epochs.
    m_AF_B = 25e-5
    m_AF_EU = 3e-5
    m_AF_AS = 1.9e-5
    m_EU_AS = 9.6e-5
    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=n_afr, initial_size=N_AF),
        msprime.PopulationConfiguration(
            sample_size=n_eur, initial_size=N_EU, growth_rate=r_EU),
        msprime.PopulationConfiguration(
            sample_size=n_eas, initial_size=N_AS, growth_rate=r_AS)
    ]
    migration_matrix = [
        [      0, m_AF_EU, m_AF_AS],
        [m_AF_EU,       0, m_EU_AS],
        [m_AF_AS, m_EU_AS,       0],
    ]
    demographic_events = [
        # CEU and CHB merge into B with rate changes at T_EU_AS
        msprime.MassMigration(
            time=T_EU_AS, source=2, destination=1, proportion=1.0),
        msprime.MigrationRateChange(time=T_EU_AS, rate=0),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
        msprime.MigrationRateChange(
            time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
        msprime.PopulationParametersChange(
            time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),
        # Population B merges into YRI at T_B
        msprime.MassMigration(
            time=T_B, source=1, destination=0, proportion=1.0),
        # Size changes to N_A at T_AF
        msprime.PopulationParametersChange(
            time=T_AF, initial_size=N_A, population_id=0)
    ]
    
    if debug:
        # Use the demography debugger to print out the demographic history
        # that we have just described.
        dd = msprime.DemographyDebugger(
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,
            demographic_events=demographic_events)
        dd.print_history()
    else:
        sim = msprime.simulate(population_configurations=population_configurations,
                               migration_matrix=migration_matrix, 
                               demographic_events=demographic_events,
                               mutation_rate=mu, 
                               recombination_rate=phi, 
                               length=length,
                               random_seed=seed)
        return sim

###############################################################################
# Fu et al., 2013
###############################################################################
def fu_model(mu=1.5e-8, phi=0, length=1e4, n_afr=0, n_eur=0, debug=False):
    # First we set out the maximum likelihood values of the various parameters
    # given in Table 1.
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    
    # 220kya:
    # African population constant with Ne~7300
    N_A = 7310
    
    # 148kya:
    # instantaneous growth to Ne~14000
    T_AF = 148e3 / generation_time
    N_AF = 14474
    
    # 51kya:
    # non-AFR pops migrate OOA; bottlenecks to Ne~1800
    # migration between AFR occurs
    N_B = 1861
    T_SPLIT = 51e3 / generation_time
    m_AF_B = 15e-5
    
    # 23kya:
    # 2nd EUR bottlenecks to Ne~1000 & starts growing with rate 0.307%
    # migration rate slows between AFR-EUR
    N_EU0 = 1032
    T_EU_B = 23e3 / generation_time
    m_AF_EU = 2.5e-5
    r_EU0 = 0.00307
    N_EU1 = N_EU0 / math.exp(-r_EU0 * T_EU_B)
    
    # 5.1kya:
    # explosive growth in both AFR & EUR
    # Fu 2013  
#     T_EG = 5.1e3 / generation_time
#     r_EU = 0.0195
#     r_AF = 0.0166
#     N_EU_start = N_EU1 / math.exp(-r_EU * T_EG)
#     m_EG = 0
#     N_AF_start = N_AF / math.exp(-r_AF * T_EG)

    # Chen 2015
    T_EG = 7.26e3 / generation_time 
    r_EU = 0.0149
    r_AF = 0.00735
    N_EU_start = N_EU1 / math.exp(-r_EU * T_EG)
    m_EG = 0
    N_AF_start = N_AF / math.exp(-r_AF * T_EG)

    # Gazave 2014
#     T_EG = 3.52e3 / generation_time 
#     r_EU = 0.034
#     r_AF = 0.00735
#     N_EU = N_EU1 / math.exp(-r_EU * T_EG)
#     m_EG = 0
#     N_AF1 = N_AF / math.exp(-r_AF * T_EG)
    
    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU initially.
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=n_afr, initial_size=N_AF_start, growth_rate=r_AF),
        msprime.PopulationConfiguration(
            sample_size=n_eur, initial_size=N_EU_start, growth_rate=r_EU)#,
    ]

    # up to 5.1kya, no migration
    migration_matrix = [
        [0, 0],
        [0, 0],
    ]
    
    demographic_events = [
        # at 5.1kya, change to slow growth rate in EUR & stop growth in AFR;
        # add migration rate
        msprime.MigrationRateChange(
            time=T_EG, rate=m_AF_EU, matrix_index=(0, 1)),
        msprime.MigrationRateChange(
            time=T_EG, rate=m_AF_EU, matrix_index=(1, 0)),
        msprime.PopulationParametersChange(
            time=T_EG, growth_rate=r_EU0, initial_size=N_EU1, population_id=1),
        msprime.PopulationParametersChange(
            time=T_EG, growth_rate=0, population_id=0),
        
        # at 23kya, EUR growth stops and migration rates increase
        msprime.MigrationRateChange(
            time=T_EU_B, rate=m_AF_B, matrix_index=(0, 1)),
        msprime.MigrationRateChange(
            time=T_EU_B, rate=m_AF_B, matrix_index=(1, 0)),
        msprime.PopulationParametersChange(
            time=T_EU_B, initial_size=N_EU0, growth_rate=0, population_id=1),
        
        # at 51kya, population B merges into AFR
        msprime.MassMigration(
            time=T_SPLIT, source=1, destination=0, proportion=1.0),
        msprime.PopulationParametersChange(
            time=T_SPLIT, initial_size=N_B, population_id=1),
        
        # At 148kya, instantaneous growth in AFR
        msprime.PopulationParametersChange(
            time=T_AF, initial_size=N_A, population_id=0)
    ]
    
    if(debug):
        # Use the demography debugger to print out the demographic history
        # that we have just described.
        dd = msprime.DemographyDebugger(
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,
            demographic_events=demographic_events)
        dd.print_history()
    else:
        sim = msprime.simulate(population_configurations=population_configurations,
                               migration_matrix=migration_matrix, 
                               demographic_events=demographic_events,
                               mutation_rate=mu, 
                               recombination_rate=phi, 
                               length=length,
                               random_seed=30)
        return sim

###############################################################################
# Chen et al., 2015
###############################################################################
def chen_model(mu=1.5e-8, phi=0, length=1e4, n_afr=0, n_eur=0, debug=False):
    # First we set out the maximum likelihood values of the various parameters
    # given in Table 1.
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    
    # 220kya:
    # African population constant with Ne~7300
    N_A = 7310
    
    # 148kya:
    # instantaneous growth to Ne~14000
    T_AF = 148e3 / generation_time
    N_AF = 14474
    
    N6_EU = 13143
    
    # 118kya:
    # non-AFR pops migrate OOA; bottlenecks to Ne~1800
    # migration between AFR occurs
    N_B = 1861
    T5 = 118e3 / generation_time
    T4 = T5
    m_AF_B = 15e-5
    N5_EU = 62
    N4_EU = N6_EU
    
    # 18kya:
    # 2nd EUR bottlenecks to Ne~1000 & starts growing with rate 0.307%
    # migration rate slows between AFR-EUR
    N_EU0 = 1032
    T3 = 18e3 / generation_time
    T2 = T3
#     m_AF_EU = 2.5e-5
    r_EU0 = 0 # 0.00307
#     N2_EU = 15829 # N_EU0 / math.exp(-r_EU0 * T_EU_B)
    N2_EU = 16178
    N2_AF = 26682
    N3_EU = 2020
    # 5.1kya:
    # explosive growth in both AFR & EUR
    # Fu 2013  
#     T_EG = 5.1e3 / generation_time
#     r_EU = 0.0195
#     r_AF = 0.0166
#     N_EU_start = N_EU1 / math.exp(-r_EU * T_EG)
#     m_EG = 0
#     N_AF_start = N_AF / math.exp(-r_AF * T_EG)

    # Chen 2015
#     T1_EU = 7.26e3 / generation_time 
    T1_EU = 4.95e3 / generation_time
    T1_AF = 10.01e3 / generation_time
#     r_EU = 0.0149
    r_EU = 0.022
    r_AF = 0.00735
#     N1_EU = 1.2e6 # N_EU1 / math.exp(-r_EU * T_EG)
    N1_EU = 1.261e6
    m_EG = 0
    N1_AF = 5.062e5 # N_AF / math.exp(-r_AF * T_EG)
    
    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU initially.
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=n_afr, initial_size=N1_AF, growth_rate=r_AF),
        msprime.PopulationConfiguration(
            sample_size=n_eur, initial_size=N1_EU, growth_rate=r_EU)#,
    ]

    # up to 5.1kya, no migration
    migration_matrix = [
        [0, 0],
        [0, 0],
    ]
    
    demographic_events = [
        # at 5.1kya, change to slow growth rate in EUR & stop growth in AFR;
        # add migration rate
#         msprime.MigrationRateChange(
#             time=T_EG, rate=m_AF_EU, matrix_index=(0, 1)),
#         msprime.MigrationRateChange(
#             time=T_EG, rate=m_AF_EU, matrix_index=(1, 0)),
        msprime.PopulationParametersChange(
            time=T1_EU, growth_rate=0, initial_size=N2_EU, population_id=1),
        msprime.PopulationParametersChange(
            time=T1_AF, growth_rate=0, initial_size=N2_AF, population_id=0),
        
        # at 18kya, bottleneck + instantaneous recovery to smaller Ne
        msprime.PopulationParametersChange(
            time=T2, initial_size=N3_EU, population_id=1),
        msprime.PopulationParametersChange(
            time=T3, initial_size=N4_EU, population_id=1),
        
        # at 118kya, bottleneck + instantaneous recovery to same Ne
        msprime.PopulationParametersChange(
            time=T4, initial_size=N5_EU, population_id=1),
        msprime.PopulationParametersChange(
            time=T5, initial_size=N6_EU, population_id=1),
#         msprime.PopulationParametersChange(
#             time=T6, initial_size=N4_EU, population_id=1),
        
        
#         msprime.MigrationRateChange(
#             time=T_EU_B, rate=m_AF_B, matrix_index=(0, 1)),
#         msprime.MigrationRateChange(
#             time=T_EU_B, rate=m_AF_B, matrix_index=(1, 0)),
#         msprime.PopulationParametersChange(
#             time=T2, initial_size=N_EU0, growth_rate=0, population_id=1),
        
        # at 51kya, population B merges into AFR
#         msprime.MassMigration(
#             time=T_SPLIT, source=1, destination=0, proportion=1.0),
#         msprime.PopulationParametersChange(
#             time=T_SPLIT, initial_size=N_B, population_id=1),
        
        # At 148kya, instantaneous growth in AFR
        msprime.PopulationParametersChange(
            time=T_AF, initial_size=N_A, population_id=0)
    ]
    
    if(debug):
        # Use the demography debugger to print out the demographic history
        # that we have just described.
        dd = msprime.DemographyDebugger(
            population_configurations=population_configurations,
            migration_matrix=migration_matrix,
            demographic_events=demographic_events)
        dd.print_history()
    else:
        sim = msprime.simulate(population_configurations=population_configurations,
                               migration_matrix=migration_matrix, 
                               demographic_events=demographic_events,
                               mutation_rate=mu, 
                               recombination_rate=phi, 
                               length=length,
                               random_seed=30)
        return sim

