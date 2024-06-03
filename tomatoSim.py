import msprime
import numpy
import os

os.chdir("/home/milesanderson/PhD/msprime/VCFs")
# %%
#Initial parameters
generation_time=5
T_SHG=56013/generation_time
T_SCG=57676/generation_time
T_CV4=70411/generation_time
T_CV3=105009/generation_time
T_CV1=106855/generation_time

TB_SHG=8466/generation_time
TB_SCG=22882/generation_time
TB_CV4=14749/generation_time
TB_CV3=18780/generation_time
TB_CV1=16988/generation_time
TB_CV2=3900

N_SCG=54554
N_SHG=6447
N_CV1=8213
N_CV2=3000
N_CV3=5403
N_CV4=1002

Nsplit_SCG=2043
Nsplit_SHG=10583
Nsplit_CV1=35154
Nsplit_CV2=25000
Nsplit_CV3=24532
Nsplit_CV4=92040

r_SCG=-(numpy.log(Nsplit_SCG/N_SCG)/TB_SCG)
r_SHG=-(numpy.log(Nsplit_SHG/N_SHG)/TB_SHG)
r_CV1=-(numpy.log(Nsplit_CV1/N_CV1)/TB_CV1)
r_CV2=-(numpy.log(Nsplit_CV2/N_CV2)/TB_CV2)
r_CV3=-(numpy.log(Nsplit_CV3/N_CV3)/TB_CV3)
r_CV4=-(numpy.log(Nsplit_CV4/N_CV4)/TB_CV4)


demography = msprime.Demography()
demography.add_population(name="SCG", initial_size=54554, growth_rate=r_SCG)
demography.add_population(name="SHG", initial_size=6447, growth_rate=r_SHG)
demography.add_population(name="CV1", initial_size=8213, growth_rate=r_CV1)
demography.add_population(name="CV2", initial_size=3000, initially_active=True, growth_rate=r_CV2)
demography.add_population(name="CV3", initial_size=5403, growth_rate=r_CV3)
demography.add_population(name="CV4", initial_size=1002, growth_rate=r_CV4)

demography.set_symmetric_migration_rate(["CV1", "CV2"], 3.79171037934382e-05)
demography.set_symmetric_migration_rate(["CV2", "CV3"], 5.00432022420341e-05)
demography.set_symmetric_migration_rate(["CV2", "CV4"], 6.12623886707099e-05)
demography.set_symmetric_migration_rate(["CV2", "SHG"], 8.95358767563108e-05)
demography.set_symmetric_migration_rate(["CV2", "SCG"], 6.38552836678913e-06)

demography.add_population_parameters_change(time=TB_SHG, population="SHG", initial_size=Nsplit_SHG, growth_rate=0)
demography.add_population_parameters_change(time=TB_SCG, population="SCG", initial_size=Nsplit_SCG, growth_rate=0)
demography.add_population_parameters_change(time=TB_CV1, population="CV1", initial_size=Nsplit_CV1, growth_rate=0)
demography.add_population_parameters_change(time=TB_CV2, population="CV2", initial_size=Nsplit_CV2, growth_rate=0)
demography.add_population_parameters_change(time=TB_CV3, population="CV3", initial_size=Nsplit_CV3, growth_rate=0)
demography.add_population_parameters_change(time=TB_CV4, population="CV4", initial_size=Nsplit_CV4, growth_rate=0)

demography.add_symmetric_migration_rate_change(time=TB_CV1, populations=["CV1","CV2"], rate=9.98239666437443e-05)
demography.add_symmetric_migration_rate_change(time=TB_CV3, populations=["CV2","CV3"], rate=0.000214905866891655)
demography.add_symmetric_migration_rate_change(time=TB_CV4, populations=["CV2","CV4"], rate=3.95890041650506e-05)
demography.add_symmetric_migration_rate_change(time=TB_SHG, populations=["CV2","SHG"], rate=0)
demography.add_symmetric_migration_rate_change(time=TB_SCG, populations=["CV2","SCG"], rate=0)

demography.add_population_split(time=T_SHG, derived=["SHG"], ancestral="CV2")
demography.add_population_split(time=T_SCG, derived=["SCG"], ancestral="CV2")
demography.add_population_split(time=T_CV4, derived=["CV4"], ancestral="CV2")
demography.add_population_split(time=T_CV3, derived=["CV3"], ancestral="CV2")
demography.add_population_split(time=T_CV1, derived=["CV1"], ancestral="CV2")

demography.sort_events()
demography

import demes
import demesdraw
import matplotlib.pyplot as plt

graph = msprime.Demography.to_demes(demography)
fig, ax = plt.subplots()  # use plt.rcParams["figure.figsize"]
demesdraw.tubes(graph, seed=1)
plt.show()
# %%
import multiprocessing

def printVCF(i):
    id = str(i)
    ts = msprime.sim_ancestry(
        {"CV1": 10, "CV2": 10, "CV3": 10, "CV4": 10, "SHG": 7, "SCG": 6}, 
        demography=demography, 
        ploidy=2,
        model="dtwf", 
        sequence_length=1000
        )
    mts = msprime.sim_mutations(ts, rate = 0.00000001)
    vcf = open(id + ".vcf", 'w')
    mts.write_vcf(vcf, contig_id=id)
    vcf.close
if __name__ == '__main__':
    with multiprocessing.Pool() as pool:
        for result in pool.map(printVCF, range(100)):
            print(result)

