import msprime

# %%
demography = msprime.Demography()
demography.add_population(name="SCG", initial_size=8790)
demography.add_population(name="SHG", initial_size=15890)
demography.add_population(name="CV1", initial_size=32142)
demography.add_population(name="CV2", initial_size=23751)
demography.add_population(name="CV3", initial_size=27009)
demography.add_population(name="CV4", initial_size=15490)
demography.add_population(name="CV23", initial_size=50000)
demography.add_population(name="CV123", initial_size=50000)
demography.add_population(name="CV123SHG", initial_size=50000)
demography.add_population(name="CV1234SHG", initial_size=50000)
demography.add_population(name="CV1234SHGSCG", initial_size=50000)
demography.add_population_split(time=5400, derived=["CV2", "CV3"], ancestral="CV23")
demography.add_population_split(time=9000, derived=["CV1", "CV23"], ancestral="CV123")
demography.add_population_split(time=10000, derived=["SHG", "CV123"], ancestral="CV123SHG")
demography.add_population_split(time=12000, derived=["CV4", "CV123SHG"], ancestral="CV1234SHG")
demography.add_population_split(time=22400, derived=["SCG", "CV1234SHG"], ancestral="CV1234SHGSCG")
demography

# %%
for i in range(1, 3):
    id = str(i)
    ts = msprime.sim_ancestry(
        {"CV1": 10, "CV2": 10, "CV3": 10, "CV4": 10, "SHG": 10, "SCG": 10}, 
        demography=demography, 
        model="dtwf", 
        sequence_length=1000
        )
    mts = msprime.sim_mutations(ts, rate = 0.00000001)
    vcf = open(id + ".vcf", 'w')
    mts.write_vcf(vcf, contig_id=id)
    vcf.close()
