all:glv_deter glv_0_0p1 glv_1_0p1 glv_all_0p1 glv_0_efs_1p0 glv_1_efs_1p0 glv_all_efs_1p0 glv_0_efs_1p0_ns_0_1p0 glv_1_efs_1p0_ns_1_1p0 glv_1_efs_1p0_nall_1p0 glv_0_efs_1p0_nall_1p0

NSPECIES=2
STEPS=2000
PSTEPS=20
F0=-4.0
STD=0.01
NTRIALS=5000

# Deterministic with NSPECIES without extft
glv_deter:
	glv-simulate run_name=glv_deter \
	n_trials=1 \
	n_species=$(NSPECIES) \
	steps=$(STEPS) \
	noise=none

# Random noise in on index=0 without extft
glv_0_0p1:
	glv-simulate run_name=glv_0_0p1 \
	n_trials=$(NTRIALS) \
	n_species=$(NSPECIES) \
	steps=$(STEPS) \
	noise=single noise.std=$(STD) noise.index=0

# Random noise in on index=1 without extf
glv_1_0p1:
	glv-simulate run_name=glv_1_0p1 \
	n_trials=$(NTRIALS) \
	n_species=$(NSPECIES) \
	steps=$(STEPS) \
	noise=single noise.std=$(STD) noise.index=1

# Random noise in all NSPECIES without extf
glv_all_0p1:
	glv-simulate run_name=glv_all_0p1 \
	n_trials=$(NTRIALS) \
	n_species=$(NSPECIES) \
	steps=$(STEPS) \
	noise=all noise.std=$(STD)

# Deterministic with external force on index=0
glv_0_efs_1p0:
	glv-simulate run_name=glv_0_efs \
	n_trials=1 \
	n_species=$(NSPECIES) \
	steps=$(STEPS) \
	extft=single extft.index=0 extft.f0=$(F0) extft.psteps=$(PSTEPS)

# Deterministic with external force on index=1
glv_1_efs_1p0:
	glv-simulate run_name=glv_1_efs \
	n_trials=1 \
	n_species=$(NSPECIES) \
	steps=$(STEPS) \
	extft=single extft.index=1 extft.f0=$(F0) extft.psteps=$(PSTEPS)

# Deterministic with external force on all NSPECIES
glv_all_efs_1p0:
	glv-simulate run_name=glv_all_efs \
	n_trials=1 \
	n_species=$(NSPECIES) \
	steps=$(STEPS) \
	extft=all extft.f0=$(F0) extft.psteps=$(PSTEPS)

# External force on index=0 with noise index=0
glv_0_efs_1p0_ns_0_1p0:
	glv-simulate run_name=glv_0_efs_ns \
	n_trials=$(NTRIALS) \
	n_species=$(NSPECIES) \
	steps=$(STEPS) \
	extft=single extft.index=0 extft.f0=$(F0) extft.psteps=$(PSTEPS) \
	noise=single noise.index=0 noise.std=$(STD)

# External force on index=1 with noise index=1
glv_1_efs_1p0_ns_1_1p0:
	glv-simulate run_name=glv_1_efs_ns \
	n_trials=$(NTRIALS) \
	n_species=$(NSPECIES) \
	steps=$(STEPS) \
	extft=single extft.index=1 extft.f0=$(F0) extft.psteps=$(PSTEPS) \
	noise=single noise.index=1 noise.std=$(STD)

# External force on index=0 with noise on all
glv_0_efs_1p0_nall_1p0:
	glv-simulate run_name=glv_0_efs_na \
	n_trials=$(NTRIALS) \
	n_species=$(NSPECIES) \
	steps=$(STEPS) \
	extft=single extft.index=0 extft.f0=$(F0) extft.psteps=$(PSTEPS) \
	noise=all noise.std=$(STD)

# External force on index=1 with noise on all
glv_1_efs_1p0_nall_1p0:
	glv-simulate run_name=glv_1_efs_na \
	n_trials=$(NTRIALS) \
	n_species=$(NSPECIES) \
	steps=$(STEPS) \
	extft=single extft.index=1 extft.f0=$(F0) extft.psteps=$(PSTEPS) \
	noise=all noise.std=$(STD)

perf:
	pytest benchmarks/test_perf_core.py \
	--benchmark-json=local_scaling.json \
	--benchmark-histogram=comparando
# 	pytest benchmarks/test_perf_core.py \
# 	--benchmark-group-by=param:n_species --benchmark-histogram=hardware_comparison

synopsis:
	quarto render Synopsis

clean:
	rm -f data/simul/glv_deter_n_2_traj.h5
	rm -f data/simul/glv_0_0p1_n_2_traj_SingleNoise.h5
	rm -f data/simul/glv_1_0p1_n_2_traj_SingleNoise.h5
	rm -f data/simul/glv_all_0p1_n_2_traj_AllNoise.h5
	rm -f data/simul/glv_0_efs_n_2_traj.h5
	rm -f data/simul/glv_1_efs_n_2_traj.h5
