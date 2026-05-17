all:glv_deter glv_0_0p1 glv_1_0p1 glv_all_0p1 
STEPS=2000
# No file dependency just for now
glv_deter:
	glv-simulate run_name=glv_deter \
	n_species=2 \
	steps=$(STEPS) \
	noise=none

glv_0_0p1:
	glv-simulate run_name=glv_0_0p1 \
	n_species=2 \
	steps=$(STEPS) \
	noise=single noise.std=0.1 noise.index=0

glv_1_0p1:
	glv-simulate run_name=glv_1_0p1 \
	n_species=2 \
	steps=$(STEPS) \
	noise=single noise.std=0.1 noise.index=1

glv_all_0p1:
	glv-simulate run_name=glv_all_0p1 \
	n_species=2 \
	steps=$(STEPS) \
	noise=all noise.std=0.1

synopsis:
	quarto render Synopsis

clean:
	rm -f data/simul/glv_deter_n_2_traj.h5
	rm -f data/simul/glv_0_0p1_n_2_traj_SingleNoise.h5
	rm -f data/simul/glv_1_0p1_n_2_traj_SingleNoise.h5
