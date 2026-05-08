all:glv_0_0p1 glv_1_0p1

glv_0_0p1:
	glv-simulate run_name=glv_0_0p1 \
	n_species=2 \
	noise=single noise.std=0.1 noise.index=0

glv_1_0p1:
	glv-simulate run_name=glv_1_0p1 \
	n_species=2 \
	noise=single noise.std=0.1 noise.index=1
