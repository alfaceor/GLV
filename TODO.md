# TODO list

- [x] Aparently the unities in Jinju's model are quite peculiar the abundance is in the order to log10(CFU/g) and the mu{kprime,k} are in 10^(-9) log10(CFU/g/day). So that is a source of error on the stability.
- [x] I need to renormalize the parameters in Jinju's estimation to avoid numerical errors
- [ ] Now that I have a program to solve any lotka-volterra system I need to add external effects on this system, like the effect of antibiotics, based on the first principles of the antibiotics
- [ ] Implement the perturbation (antibiotic effect) in program
- [ ] Implement the calculation of beta-diversity, particularly Bray-Curtis and Jensen-Shannon to follow the dynamics
- [ ] Write Methods in paper
- [ ] How can I make something similar to a Mori projection for the effect of the other species?
- [ ] If I add a stochastic noise to each species ODE, should this be drawn from the same distribution? Some particles will feel a different thermal bath just for size differences so probably not the same spectral properties, rigth?
- [ ] Modify perturbation function to get antiobitic concentrationss
