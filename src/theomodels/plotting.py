import torch
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter


def plot_trajectories(
    traj: torch.Tensor,
    ax: plt.Axes,
    linestyle: str = "-",
    label: str | None = None,
) -> None:
    """Plot the trajectories of all species."""
    tt = torch.arange(traj.size(0))
    ax.plot(tt.cpu().numpy(), traj.cpu().numpy(), linestyle=linestyle, label=label)
    ax.set_xlabel("Time")
    ax.set_ylabel("Abundance")
    ax.set_title("Generalized Lotka-Volterra Simulation")


def plot_matrix_A(A: torch.Tensor, ax: plt.Axes) -> None:
    """Plot the interaction matrix."""
    A_np = A.cpu().numpy()
    n_i, n_j = A_np.shape
    im = ax.imshow(A_np, cmap="coolwarm", interpolation="none")
    ax.set_title("Interaction matrix A")
    ax.set_xlabel("Species j")
    ax.set_ylabel("Species i")
    cbar = plt.colorbar(im, ax=ax, orientation="vertical", shrink=0.8)
    cbar.formatter = FormatStrFormatter("%.2g")
    cbar.update_ticks()

    for i in range(n_i):
        for j in range(n_j):
            ax.text(j, i, f"{A_np[i, j]:.2g}", ha="center", va="center", color="black")


