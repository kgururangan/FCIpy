import numpy as np
import datetime
from fcipy.utilities import get_memory_usage

WHITESPACE = "  "

ITERATION_HEADER_FMT = "{:>10} {:>12} {:>14} {:>17} {:>19} {:>12}"
ITERATION_FMT = "{:>8} {:>17.10f} {:>17.10f} {:>17.10f} {:>15} {:>12}"

DAV_ITERATION_HEADER = ITERATION_HEADER_FMT.format(
    "Iter.", "Residuum", "E", "dE", "Wall time", "Memory"
)

def print_ci_amplitudes(system, det, coef, thresh):
    from fcipy.lib import ci
    noa = system.noccupied_alpha
    nob = system.noccupied_beta
    idx = np.argsort(np.abs(coef))
    print("\n   Largest CI Amplitudes")
    for n, i in enumerate(reversed(idx)):
        if abs(coef[i]) < thresh: continue
        occ = ci.ci.get_occupied(det[:, :, i], noa)
        oa = [str(p) + "a" for p in occ[:noa, 0]]
        ob = [str(p) + "b" for p in occ[:nob, 1]]
        print(f"   [{n + 1}]    {oa}    {ob}     {coef[i]}")
    return

def print_ci_amplitudes_to_file(file, system, det, coef, thresh):
    from fcipy.lib import ci
    noa = system.noccupied_alpha
    nob = system.noccupied_beta
    idx = np.argsort(np.abs(coef))
    with open(file, "w") as f:
        for n, i in enumerate(reversed(idx)):
            if abs(coef[i]) < thresh: continue
            occ = ci.ci.get_occupied(det[:, :, i], noa)
            oa = [str(p) + "a" for p in occ[:noa, 0]]
            ob = [str(p) + "b" for p in occ[:nob, 1]]
            f.write(f"   [{n + 1}]    {oa}    {ob}     {coef[i]}\n")
    return

def get_timestamp():
        return datetime.datetime.strptime(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "%Y-%m-%d %H:%M:%S")

def dav_calculation_summary(energy, det, coef, is_converged, istate, system, print_thresh):
    DATA_FMT = "{:<30} {:>20.8f}"
    if is_converged:
        convergence_label = 'converged'
    else:
        convergence_label = 'not converged'
    print("\n   CI Calculation Summary (%s) - Root %i" % (convergence_label, istate))
    print("  --------------------------------------------------")
    print(DATA_FMT.format("   Total energy", energy))
    print_ci_amplitudes(system, det, coef, print_thresh)
    print("")

def print_dav_iteration_header():
    print("\n", DAV_ITERATION_HEADER)
    print('    '+(len(DAV_ITERATION_HEADER)) * "-")

def print_dav_iteration(
    iteration_idx, energy, residuum, delta_energy, elapsed_time
):
    minutes, seconds = divmod(elapsed_time, 60)
    time_str = f"({minutes:.1f}m {seconds:.1f}s)"
    memory = f"{round(get_memory_usage(), 2)} MB"
    print(
        ITERATION_FMT.format(
            iteration_idx, residuum, energy, delta_energy, time_str, memory
        )
    )

def print_block_eomcc_iteration(
    iteration_idx, curr_size, omega, residuum, delta_energy, elapsed_time, state_index
):
    minutes, seconds = divmod(elapsed_time, 60)
    time_str = f"({minutes:.1f}m {seconds:.1f}s)"
    memory = f"{round(get_memory_usage(), 2)} MB"
    for j, istate in enumerate(state_index):
        if j == 0:
            print(
                ITERATION_FMT.format(
                    iteration_idx, residuum[j], omega[istate], delta_energy[j], time_str, memory
                )
            )
        else:
            print(
                ITERATION_FMT.format(
                    "", residuum[j], omega[istate], delta_energy[j], time_str
                )
            )
    print("      Current subspace size = ", curr_size)
    print("      ............................................................")

class SystemPrinter:
    def __init__(self, system):
        self.system = system

    def header(self):

        print(WHITESPACE, "System Information:")
        print(WHITESPACE, "----------------------------------------------------")
        print(WHITESPACE, "  Number of correlated electrons =", self.system.nelectrons)
        print(WHITESPACE, "  Number of correlated orbitals =", self.system.norbitals)
        print(WHITESPACE, "  Number of frozen orbitals =", self.system.nfrozen)
        print(
            WHITESPACE,
            "  Number of alpha occupied orbitals =",
            self.system.noccupied_alpha,
        )
        print(
            WHITESPACE,
            "  Number of alpha unoccupied orbitals =",
            self.system.nunoccupied_alpha,
        )
        print(
            WHITESPACE,
            "  Number of beta occupied orbitals =",
            self.system.noccupied_beta,
        )
        print(
            WHITESPACE,
            "  Number of beta unoccupied orbitals =",
            self.system.nunoccupied_beta,
        )
        print(WHITESPACE, "  Charge =", self.system.charge)
        print(WHITESPACE, "  Point group =", self.system.point_group)
        print(WHITESPACE, "  Symmetry of reference =", self.system.reference_symmetry)
        print(
            WHITESPACE, "  Spin multiplicity of reference =", self.system.multiplicity
        )
        print("")

        HEADER_FMT = "{:>10} {:>20} {:>13} {:>13}"
        MO_FMT = "{:>10} {:>20.6f} {:>13} {:>13.1f}"

        header = HEADER_FMT.format("MO #", "Energy (a.u.)", "Symmetry", "Occupation")
        print(header)
        print(len(header) * "-")
        for i in range(self.system.norbitals + self.system.nfrozen):
            print(
                MO_FMT.format(
                    i + 1,
                    self.system.mo_energies[i],
                    self.system.orbital_symmetries_all[i],
                    self.system.mo_occupation[i],
                )
            )
        print("")
        print(WHITESPACE, "Memory Usage =", get_memory_usage(), "MB")
        print(WHITESPACE, "Frozen Core Energy =", self.system.frozen_energy)
        print(WHITESPACE, "Nuclear Repulsion Energy =", self.system.nuclear_repulsion)
        print(WHITESPACE, "Reference Energy =", self.system.reference_energy)
        print("")
        return

