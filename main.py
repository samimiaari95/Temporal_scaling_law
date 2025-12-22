import sys
import argparse
import simulation_setup
import data_processing
import plotting

# ==============================================================================
# MAIN SCRIPT
# Orchestrates the workflow: setup, run, process, analyze, and plot
# ==============================================================================

def run_step(step):
    """Run a specific step or all steps."""
    if step in [1, "all"]:
        print("Running Step 1: Setup simulations...")
        simulation_setup.parflow_namelist_inex()
        simulation_setup.simulation_bash()
        simulation_setup.parflow_namelist_inex_fixed()
    
    if step in [2, "all"]:
        print("Step 2: Run simulations (assumed to be done externally)")
        print("ParFlow should be run using the generated bash script.")
    
    if step in [3, "all"]:
        print("Running Step 3: Process simulation results...")
        data_processing.steadystate_kinsol_infiltrationdrainage()
    
    if step in [4, "all"]:
        print("Running Step 4: Analyze data and plot results...")
        plotting.plot_pressure_profile()
        plotting.infexf_lambda_dependence()
        plotting.dexf_dependence()
        plotting.dinf_dependence()
        plotting.inf_solution_plots()
        plotting.exf_solution_plots()
        plotting.qfit()
        plotting.q_vs_slope()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run workflow steps for temporal scaling law analysis"
    )
    parser.add_argument(
        "--step",
        type=int,
        choices=[1, 2, 3, 4],
        default="all",
        nargs="?",
        const=1,
        help="Specify which step to run (1-4) or 'all' for all steps"
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Run all steps"
    )
    
    args = parser.parse_args()
    
    if args.all:
        run_step("all")
    else:
        run_step(args.step)