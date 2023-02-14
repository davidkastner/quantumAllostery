"""Command-line interface (CLI) entry point."""

# Print first to welcome the user while it waits to load the modules
print("\n.-------------------------------------------.")
print("| WELCOME TO THE QUANTUM ALLOSTERY (QA) CLI |")
print(".-------------------------------------------.")
print("Default programmed actions for the quantumAllostery package.")
print("GitHub: https://github.com/davidkastner/quantumAllostery")
print("Documenation: https://quantumallostery.readthedocs.io\n")

import click

@click.command()
@click.option("--combine_restarts", "-a", is_flag=True, help="Combines restarts within a single replicate.")
@click.option("--combine_replicates", "-c", is_flag=True, help="Combines combined replicates trajectories.")
@click.option("--xyz2pdb", "-d", is_flag=True, help="Converts an xyz to a pdb.")
@click.option("--clean_frames", "-e", is_flag=True, help="Cleans incomplete frames.")
@click.option("--charge_coupling_plot", "-f", is_flag=True, help="Charge coupling between two residues plot.")
@click.option("--find_stalled", "-g", is_flag=True, help="Find TeraChem jobs stalled.")
@click.option("--get_heatmap", "-i", is_flag=True, help="Heat map of amino acid correlations.")
@click.option("--cpptraj_covars", "-j", is_flag=True, help="Use CPPTraj to calculate covariance.")
@click.option("--charge_matrix_analysis", "-k", is_flag=True, help="Create a matrix of charge couplings.")
@click.option("--clean_qm", "-l", is_flag=True, help="Cleans QM single point jobs.")
@click.option("--combine_qm_charges", "-m", is_flag=True, help="Combine charge data across QM single points.")
@click.help_option('--help', '-h', is_flag=True, help='Exiting quantumAllostery.')
def cli(
    combine_restarts,
    combine_replicates,
    xyz2pdb,
    clean_frames,
    charge_coupling_plot,
    find_stalled,
    get_heatmap,
    cpptraj_covars,
    charge_matrix_analysis,
    clean_qm,
    combine_qm_charges,
    ):
    """
    The overall command-line interface (CLI) entry point.
    The CLI interacts with the rest of the package.

    A complete reference of quantumAllostery functionality.
    This is advantagous because it quickly introduces so quantumAllostery.
    Specificaly, to the complete scope of available functionality.
    It also improves long-term maintainability and readability.

    """
    
    if combine_restarts:
        click.echo("> Combine restarts:")
        click.echo("> Loading...")
        import qa.process
        compute_replicates = input("> Would you like this performed across replicates (y/n)? ")
        atom_count = qa.process.get_atom_count()
        
        if compute_replicates == "y":
            qa.manage.run_all_replicates(lambda: qa.process.combine_restarts(atom_count))
        elif compute_replicates == "n":
            qa.process.combine_restarts(atom_count)
        else:
            print(f"> {compute_replicates} is not a valid response.")

    elif combine_replicates:
        click.echo("> Combine trajectories from multiple replicates:")
        click.echo("> Loading...")
        import qa.process
        qa.process.combine_replicates()
    
    elif xyz2pdb:
        click.echo("> Convert an xyz to a pdb trajectory:")
        click.echo("> Loading...")
        import qa.process
        qa.process.xyz2pdb_traj()
    
    elif clean_frames:
        click.echo("> Clean frames with problems:")
        click.echo("> Loading...")
        import qa.process
        qa.process.remove_incomplete_xyz()
    
    elif charge_coupling_plot:
        click.echo("> Generate a charge coupling plot for two residues:")
        click.echo("> Loading...")
        import qa.process
        import qa.analyze
        import qa.plot
        import qa.library
        res_x = input("> What is the first residue (Asp1)? ")
        res_y = input("> What is the second residue (Gly2)? ")
        charge_df = qa.analyze.get_joint_qres(res_x, res_y)

    elif find_stalled:
        click.echo("> Checking for stalled TeraChem jobs:")
        click.echo("> Loading...")
        import qa.manage
        qa.manage.find_stalled()
    
    elif get_heatmap:
        click.echo("> Create heatmap of charge correlations:")
        click.echo("> Loading...")
        import qa.plot
        import qa.library
        data_file = input("> What file would you like to plot? ")
        qa.plot.heatmap(data=data_file, protein=protein, out_file="matrix_geom.png")
    
    elif cpptraj_covars:
        click.echo("> Generate geometric covariance using CPPTraj:")
        click.echo("> Loading...")
        import qa.manage
        import qa.plot
        import qa.analyze

        compute_replicates = input("> Would you like this performed across replicates (y/n)? ")
        delete = [[0,15,16,27,28,29],[]]
        
        if compute_replicates == "y":
            qa.manage.run_all_replicates(lambda: qa.analyze.cpptraj_covars(delete, recompute=False))
        elif compute_replicates == "n":
            qa.analyze.cpptraj_covars(delete, recompute=False)
        else:
            print(f"> {compute_replicates} is not a valid response.")

    elif charge_matrix_analysis:
        click.echo("> Perform a complete charge matrix analysis:")
        click.echo("> Loading...")
        import qa.manage
        import qa.plot
        import qa.analyze

        compute_replicates = input("> Would you like this performed across replicates (y/n)? ")
        delete = [[],[]]
        
        # Perform the charge analysis and generate matrix files
        if compute_replicates == "y":
            qa.manage.run_all_replicates(lambda: qa.analyze.charge_matrix_analysis(delete, recompute=False))
        elif compute_replicates == "n":
            qa.analyze.charge_matrix_analysis(delete, recompute=False)
        else:
            print(f"> {compute_replicates} is not a valid response.")

    elif clean_qm:
        click.echo("> Cleaning QM single point jobs:")
        click.echo("> Loading...")
        import qa.manage
        import qa.process
        qa.process.clean_qm_jobs(0, 39900, 100)

    elif combine_qm_charges:
        click.echo("> Combining the QM charge data across single points:")
        click.echo("> Loading...")
        import qa.manage
        import qa.process
        qa.process.combine_qm_charges(0, 39900, 100)

    else:
        click.echo("No functionality was requested.\nTry --help.")


if __name__ == "__main__":
    # Run the command-line interface when this script is executed
    cli()
