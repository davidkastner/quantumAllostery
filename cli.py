"""Command-line interface (CLI) entry point."""

# Print first to welcome the user while it waits to load the modules
print("\n.-------------------------------------------.")
print("| WELCOME TO THE QUANTUM ALLOSTERY (QA) CLI |")
print(".-------------------------------------------.")
print("Default programmed actions for the quantumAllostery package.")
print("GitHub: https://github.com/davidkastner/quantumAllostery")
print("Documenation: https://quantumallostery.readthedocs.io\n")

import sys
import click

@click.command()
@click.option("--combine_restarts", "-a", is_flag=True, help="Combines restarts within a single replicate.")
@click.option("--combine_restart_replicates", "-b", is_flag=True, help="Combines restarts across all replicates.")
@click.option("--combine_replicates", "-c", is_flag=True, help="Combines combined replicates trajectories.")
@click.option("--xyz2pdb", "-d", is_flag=True, help="Converts an xyz to a pdb.")
@click.option("--clean_incomplete_frames", "-e", is_flag=True, help="Cleans incomplete frames.")
@click.option("--residue_charge_coupling_plot", "-f", is_flag=True, help="Charge coupling between two residues plot.")
@click.option("--find_stalled", "-g", is_flag=True, help="Find TeraChem jobs stalled at the COSMO step.")
@click.help_option('--help', '-h', is_flag=True, help='Exiting quantumAllostery.')
def cli(
    combine_restarts,
    combine_restart_replicates,
    combine_replicates,
    xyz2pdb,
    clean_incomplete_frames,
    residue_charge_coupling_plot,
    find_stalled
    ):
    """
    The overall command-line interface (CLI) entry point.
    The CLI interacts with the rest of the package.

    It will prompt the user with multiple levels of options.
    This is advantagous as quantumAllostery is moderate in scope,
    and because it introduces the user to the available functionality.
    Improves long-term maintainability.

    """
    if combine_restarts:
        click.echo("> Combine restarts:")
        click.echo("> Loading...")
        import qa.process
        atom_count = qa.process.get_atom_count()
        qa.process.combine_restarts(atom_count)
    
    elif combine_restart_replicates:
        click.echo("> Combine restarts from multiple replicates:")
        click.echo("> Loading...")
        import qa.process
        atom_count = qa.process.get_atom_count()
        qa.process.run_all_replicates(lambda: qa.process.combine_restarts(atom_count))
    
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
    
    elif clean_incomplete_frames:
        click.echo("> Clean frames with problems:")
        click.echo("> Loading...")
        import qa.process
        qa.process.remove_incomplete_xyz()
    
    elif residue_charge_coupling_plot:
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
        import qa.process
        qa.process.find_stalled()
    
    else:
        click.echo("No functionality was requested.\nTry --help.")

if __name__ == "__main__":
    # Run the command-line interface when this script is executed
    cli()
