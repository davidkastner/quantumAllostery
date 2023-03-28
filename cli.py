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
@click.option("--predict", "-p", is_flag=True, help="Uses simple ML models to predict key residues.")
@click.option("--multiwfn_charges", "-q", is_flag=True, help="Calculates charge schemes from Multiwfn.")
@click.argument("multiwfn_charge_args", nargs=4, type=int, required=False)
@click.option("--calc_esp", "-r", is_flag=True, help="Calculates ESP from Multiwfn output.")
@click.option("--check_esp_failed", "-s", is_flag=True, help="Checks for unfinished ESP jobs.")
@click.option("--plot_esp", "-t", is_flag=True, help="Plot the ESP of each scheme and component.")
@click.option("--combine_sp_xyz", "-u", is_flag=True, help="Combine single point xyz's.")
@click.option("--plot_heme_distortion", "-v", is_flag=True, help="Plots heme distortion across replicates.")
@click.option("--td_coupling", "-w", is_flag=True, help="Time-dependent charge coupling.")
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
    predict,
    multiwfn_charges,
    multiwfn_charge_args,
    calc_esp,
    check_esp_failed,
    plot_esp,
    combine_sp_xyz,
    plot_heme_distortion,
    td_coupling,
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
        xyz_type = input("> Process a trajectory (t) or frames (f)? ")
        if xyz_type == "t":
            
            qa.process.xyz2pdb_traj()
        elif xyz_type == "f":
            qa.process.xyz2pdb(["0.xyz","10000.xyz","20000.xyz","30000.xyz","39900.xyz"])
    
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
        import qa.manage
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
        import qa.manage
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
        recompute = True
        
        # Perform the charge analysis and generate matrix files
        if compute_replicates == "y":
            qa.manage.run_all_replicates(lambda: qa.analyze.charge_matrix_analysis(delete, recompute=recompute))
        elif compute_replicates == "n":
            qa.analyze.charge_matrix_analysis(delete, recompute=recompute)
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

        compute_replicates = input("> Would you like this performed across replicates (y/n)? ")
        if compute_replicates == "n":
            qa.process.combine_qm_replicates()
        elif compute_replicates == "y":
            qa.process.combine_qm_charges(0, 39900, 100)
        else:
            print(f"> {compute_replicates} is not a valid response.")

    elif predict:
        click.echo("> Making predictions:")
        click.echo("> Loading...")
        import qa.process
        import qa.predict
        import qa.plot

        charge_files = ["mc6.xls","mc6s.xls","mc6sa.xls"]
        templates = ["mc6.pdb","mc6s.pdb","mc6sa.pdb"]
        mutations = [0,2,15,16,19,22,27] # Which amino acids to remove
        models = ["RF", "MLP"]
        n_frames = 1 # Get every nth frame to cut down the training time
        shuffle = False # Shuffle the data to get a better estimate of the error

        charges_df, labels_df = qa.predict.create_combined_csv(charge_files, templates, mutations)
        charges_mat, labels_mat = qa.predict.data_processing(charges_df, labels_df, n_frames = n_frames)

        if shuffle:
            charges_mat, labels_mat = qa.predict.shuffle_data(charges_mat, labels_mat)
        
        qa.predict.run_ml(charges_mat, labels_mat, models=models, recompute=True)
        # Control whether you want a by atom or by residue analysis with by_atom
        qa.plot.plot_feature_importance(models, templates[0], mutations, by_atom = False)

    elif multiwfn_charges:
        # I had to pass multiwfn_charge_args as an arg although it was ugly
        click.echo("> Computed charge schemes with Multiwfn:")
        click.echo("> Loading...")
        import qa.analyze
        import qa.manage

        # replicate = 1, first = 0, last = 39901, step = 100
        replicate = str(multiwfn_charge_args[0])
        first = multiwfn_charge_args[1]
        last = multiwfn_charge_args[2]
        step = multiwfn_charge_args[3]
        get_charges = lambda: qa.analyze.calculate_charge_schemes()
        qa.manage.replicate_interval_submit(replicate, first, last, step, get_charges)

    elif calc_esp:
        click.echo("> Computed charge schemes with Multiwfn:")
        click.echo("> Loading...")
        import qa.analyze
        import qa.manage
        first = 0
        last = 39901
        step = 100
        qa.manage.collect_esp_components(first, last, step)

    elif check_esp_failed:
        click.echo("> Checking for unfinished ESP jobs:")
        click.echo("> Loading...")
        import qa.manage
        qa.manage.check_esp_failed()

    elif plot_esp:
        click.echo("> Plotting the ESP results:")
        click.echo("> Loading...")
        import qa.plot
        qa.plot.esp_combined_barchart()

    elif combine_sp_xyz:
        click.echo("> Combine the xyz files from all the single points:")
        click.echo("> Loading...")
        import qa.manage
        qa.manage.combine_sp_xyz()

    elif plot_heme_distortion:
        click.echo("> Computes the RMSD for a structure across the trajectory:")
        click.echo("> Loading...")

        import qa.analyze
        import qa.process
        import qa.plot

        # User enters the atom sets here
        ref_atoms = "1-20,23,25,27,29,31,35,38,42,46,50,53"
        mc6_atoms = "435-454,457,459,461,463,465,469,472,476,480,484,487"
        mc6s_atoms = "439-458,461,453,465,467,469,474,476,480,484,488,491"
        mc6sa_atoms = "437-456,459,461,463,465,467,471,474,478,482,486,489"

        # Convert the strings to lists
        atoms_list = [ref_atoms, mc6_atoms, mc6s_atoms, mc6sa_atoms]
        atoms_list = qa.process.string_to_list(atoms_list)

        # Compute the RMSDs
        rmsd_list = qa.analyze.get_rmsd(atoms_list[0], atoms_list[1:])
        qa.plot.plot_rmsd(rmsd_list, ["MC6", "MC6*", "MC6*a"])

    elif td_coupling:
        click.echo("> Plot the fluctuations of two amino acids by time:")
        click.echo("> Do yo need to update your replicate?")
        click.echo("> Loading...")
        import qa.process
        import qa.analyze
        import qa.plot
        import qa.manage
        res_x = input("> What is the first residue (Asp1)? ")
        res_y = input("> What is the second residue (Gly2)? ")
        charge_df = qa.analyze.td_coupling(res_x, res_y, replicate_dir="6")

    else:
        click.echo("No functionality was requested.\nTry --help.")


if __name__ == "__main__":
    # Run the command-line interface when this script is executed
    cli()
