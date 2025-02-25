"""Command-line interface (CLI) entry point and complete reference of qa functionality."""

import click

def welcome():
    """Print first to welcome the user while it waits to load the modules"""
    print("\n")
    print("              ╔═══════════════════════════════╗           ")
    print("              ║    _________  ________        ║           ")
    print("              ║   |\   ___  \|\   __  \       ║           ")
    print("              ║   \ \  \ |\  \ \  \|\  \      ║           ")
    print("              ║    \ \  \  \  \ \   __  \     ║           ")
    print("              ║     \ \  \__\  \ \  \ \  \    ║           ")
    print("              ║      \ \______  \ \__\ \__\   ║           ")
    print("              ║       \|____| \__\|__|\|__|   ║           ")
    print("              ║              \|__|            ║           ")
    print("              ║                               ║           ")
    print("              ║       QUANTUMALLOSTERY        ║           ")
    print("              ║  [quantumallostery.rtfd.io]   ║           ")
    print("              ╚═══════════════╗╔══════════════╝           ")
    print("                      ╔═══════╝╚═══════╗                  ")
    print("                      ║ THE KULIK LAB  ║                  ")
    print("                      ╚═══════╗╔═══════╝                  ")
    print("  ╔═══════════════════════════╝╚═══════════════════════╗ ")
    print("  ║   Code: github.com/davidkastner/quantumAllostery   ║ ")
    print("  ║   Docs: quantumallostery.readthedocs.io            ║ ")
    print("  ║      - Usage: qa --help                            ║ ")
    print("  ╚════════════════════════════════════════════════════╝ \n")

# Greet the user with the welcome screen
welcome()

@click.command()
@click.option("--combine_restarts", "-a", is_flag=True, help="Combines restarts within a single replicate",)
@click.option("--combine_replicates","-b",is_flag=True,help="Combines combined replicates trajectories",)
@click.option("--simple_xyz_combine", "-c", is_flag=True, help="Combines xyz files")
@click.option("--xyz2pdb", "-d", is_flag=True, help="Converts an xyz to a pdb")
@click.option("--clean_frames", "-e", is_flag=True, help="Cleans incomplete frames.")
@click.option("--charge_coupling_plot","-f",is_flag=True,help="Charge coupling between two residues plot",)
@click.option("--find_stalled", "-g", is_flag=True, help="Find TeraChem jobs stalled")
@click.help_option("--help", "-h", is_flag=True, help="List of all quantumAllostery key words")
@click.option("--get_heatmap", "-i", is_flag=True, help="Heat map of amino acid correlations")
@click.option("--cpptraj_covars", "-j", is_flag=True, help="Use CPPTraj to calculate covariance")
@click.option("--charge_matrix_analysis","-k",is_flag=True,help="Create a matrix of charge couplings",)
@click.option("--clean_qm", "-l", is_flag=True, help="Cleans QM single point jobs")
@click.option("--combine_qm_charges","-m",is_flag=True,help="Combine charge data across QM single points",)
@click.option("--combine_rep_qm_charges","-ma",is_flag=True,help="Combine charge data across QM replicates",)
@click.option("--predict", "-n",is_flag=True,help="Uses ML models to predict key residues",)
@click.option("--multiwfn_charges","-o",is_flag=True,help="Calculates charge schemes from Multiwfn",)
@click.argument("multiwfn_charge_args", nargs=4, type=int, required=False)
@click.option("--calc_esp", "-p", is_flag=True, help="Calculates ESP from Multiwfn output")
@click.option("--check_esp_failed", "-q", is_flag=True, help="Checks for unfinished ESP jobs")
@click.option("--plot_esp", "-r", is_flag=True, help="Plot the ESP of each scheme and component")
@click.option("--combine_sp_xyz", "-s", is_flag=True, help="Combine single point xyz's")
@click.option("--plot_heme_distortion","-t",is_flag=True,help="Plots heme distortion across replicates",)
@click.option("--td_coupling", "-u", is_flag=True, help="Time-dependent charge coupling")
@click.option("--centroid_distance", "-v", is_flag=True, help="Distance between two centroids")
@click.option("--distance_esp_plot", "-w", is_flag=True, help="Plot distance vs. ESP")
@click.option("--kde_dist_esp_plot", "-x", is_flag=True, help="Plot distance vs. ESP as a KDE")
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
    combine_rep_qm_charges,
    predict,
    multiwfn_charges,
    multiwfn_charge_args,
    calc_esp,
    check_esp_failed,
    plot_esp,
    combine_sp_xyz,
    plot_heme_distortion,
    td_coupling,
    centroid_distance,
    distance_esp_plot,
    simple_xyz_combine,
    kde_dist_esp_plot,
):

    if combine_restarts:
        click.echo("> Combine restarts:")
        click.echo("> Loading...")
        import qa.process
        import qa.manage

        compute_replicates = input(
            "> Would you like this performed across replicates (y/n)? "
        )
        atom_count = qa.process.get_atom_count()

        if compute_replicates == "y":
            qa.manage.run_all_replicates(
                lambda: qa.process.combine_restarts(atom_count)
            )
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
            pdb_template = "template.pdb"
            orig_xyz_name = "mc6sa_geometry.xyz"
            new_pdb_name = "mc6sa_geometry.pdb"
            qa.process.xyz2pdb_traj(orig_xyz_name, new_pdb_name, pdb_template)
        elif xyz_type == "f":
            # You can process multiple frames at once
            qa.process.xyz2pdb(
                ["0.xyz", "10000.xyz", "20000.xyz", "30000.xyz", "39900.xyz"]
            )

    elif clean_frames:
        click.echo("> Clean frames with problems:")
        click.echo("> Loading...")
        import qa.process

        qa.process.remove_incomplete_xyz()

    elif charge_coupling_plot:
        click.echo("> Generate a charge coupling plot for two residues:")
        click.echo("> Run this from the replicate you are interested in.")
        click.echo("> Loading...")
        import qa.process
        import qa.analyze
        import qa.plot
        import qa.manage

        mimochrome = input("> Which mimochrome (mc6, mc6s, mc6sa, mc7)? ")
        res_x = input("> What is the first residue (Asp1)? ")
        res_y = input("> What is the second residue (Gly2)? ")

        # Create a dictionary to map the pair of residues to their ranges
        if mimochrome == "mc6":
            residue_ranges = {
                ("Asp18", "Gln20"): [[-1.16, -0.52], [-0.35, 0.29]], # repl-1
                ("Asp18", "Gln21"): [[-1.11, -0.45], [-0.41, 0.25]], # repl-8
                ("Lys12", "Gln9"): [[0.71, 1.22], [-0.26, 0.25]], # repl-1
                ("Glu19", "Arg11"): [[-1.05, -0.55], [0.58, 1.08]], # repl-1
                ("Ser23", "Hm129"): [[-0.20, 0.41], [-1.21, -0.61]], # repl-5
                ("Arg27", "Ser23"): [[0.64, 1.28], [-0.28, 0.36]], # repl-5
            }
        elif mimochrome == "mc6s":
            residue_ranges = {
                ("Asp18", "Gln20"): [[-1.16, -0.52], [-0.35, 0.29]], # repl-1
                ("Asp18", "Gln21"): [[-1.11, -0.45], [-0.41, 0.25]], # repl-5
                ("Lys12","Gln9"): [[0.71, 1.22], [-0.26, 0.25]], # repl-5
                ("Glu19", "Arg11"): [[-1.05, -0.55], [0.58, 1.08]], # repl-1
                ("Ser23", "Hm129"): [[-0.20, 0.41], [-1.21, -0.61]], # repl-3
                ("Arg27", "Ser23"): [[0.64, 1.28], [-0.28, 0.36]], # repl-7
            }
        elif mimochrome == "mc6sa":
            residue_ranges = {
                ("Asp18", "Aib20"): [[-1.16, -0.52], [-0.35, 0.29]], # repl-3
                ("Asp18", "Gln21"): [[-1.11, -0.45], [-0.41, 0.25]], # repl-1
                ("Lys12", "Gln9"): [[0.71, 1.22], [-0.26, 0.25]], # repl-4
                ("Glu19", "Arg11"): [[-1.05, -0.55], [0.58, 1.08]], # repl-1
                ("Aib23", "Hm129"): [[-0.20, 0.41], [-1.21, -0.61]], # repl-4   
                ("Arg27", "Aib23"): [[0.64, 1.28], [-0.28, 0.36]], # repl-5
            }

        # Retrieve the range from the dictionary
        selected_range = residue_ranges.get((res_x, res_y))
        if selected_range:
            df = qa.analyze.get_joint_qres(res_x, res_y, selected_range)
        else:
            print(f"No range found for the pair: {res_x}, {res_y}")
        df.to_csv(f"{res_x}{res_y}.csv")

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
        data = input("> What file would you like to plot (default = chargematbb.csv)? ")
        data = data if data else "chargematbb.csv"

        out_file = input("> What would you like to name the output file (default = coupling)? ")
        out_file = out_file if out_file else "coupling"

        # delete = [0,15,16,27,28,29]
        delete = [0,15,16,27]
        mimochrome = input("> Which mimochrome (e.g., mc6, mc6s, mc6sa, mc7)? ")

        low_input = input("> What is your low value (default = -0.4) ")
        low = float(low_input) if low_input else -0.4

        high_input = input("> What is your high value (default = 0.4) ")
        high = float(high_input) if high_input else 0.4

        mimochrome_mapping = {"mc6": mc6, "mc6s": mc6s, "mc6sa": mc6sa, "mc7": mc7}
        if mimochrome in mimochrome_mapping:
            qa.plot.heatmap(data, mimochrome_mapping[mimochrome], delete = delete, out_file=out_file, v=[low, high])
        else:
            click.echo(f"> Unknown mimochrome: {mimochrome}")

    elif cpptraj_covars:
        click.echo("> Generate geometric covariance using CPPTraj:")
        click.echo("> Loading...")
        import qa.manage
        import qa.plot
        import qa.analyze
        compute_replicates = input("> Would you like this performed across replicates (y/n)? ")
        delete = [[0, 15, 16, 27, 28, 29], []]

        if compute_replicates == "y":
            qa.manage.run_all_replicates(
                lambda: qa.analyze.cpptraj_covars(delete, recompute=False)
            )
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

        compute_replicates = input(
            "> Would you like this performed across replicates (y/n)? "
        )
        delete = []
        recompute = True

        # Perform the charge analysis and generate matrix files
        if compute_replicates == "y":
            qa.manage.run_all_replicates(
                lambda: qa.analyze.charge_matrix_analysis(delete, recompute=recompute)
            )
        elif compute_replicates == "n":
            qa.analyze.charge_matrix_analysis(delete, recompute=recompute)
        else:
            print(f"> {compute_replicates} is not a valid response.")

    elif clean_qm:
        click.echo("> Cleaning QM single point jobs:")
        click.echo("> Loading...")
        import qa.manage
        import qa.process

        qa.process.clean_qm_jobs(0, 39901, 100)

    elif combine_qm_charges:
        click.echo("> Combining the QM charge data across single points:")
        click.echo("> Loading...")
        import qa.manage
        import qa.process

        compute_replicates = input(
            "> Would you like this performed across replicates (y/n)? "
        )
        if compute_replicates == "n":
            qa.process.combine_qm_replicates()
        elif compute_replicates == "y":
            qa.process.combine_qm_charges(0, 39901, 100)
        else:
            print(f"> {compute_replicates} is not a valid response.")


    elif combine_rep_qm_charges:
        import qa.analyze
        qa.analyze.combine_qm_charges_replicates()
        

    elif predict:
        click.echo("> Making predictions:")
        click.echo("> Loading...")
        import qa.process
        import qa.predict
        import qa.plot

        charge_files = ["mc6.xls", "mc6s.xls", "mc6sa.xls"]
        templates = ["mc6.pdb", "mc6s.pdb", "mc6sa.pdb"]
        mutations = [0, 2, 15, 16, 19, 22, 27]  # Which amino acids to remove
        models = ["RF", "MLP"]
        n_frames = 1  # Get every nth frame to cut down the training time
        shuffle = False  # Shuffle the data to get a better estimate of the error

        charges_df, labels_df = qa.predict.create_combined_csv(
            charge_files, templates, mutations
        )
        charges_mat, labels_mat = qa.predict.data_processing(
            charges_df, labels_df, n_frames=n_frames
        )

        if shuffle:
            charges_mat, labels_mat = qa.predict.shuffle_data(charges_mat, labels_mat)

        qa.predict.run_ml(charges_mat, labels_mat, models=models, recompute=True)
        # Control whether you want a by atom or by residue analysis with by_atom
        qa.plot.plot_feature_importance(models, templates[0], mutations, by_atom=False)

    elif multiwfn_charges:
        # I had to pass multiwfn_charge_args as an arg although it was ugly
        click.echo("> Computed charge schemes with Multiwfn:")
        click.echo("> Example usage: qa -o 1 0 39901 100")
        click.echo("> Load Multiwfn: module load multiwfn/noGUI_3.7")
        click.echo("> Loading...")
        
        import qa.analyze
        import qa.manage

        #replicate = 1, first = 0, last = 39901, step = 100
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

        which_mimochrome = input("   > Which mimochrome would you like this computed for (e.g., mc6, mc6s, mc6sa)? ")

        if which_mimochrome == "mc6":
            components = {
            "all": "1-487",
            "lower": "1-252",
            "upper": "253-424",
            "lower-his": "1-86,104-252",
            "heme": "425-486",
            "his": "87-103",
            "asp-glu": "259-285",
            "asp-glu-arg": "259-285,398-421"}
        elif which_mimochrome == "mc6s":
            components = {
            "all": "1-491", 
            "lower": "1-256",
            "upper": "257-428",
            "lower-his": "1-90,108-256",
            "heme": "429-490",
            "his": "91-107",
            "asp-glu": "263-289",
            "asp-glu-arg": "263-289,402-425"}
        elif which_mimochrome == "mc6sa":
            components = {
            "all": "1-489",
            "lower": "1-256",
            "upper": "257-426",
            "lower-his": "1-90,108-256",
            "heme": "427-488",
            "his": "91-107",
            "asp-glu": "263-289",
            "asp-glu-arg": "263-289,400-423",
            "asp-glu-aib-arg": "263-289,337-353,400-423"}
        qa.manage.collect_esp_components(components, first, last, step)

    elif check_esp_failed:
        click.echo("> Checking for unfinished ESP jobs:")
        click.echo("> Loading...")
        import qa.manage

        qa.manage.check_esp_failed()

    elif plot_esp:
        click.echo("> Plotting the ESP results:")
        click.echo("> Loading...")
        import qa.plot

        # schemes = ["ADCH_esp.csv", "Hirshfeld_esp.csv", "Mulliken_esp.csv", "Voronoi_esp.csv"]
        schemes = ["mc6.csv", "mc6*.csv", "mc6*a.csv"]
        qa.plot.esp_combined_barchart(schemes)

    elif combine_sp_xyz:
        click.echo("> Combine the xyz files from all the single points:")
        click.echo("> Loading...")
        import qa.process

        qa.process.combine_sp_xyz()

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
        click.echo("> Run from outside all replicates.")
        click.echo("> Loading...")
        import qa.process
        import qa.analyze
        import qa.plot
        import qa.manage

        res_x = input("> What is the first residue (Asp1)? ")
        res_y = input("> What is the second residue (Gly2)? ")
        replicate = int(input("> What replicate (e.g., 1,2,3)? "))
        qa.analyze.td_coupling(res_x, res_y, replicate_dir=f"{replicate}/")

    elif centroid_distance:
        click.echo("> Calculate the distance between two centroids:")
        click.echo("> Loading...")
        import qa.analyze

        print("\nExample mimochrome values:")
        print("MC6: {upper-helix: 253-424, Fe: 487, Asp18: 259-270, Arg27: 398-420}")
        print("MC6*: {upper-helix: 257-428, Fe: 491, Asp18: 263-274, Arg27: 402-425}")
        print("MC6*a: {upper-helix: 257-426, Fe: 489, Asp18: 263-274, Arg27: 400-423}\n")

        centroid1_atoms = input("What atoms are in component 1 (e.g., 487)? ")
        centroid2_atoms = input("What atoms are in component 2 (e.g., 253-424)? ")
        centroids = [[centroid1_atoms], [centroid2_atoms]]
        qa.analyze.centroid_distance(centroids)

    elif distance_esp_plot:
        click.echo("> Plot the distance between two componenets vs. their ESP:")
        click.echo("> Loading...")
        import qa.plot

        esp_choice = int(input("   > 0-all, 1-lower, 2-upper, 3-lower-his, 4-heme, 5-his? "))
        colormap = input("   > What colormap would you like (e.g., inferno, jet, turbo, viridis)? ")
        residues = input("   > Which residues (e.g., Asp18, Arg27, D-chain)? ")
        

        if residues == "Asp18":
            qa.plot.esp_dist_plot(
                esp_choice,
                xlim=(7.6, 23),
                ylim=(-230, 375),
                color_map=colormap,
                custom_colors=None,
            )
        if residues == "Arg27":
            qa.plot.esp_dist_plot(
                esp_choice,
                xlim=(3.1, 12.7),
                ylim=(-230, 375),
                color_map=colormap,
                custom_colors=None,
            )
        if residues == "D-chain":
            qa.plot.esp_dist_plot(
                esp_choice,
                xlim=(4.1, 13.5),
                ylim=(-225, 375),
                color_map=colormap,
                custom_colors=None,
            )

    elif kde_dist_esp_plot:
        click.echo("> Plot the distance between two componenets vs. their ESP:")
        click.echo("> Loading...")
        import qa.plot

        esp_choice = int(input("   > 0-all, 1-lower, 2-upper, 3-lower-his, 4-heme, 5-his? "))
        colormap = input("   > What colormap would you like (e.g., inferno, jet, turbo, viridis, Blues)? ")
        residues = input("   > Which residues (e.g., Asp18, Arg27, D-chain)? ")

        if residues == "Asp18":
            qa.plot.esp_kde_dist_plot(
                esp_choice,
                xlim=(7.6, 23),
                ylim=(-230, 375),
                color_map=colormap,
            )
        if residues == "Arg27":
            qa.plot.esp_kde_dist_plot(
                esp_choice,
                xlim=(3.1, 12.7),
                ylim=(-230, 375),
                color_map=colormap,
            )
        if residues == "D-chain":
            qa.plot.esp_kde_dist_plot(
                esp_choice,
                xlim=(4.1, 13.5),
                ylim=(-225, 375),
                color_map=colormap,
            )

    elif simple_xyz_combine:
        click.echo(
            "> Combine xyz files in the current working directory to create a trajectory:"
        )
        click.echo("> Loading...")
        import qa.process

        qa.process.simple_xyz_combine()

    else:
        click.echo("No functionality was requested.\nTry --help.")



mc6 = ["ACE","ASP","GLU","GLN","GLN","LEU","HIS","SER","GLN","LYS","ARG","LYS","ILE","THR","LEU","NHE","ACE","ASP","GLU","GLN","GLN","LEU","SER","SER","GLN","LYS","ARG","NHE","HEME","FE"]
mc6s = ["ACE","ASP","LEU","GLN","GLN","LEU","HIS","SER","GLN","LYS","ARG","LYS","ILE","THR","LEU","NHE","ACE","ASP","GLU","GLN","GLN","LEU","SER","SER","GLN","LYS","ARG","NHE","HEME","FE"]
mc6sa = ["ACE","ASP","LEU","GLN","GLN","LEU","HIS","SER","GLN","LYS","ARG","LYS","ILE","THR","LEU","NHE","ACE","ASP","GLU","AIB","GLN","LEU","SER","AIB","GLN","LYS","ARG","NHE","HEME","FE"]
mc7 = ["ACE","ASP","LEU","GLN","GLN","LEU","HIS","SER","GLN","LYS","ARG","LYS","ILE","THR","LEU","NHE","ACE","ASP","GLU","AIB","GLN","LEU","AIB","SER","GLN","LYS","ARG","NHE","HEME","FE"]

if __name__ == "__main__":
    # Run the command-line interface when this script is executed
    welcome()
    cli()
