Pipeline usage








Prediction_boltz_pdb:

Following the designpipline.ipynb
1. Before prediction
generate yaml file, boltzutils.fasta_to_boltz_yaml
2. Clean Pdb OR ADD Hydrogen on, using boltzutils.do_boltzpdb_clean
Especially, for Gnina docking.
3. Parse docking score - parse_gnina.py
4. Parse interface score using - parse_rosetta_score.py
5. Parse confidence_score from boltz - boltzutils.prase_boltz_confidence_results
Parse predicted pdb rmsd using - parse_boltzpdb_rmsd


Culling















eval_mut_comb:



  - do_mut_pipeline:


     - extract a pose, and relax it using pyrosetta.pose_from_pdb()
     - df = read_from cotimed
     - filtered_AAs (prob > 0.1), output: e.g. {'MET1': ['S', 'T'], 'LYS2': ['Q', 'R']}
     - do_filtered_aas()
        the output should be a json document recorded ddG/wt_dG/mutant_dG/mutant_seq




     go to the scan of 