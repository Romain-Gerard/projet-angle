import MDAnalysis as mda

#we will load the trajectory using start.gro as our topology file
u=mda.Universe("../data/start.gro", "../data/md_OK_dt100.xtc")

#now we will create .gro files for each frame in the xtc
for ts in u.trajectory:
  frame_number = ts.frame #getting the frame index
  output_filename=f"frame_{frame_number:04d}.gro" #eg: frame_0001.gro

  with mda.Writer(output_filename, n_atoms = u.atoms.n_atoms) as w:
    w.write(u.atoms)

print(f"Saved: {output_filename}")
