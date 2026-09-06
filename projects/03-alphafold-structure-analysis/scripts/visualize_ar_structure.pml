# PyMOL Visualization Script: Human Androgen Receptor Structure & AlphaFold Analysis
# Usage in PyMOL: @visualize_ar_structure.pml

reinitialize

# Load experimental crystal structure (1E3G) and AlphaFold model
load ../data/1e3g_human_ar_lbd.pdb, ar_crystal
load ../data/af_p10275_human_ar.pdb, ar_alphafold

# Hide all default representations
hide everything, all

# Display crystal structure as cartoon
show cartoon, ar_crystal
color gray80, ar_crystal

# Display AlphaFold model as cartoon
show cartoon, ar_alphafold

# Color AlphaFold structure by pLDDT confidence (stored in b-factor column)
# Very High (>90): Deep Blue
color 0x0053D6, ar_alphafold and b > 90
# Confident (70-90): Light Blue / Cyan
color 0x65CBF3, ar_alphafold and b >= 70 and b <= 90
# Low (50-70): Yellow
color 0xFFDB13, ar_alphafold and b >= 50 and b < 70
# Very Low / Disordered (<50): Orange
color 0xFF7D45, ar_alphafold and b < 50

# Superimpose AlphaFold LBD (residues 676-919) onto experimental crystal structure
align ar_alphafold and resi 676-919 and name CA, ar_crystal and name CA

# Display DHT ligand in 1E3G
show sticks, resn DHT
color green, resn DHT and elem C
color red, resn DHT and elem O

# Highlight key active site pocket residues (Asn705, Thr877, Phe764, Met745)
select pocket_residues, (ar_crystal and resi 705+877+764+745+704+711)
show sticks, pocket_residues
color magenta, pocket_residues and elem C

# Ray-tracing settings
set ray_shadow, 1
set depth_cue, 1
set cartoon_transparency, 0.2, ar_alphafold
zoom resn DHT, 12

print("PyMOL visualization script executed successfully.")
