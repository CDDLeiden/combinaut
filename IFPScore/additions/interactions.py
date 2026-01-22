# hbonds
interactions_essential_hbond = [
    "hbondd_LYS_311_A",
    "hbondd_PHE_312_A",
]
interactions_interesting_hbond = [
    "hbondd_ARG_138_A",
    "hbonda_LYS_237_A",
    "hbondd_TYR_315_A",
    "hbondd_THR_77_A",
    "hbonda_TYR_315_A",
    "hbonda_THR_77_A"
]

# pi stacking
interactions_essential_pistack = [
    "pistack_TYR_305_A",
]
interactions_interesting_pistack = [
    "pistack_TYR_315_A",
    "pistack_PHE_312_A"
]

# lipophilic
interactions_essential_lipophilic = [  # the upper lipohilic pocket
    # top binding site
    "hydroph_VAL_244_A",  # this one could be specific to CCR2 only (CCR5 has leucine)
    "hydroph_LEU_81_A",

    # aromatic residues
    "hydroph_PHE_312_A",
    "hydroph_TYR_315_A",
]
interactions_interesting_lipophilic = [  # probably will turn up if pi-stacking not at optimal angle
    "hydroph_TYR_305_A",
    "hydroph_THR_77_A"  # from Lisa's review
    
    # top binding site
    "hydroph_ILE_245_A",
    "hydroph_ALA_241_A",
    "hydroph_LEU_67_A",
    "hydroph_VAL_63_A",
    "hydroph_LEU_134_A",
]
# interactions_maybe_lipohilic = [
#     "hydroph_LYS_311_A",  # questionable, with only part of the lysine
#     "hydroph_LYS_237_A",  # questionable, with only part of the lysine
#     "hydroph_ARG_138_A"  # questionable, with only part of the arginine
# ]

# other
interactions_essential_other = [
    # halogen bonds
    "halogenbond_VAL_63_A",
    "hbondd_GLU_310_A"
]
interactions_interesting_other = [
    # pi-cation
    "pication_LYS_311_A",
    "hbondd_LYS_311_A",
    "saltbridge_GLU_310_A",
    "hbonda_GLU_310_A"
]

# essential polar interactions
required = interactions_essential_hbond + interactions_essential_pistack
essential = interactions_essential_lipophilic + interactions_essential_other
interesting = interactions_interesting_lipophilic + interactions_interesting_hbond + interactions_interesting_pistack + interactions_interesting_other
