import argparse
from collections import defaultdict
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import pyrosetta
from Bio.PDB.Polypeptide import one_to_three


def mutate(
    pdb: str, mut_list: List[str], n_struct: int = 1_000, repack_shell: float = 6.0
) -> None:

    """mutate template model and repack local environment"""

    # output filename
    name = f"{Path(pdb).stem}_mut"

    # read input
    pose = pyrosetta.pose_from_pdb(pdb)

    # setup scoring
    sfxn = pyrosetta.get_fa_scorefxn()
    soft_rep = pyrosetta.create_score_function("soft_rep")

    # task ops to repack around each mutation
    around = pyrosetta.rosetta.protocols.task_operations.DesignAroundOperation()
    around.allow_design(False)
    around.resnums_allow_design(False)
    around.repack_shell(repack_shell)

    # mutate
    for mutation in mut_list:

        # parse mutation string
        start, pos, final = mutation[0], int(mutation[1:-1]), mutation[-1]
        final = one_to_three(final)

        # make each individual mutation and set repack around
        mut = pyrosetta.rosetta.protocols.simple_moves.MutateResidue(pos, final)
        mut.apply(pose)
        around.include_residue(pos)

    # set taskfactory
    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
    tf.push_back(around)

    # repack
    enzrepackmin = pyrosetta.rosetta.protocols.enzdes.EnzRepackMinimize()
    enzrepackmin.set_scorefxn_minimize(sfxn)
    enzrepackmin.set_scorefxn_repack(soft_rep)
    enzrepackmin.set_min_sc(True)
    enzrepackmin.set_design(False)
    enzrepackmin.set_min_bb(True)
    enzrepackmin.task_factory(tf)

    # setup mc sampling
    mc = pyrosetta.MonteCarlo(pose, sfxn, 0.6)
    trial = pyrosetta.TrialMover(enzrepackmin, mc)
    rep = pyrosetta.RepeatMover(trial, 10)

    jd = pyrosetta.PyJobDistributor(name, n_struct, sfxn)
    while not jd.job_complete:

        # work on temp pose
        temp = pyrosetta.Pose()
        temp.assign(pose)

        # make moves and save final structure
        mc.reset(temp)
        rep.apply(temp)
        mc.recover_low(temp)
        jd.output_decoy(temp)


def parse_sc(sc: str) -> None:

    """rosetta scorefile to dataframe"""

    scores = defaultdict(list)

    with open(sc, "r") as f:
        # skip header
        f.readline()
        lines = f.readlines()

    # each output stored in row of label: value
    for line in lines:
        split = line.replace(":", "").split()

        labels = split[::2]
        vals = split[1::2]

        for label, val in zip(labels, vals):
            scores[label].append(val)

    # save dataframe
    df = pd.DataFrame.from_dict(scores).sort_values(by="total_score", ascending=False)
    df.to_csv(f"{sc.split('.')[0]}.csv", index=False)


def clean_pdb(pdb: str) -> None:

    """clean pdb to use in amber"""

    name = Path(pdb).stem

    with open(pdb, "r") as f:
        lines = f.readlines()

    res = []
    for line in lines:

        # discard the rosetta scores after TER
        if "TER" in line:
            break

        # skip hydrogens, tleap will add
        if line.startswith("ATOM") and line.split()[-1] == "H":
            continue

        res.append(line)

    out = "".join(res)
    Path(f"{name}_clean.pdb").write_text(out)


if __name__ == "__main__":

    pyrosetta.init()

    # examples, change the mutation and top scoring output
    TEMPLATE = "../00_raw_data/lp_nox_template.pdb"
    MUT_LIST = ["D177W", "G178E"]
    N_STRUCT = 5
    SC = "lp_nox_template_mut.fasc"
    TOP = "lp_nox_template_mut_2.pdb"

    mutate(TEMPLATE, MUT_LIST, N_STRUCT)
    parse_sc(SC)
    clean_pdb(TOP)
