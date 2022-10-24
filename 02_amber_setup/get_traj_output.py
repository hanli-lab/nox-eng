import argparse

import numpy as np
import pandas as pd
import pytraj as pt


def parse_args() -> argparse.Namespace:

    parser = argparse.ArgumentParser()

    parser.add_argument("--traj", type=str, required=True, help="path to trajectory .nc")
    parser.add_argument( "--parm", type=str, required=True, help="path to topology .parm")
    parser.add_argument("--label", type=str, required=True, help="run identifier")

    return parser.parse_args()


def get_hydride_dist(traj_path: str, parm_path: str, label: str) -> pd.DataFrame:

    # load
    traj = pt.load(traj_path, parm_path)

    # strip
    traj = pt.strip(traj, ":WAT,Na+,Cl-")

    # autoimage and align
    traj = pt.autoimage(traj)
    rmsd = pt.rmsd(traj, mask="@CA")

    # calc hydride distance
    cofa = "@C4" if "mna" in traj else "@C4N"
    fad = "@N5"
    dist = pt.distance(traj, f"{cofa} {fad}", image=True)

    # save output
    df = pd.DataFrame(dist, columns=["dist"])
    df["label"] = label
    df.to_csv(f"./{label}_hydride.csv", index=False)

    return df


def get_loop_rmsd(traj_path: str, parm_path: str, label: str) -> pd.DataFrame:

    # load
    traj = pt.load(traj_path, parm_path)

    # reference to align onto (starting structure)
    template_path = f"./00_templates/{label}.pdb"
    template = pt.load(template_path)

    # align
    rmsd = pt.rmsd(traj, mask="@CA", ref=template, ref_mask="@CA")

    # loop ca rmsd
    loop_rmsd = pt.rmsd(traj, mask=":153-157,177-188,238-242@CA")
    df = pd.DataFrame({"loop_rmsd": loop_rmsd})
    df["label"] = label
    df.to_csv(f"./{label}_loop_rmsd.csv", index=False)

    return df


def get_dihedrals(traj_path: str, parm_path: str, label: str) -> pd.DataFrame:

    # load
    traj = pt.load(traj_path, parm_path)

    # strip
    traj = pt.strip(traj, ":WAT,Na+,Cl-")

    # autoimage and align
    traj = pt.autoimage(traj)
    rmsd = pt.rmsd(traj, mask="@CA")

    # dihedrals
    dihe = pt.multidihedral(traj, dihedral_types="phi psi", dtype="dataframe")
    dihe.to_csv(f"./{label}_dihe.csv", index=False)


if __name__ == "__main__":

    args = parse_args()
    get_hydride_dist(args.traj, args.parm, args.label)
    get_loop_rmsd(args.traj, args.parm, args.label)
    get_dihedrals(args.traj, args.parm, args.label)
