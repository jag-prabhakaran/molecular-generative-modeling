import os
import numpy as np
import pandas as pd
import scipy.sparse
import torch
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect as Morgan
from rdkit.Chem.QED import qed
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import Descriptors
from collections import Counter
from functools import partial


def get_mol(smiles_or_mol):
    """
    Loads SMILES/molecule into RDKit's object
    """
    if isinstance(smiles_or_mol, str):
        if len(smiles_or_mol) == 0:
            return None
        mol = Chem.MolFromSmiles(smiles_or_mol)
        if mol is None:
            return None
        try:
            Chem.SanitizeMol(mol)
        except ValueError:
            return None
        return mol
    return smiles_or_mol


def get_n_rings(mol):
    """
    Computes the number of rings in a molecule
    """
    return mol.GetRingInfo().NumRings()


def fragmenter(mol):
    """
    fragment mol using BRICS and return smiles list
    """
    fgs = AllChem.FragmentOnBRICSBonds(get_mol(mol))
    fgs_smi = Chem.MolToSmiles(fgs).split(".")
    return fgs_smi


def compute_scaffold(mol, min_rings=2):
    mol = get_mol(mol)
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    except (ValueError, RuntimeError):
        return None
    n_rings = get_n_rings(scaffold)
    scaffold_smiles = Chem.MolToSmiles(scaffold)
    if scaffold_smiles == "" or n_rings < min_rings:
        return None
    return scaffold_smiles


def average_agg_tanimoto(
    stock_vecs, gen_vecs, batch_size=5000, agg="max", device="cpu", p=1
):
    """
    For each molecule in gen_vecs finds closest molecule in stock_vecs.
    Returns average tanimoto score for between these molecules

    Parameters:
        stock_vecs: numpy array <n_vectors x dim>
        gen_vecs: numpy array <n_vectors' x dim>
        agg: max or mean
        p: power for averaging: (mean x^p)^(1/p)
    """
    assert agg in ["max", "mean"], "Can aggregate only max or mean"
    agg_tanimoto = np.zeros(len(gen_vecs))
    total = np.zeros(len(gen_vecs))
    for j in range(0, stock_vecs.shape[0], batch_size):
        x_stock = torch.tensor(stock_vecs[j : j + batch_size]).to(device).float()
        for i in range(0, gen_vecs.shape[0], batch_size):
            y_gen = torch.tensor(gen_vecs[i : i + batch_size]).to(device).float()
            y_gen = y_gen.transpose(0, 1)
            tp = torch.mm(x_stock, y_gen)
            jac = (
                (tp / (x_stock.sum(1, keepdim=True) + y_gen.sum(0, keepdim=True) - tp))
                .cpu()
                .numpy()
            )
            jac[np.isnan(jac)] = 1
            if p != 1:
                jac = jac**p
            if agg == "max":
                agg_tanimoto[i : i + y_gen.shape[1]] = np.maximum(
                    agg_tanimoto[i : i + y_gen.shape[1]], jac.max(0)
                )
            elif agg == "mean":
                agg_tanimoto[i : i + y_gen.shape[1]] += jac.sum(0)
                total[i : i + y_gen.shape[1]] += jac.shape[0]
    if agg == "mean":
        agg_tanimoto /= total
    if p != 1:
        agg_tanimoto = (agg_tanimoto) ** (1 / p)
    return np.mean(agg_tanimoto)


def cos_similarity(ref_counts, gen_counts):
    """
    Computes cosine similarity between
     dictionaries of form {name: count}. Non-present
     elements are considered zero:

     sim = <r, g> / ||r|| / ||g||
    """
    if len(ref_counts) == 0 or len(gen_counts) == 0:
        return np.nan
    keys = np.unique(list(ref_counts.keys()) + list(gen_counts.keys()))
    ref_vec = np.array([ref_counts.get(k, 0) for k in keys])
    gen_vec = np.array([gen_counts.get(k, 0) for k in keys])
    return 1 - cos_distance(ref_vec, gen_vec)


def fingerprint(
    smiles_or_mol,
    fp_type="maccs",
    dtype=None,
    morgan__r=2,
    morgan__n=1024,
    *args,
    **kwargs
):
    """
    Generates fingerprint for SMILES
    If smiles is invalid, returns None
    Returns numpy array of fingerprint bits

    Parameters:
        smiles: SMILES string
        type: type of fingerprint: [MACCS|morgan]
        dtype: if not None, specifies the dtype of returned array
    """
    fp_type = fp_type.lower()
    molecule = get_mol(smiles_or_mol, *args, **kwargs)
    if molecule is None:
        return None
    if fp_type == "maccs":
        keys = MACCSkeys.GenMACCSKeys(molecule)
        keys = np.array(keys.GetOnBits())
        fingerprint = np.zeros(166, dtype="uint8")
        if len(keys) != 0:
            fingerprint[keys - 1] = 1  # We drop 0-th key that is always zero
    elif fp_type == "morgan":
        fingerprint = np.asarray(
            Morgan(molecule, morgan__r, nBits=morgan__n), dtype="uint8"
        )
    else:
        raise ValueError("Unknown fingerprint type {}".format(fp_type))
    if dtype is not None:
        fingerprint = fingerprint.astype(dtype)
    return fingerprint


def fingerprints(mols, fp_type):
    fingerprints = []
    for mol in mols:
        fingerprints.append(fingerprint(mol, fp_type))
    return fingerprints


def compute_scaffolds(mols):
    scaffolds = []
    for mol in mols:
        scaffolds.append(compute_scaffold(mol))
    return scaffolds


def compute_fragments(mols):
    fragments = []
    for mol in mols:
        fragments.append(fragmenter(mol))
    return fragments


# def compute_scaffolds(mol_list, n_jobs=1, min_rings=2):
#     """
#     Extracts a scafold from a molecule in a form of a canonic SMILES
#     """
#     scaffolds = Counter()
#     scaffolds = Counter(partial(compute_scaffold, min_rings=min_rings), mol_list)
#     if None in scaffolds:
#         scaffolds.pop(None)
#     return scaffolds


class SNNMetric:
    """
    Computes average max similarities of gen SMILES to ref SMILES
    """

    def __init__(self, fp_type="morgan", **kwargs):
        self.fp_type = fp_type
        super().__init__(**kwargs)

    def precalc(self, mols):
        return {"fps": fingerprints(mols, fp_type=self.fp_type)}

    def metric(self, pref, pgen):
        return average_agg_tanimoto(pref["fps"], pgen["fps"], device=self.device)


class FragMetric:
    def precalc(self, mols):
        return {"frag": compute_fragments(mols)}

    def metric(self, pref, pgen):
        return cos_similarity(pref["frag"], pgen["frag"])


class ScafMetric:
    def precalc(self, mols):
        return {"scaf": compute_scaffolds(mols)}

    def metric(self, pref, pgen):
        return cos_similarity(pref["scaf"], pgen["scaf"])


def get_all_metrics(
    gen,
    k=None,
    n_jobs=1,
    device="cpu",
    batch_size=512,
    pool=None,
    test=None,
    test_scaffolds=None,
    ptest=None,
    ptest_scaffolds=None,
    train=None,
):
    """
    Computes all available metrics between test (scaffold test)
    and generated sets of SMILES.
    Parameters:
        gen: list of generated SMILES
        k: int or list with values for unique@k. Will calculate number of
            unique molecules in the first k molecules. Default [1000, 10000]
        n_jobs: number of workers for parallel processing
        device: 'cpu' or 'cuda:n', where n is GPU device number
        batch_size: batch size for FCD metric
        pool: optional multiprocessing pool to use for parallelization

        test (None or list): test SMILES. If None, will load
            a default test set
        test_scaffolds (None or list): scaffold test SMILES. If None, will
            load a default scaffold test set
        ptest (None or dict): precalculated statistics of the test set. If
            None, will load default test statistics. If you specified a custom
            test set, default test statistics will be ignored
        ptest_scaffolds (None or dict): precalculated statistics of the
            scaffold test set If None, will load default scaffold test
            statistics. If you specified a custom test set, default test
            statistics will be ignored
        train (None or list): train SMILES. If None, will load a default
            train set
    Available metrics:
        * %valid
        * %unique@k
        * Frechet ChemNet Distance (FCD)
        * Fragment similarity (Frag)
        * Scaffold similarity (Scaf)
        * Similarity to nearest neighbour (SNN)
        * Internal diversity (IntDiv)
        * Internal diversity 2: using square root of mean squared
            Tanimoto similarity (IntDiv2)
        * %passes filters (Filters)
        * Distribution difference for logP, SA, QED, weight
        * Novelty (molecules not present in train)
    """
    if test is None:
        if ptest is not None:
            raise ValueError(
                "You cannot specify custom test " "statistics for default test set"
            )
        test = get_dataset("test")
        ptest = get_statistics("test")

    # if test_scaffolds is None:
    #     if ptest_scaffolds is not None:
    #         raise ValueError(
    #             "You cannot specify custom scaffold test "
    #             "statistics for default scaffold test set")
    #     test_scaffolds = get_dataset('test_scaffolds')
    #     ptest_scaffolds = get_statistics('test_scaffolds')

    train = train or get_dataset("train")

    if k is None:
        k = [1000, 10000]
    disable_rdkit_log()
    metrics = {}
    close_pool = False
    # if pool is None:
    #     if n_jobs != 1:
    #         pool = Pool(n_jobs)
    #         close_pool = True
    #     else:
    #         pool = 1
    # metrics['valid'] = fraction_valid(gen, n_jobs=pool)
    # gen = remove_invalid(gen, canonize=True)
    # if not isinstance(k, (list, tuple)):
    #     k = [k]
    # for _k in k:
    #     metrics['unique@{}'.format(_k)] = fraction_unique(gen, _k, pool)

    # if ptest is None:
    #     ptest = compute_intermediate_statistics(test, n_jobs=n_jobs,
    #                                             device=device,
    #                                             batch_size=batch_size,
    #                                             pool=pool)
    # if test_scaffolds is not None and ptest_scaffolds is None:
    #     ptest_scaffolds = compute_intermediate_statistics(
    #         test_scaffolds, n_jobs=n_jobs,
    #         device=device, batch_size=batch_size,
    #         pool=pool
    #     )
    # mols = mapper(pool)(get_mol, gen)
    # kwargs = {'n_jobs': pool, 'device': device, 'batch_size': batch_size}
    # kwargs_fcd = {'n_jobs': n_jobs, 'device': device, 'batch_size': batch_size}
    metrics["FCD/Test"] = FCDMetric(**kwargs_fcd)(gen=gen, pref=ptest["FCD"])
    metrics["SNN/Test"] = SNNMetric(**kwargs)(gen=mols, pref=ptest["SNN"])
    metrics["Frag/Test"] = FragMetric(**kwargs)(gen=mols, pref=ptest["Frag"])
    metrics["Scaf/Test"] = ScafMetric(**kwargs)(gen=mols, pref=ptest["Scaf"])
    if ptest_scaffolds is not None:
        metrics["FCD/TestSF"] = FCDMetric(**kwargs_fcd)(
            gen=gen, pref=ptest_scaffolds["FCD"]
        )
        metrics["SNN/TestSF"] = SNNMetric(**kwargs)(
            gen=mols, pref=ptest_scaffolds["SNN"]
        )
        metrics["Frag/TestSF"] = FragMetric(**kwargs)(
            gen=mols, pref=ptest_scaffolds["Frag"]
        )
        metrics["Scaf/TestSF"] = ScafMetric(**kwargs)(
            gen=mols, pref=ptest_scaffolds["Scaf"]
        )

    metrics["IntDiv"] = internal_diversity(mols, pool, device=device)
    metrics["IntDiv2"] = internal_diversity(mols, pool, device=device, p=2)
    metrics["Filters"] = fraction_passes_filters(mols, pool)

    # Properties
    for name, func in [("logP", logP), ("SA", SA), ("QED", QED), ("weight", weight)]:
        metrics[name] = WassersteinMetric(func, **kwargs)(gen=mols, pref=ptest[name])

    if train is not None:
        metrics["Novelty"] = novelty(mols, train, pool)
    enable_rdkit_log()
    if close_pool:
        pool.close()
        pool.join()
    return metrics
