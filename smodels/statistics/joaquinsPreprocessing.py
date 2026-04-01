import numpy as np
import torch
from scipy.stats import boxcox
from sklearn.preprocessing import QuantileTransformer


def preprocess_nLLs(nLLs, trafos=None):
    mean, std = 0.0, 1.0

    if trafos:
        for fn_str in trafos:
            print("Applying nLL transformation:", fn_str)

            if fn_str == "standardization":
                nLLs, mean, std = standardization(
                    nLLs, return_mean_std=True, clip=False
                )

            elif fn_str == "log_w_negatives":
                nLLs = log_with_negatives(nLLs)

            else:
                fn = get_fn(fn_str)
                nLLs = fn(nLLs, None)

    assert np.isfinite(nLLs).all()
    return nLLs, mean, std


def log_with_negatives(nLLs):
    """
    Signed log1p transform.
    Works for negative values and is exactly invertible.
    """
    return np.sign(nLLs) * np.log1p(np.abs(nLLs))

def undo_log_with_negatives(nLLs):
    """
    Exact inverse of log_with_negatives.
    """
    return np.sign(nLLs) * (np.expm1(np.abs(nLLs)))



def undo_preprocess_nLLs(nLLs, mean, std, trafos=None):
    if trafos:
        for fn_str in reversed(trafos):

            if fn_str == "standardization":
                nLLs = nLLs * std + mean

            elif fn_str == "log_w_negatives":
                nLLs = undo_log_with_negatives(nLLs)

            else:
                inv_fn = get_inv_fn(fn_str)
                nLLs = inv_fn(nLLs, None)

    assert np.isfinite(nLLs).all()
    return nLLs



def preprocess_features(
    features_raw,
    type_tokens=None,
    trafos=None,
    incl_fvs=True,
    mean=None,
    std=None,
    eps_std=1e-2,
    return_dict=False,
):
    assert np.isfinite(features_raw).all()
    mean, std = 0.0, 1.0

    feature_sets = {}
    if trafos:
        for t_name, t_fns in trafos.items():
            transformed_features = features_raw
            for fn_str in t_fns:
                if fn_str == "standardization":
                    transformed_features, mean, std = standardization(
                        transformed_features, return_mean_std=True, clip=False
                    )
                else:
                    fn = get_fn(fn_str)
                    transformed_features = fn(transformed_features, type_tokens)
            feature_sets[t_name] = transformed_features

    # if incl_fvs:
    #     feature_sets["fvs_raw"] = sort_features(features_raw, type_tokens, "E")
    #     feature_sets["fvs_raw"] /= np.max(
    #         np.abs(feature_sets["fvs_raw"][:, ::4])
    #     )  # normalize by max energy of features
    #     # feature_sets["fvs_raw"] = features_raw

    if return_dict:
        return feature_sets, mean, std
    else:
        return np.concatenate(list(feature_sets.values()), axis=1), mean, std


def standardization(features, return_mean_std=False, clip=True):
    """standardize features"""
    mean = features.mean(axis=0)
    std = features.std(axis=0)
    if clip:
        std = np.clip(std, a_min=1e-6, a_max=None)  # avoid std=0.
    features = (features - mean) / std
    assert np.isfinite(features).all()
    if return_mean_std:
        return features, mean, std
    else:
        return features


def compute_invariants(features, eps=1e-4, incl_diag_invariants=False, reshape=True):
    """compute Lorentz invariants out of four-vectors"""
    if reshape:
        ps = features.reshape(features.shape[0], features.shape[1] // 4, 4)
    else:
        ps = features

    # compute matrix of all inner products
    def inner_product(p1, p2):
        return p1[..., 0] * p2[..., 0] - (p1[..., 1:] * p2[..., 1:]).sum(axis=-1)

    if incl_diag_invariants:
        offset = 0
    else:
        offset = 1

    idxs = np.triu_indices(ps.shape[-2], k=offset)
    invariants = inner_product(ps[..., idxs[0], :], ps[..., idxs[1], :])
    if type(invariants) == np.ndarray:
        invariants = np.clip(invariants, a_min=eps, a_max=None)
    else:
        invariants = torch.clip(invariants, eps, None)
    return invariants


def sort_features(features, type_tokens, sort_key):
    """sort features according to some property"""
    ps = features.reshape(features.shape[0], features.shape[1] // 4, 4)
    particle_types = np.unique(type_tokens)
    sorted_ps = []
    for particle_type in particle_types:
        mask = type_tokens == particle_type
        # sort according to energy, which is the first entry of the four-vector
        match sort_key:
            case "E":
                sorted_indices = np.argsort(ps[:, mask][:, :, 0])
            case "pt":
                pts = np.linalg.norm(ps[:, mask][:, :, 1:3], axis=-1)
                sorted_indices = np.argsort(pts)
            case _:
                raise ValueError(f"Unknown sort key {sort_key}")
        sorted_ps.append(
            np.take_along_axis(ps[:, mask], sorted_indices[:, :, np.newaxis], axis=1)
        )
    return np.concatenate(sorted_ps, axis=1).reshape(features.shape)


def apply_boxcox(features):
    ps = features.transpose()
    res = np.zeros(ps.shape)
    for i in range(ps.shape[0]):
        try:
            res[i] = boxcox(ps[i])[0]
        except:
            res[i] = ps[i]
    return res.transpose()


def apply_quantile_transform(features):
    return QuantileTransformer(output_distribution="normal").fit_transform(features)


def get_fn(fn_str):
    # get functions from string
    match fn_str:
        case "log":
            return lambda p, t: np.log(p)
        case "exp":
            return lambda p, t: np.exp(p)
        case "sqrt":
            return lambda p, t: np.sqrt(p)
        case "inverse":
            return lambda p, t: 1 / p
        case "boxcox":
            return lambda p, t: apply_boxcox(p)
        case "quantile_transform":
            return lambda p, t: apply_quantile_transform(p)
        case "invs":
            return lambda p, t: compute_invariants(p)
        case "standardization":
            return lambda p, t: standardization(p)
        case "log_w_negatives":
            return lambda p, t: log_with_negatives(p)
        case "sort_E":
            return lambda p, t: sort_features(p, t, "E")
        case "sort_pt":
            return lambda p, t: sort_features(p, t, "pt")
        case _:
            raise ValueError(f"Unknown transformation function {fn_str}")


def get_inv_fn(fn_str):
    # get inverse functions from string
    match fn_str:
        case "log":
            return lambda p, t: np.exp(p)
        case "exp":
            return lambda p, t: np.log(p)
        case "sqrt":
            return lambda p, t: p**2
        case "inverse":
            return lambda p, t: 1 / p
        case _:
            raise ValueError(f"Unknown inverse transformation function {fn_str}")
