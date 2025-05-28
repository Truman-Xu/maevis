import numpy as np
import torchani
import torch
from rdkit import Chem
from rdkit.Chem import AllChem

def map_species(mol, elements):
    species = []
    for atom in mol.GetAtoms():
        specie = elements.get(atom.GetSymbol())
        if specie is None:
            return
        species.append(specie)
    return np.array(species)

def get_mean_aev(mol3d, aev_comp, device='cpu'):
    mol_nh = AllChem.RemoveAllHs(mol3d)
    mol_species = map_species(mol_nh)
    if mol_species is None:
        return

    conf_coords = np.array(
        [conf.GetPositions() for conf in mol_nh.GetConformers()]
    )
    mol_species = torch.tensor(mol_species, device=device, dtype=torch.int)
    mol_species = mol_species.repeat(conf_coords.shape[0], 1)
    mol_coords = torch.tensor(conf_coords, device=device, dtype=torch.float)
    species, aev = aev_comp((mol_species, mol_coords))
    return species, aev
    # return aev.mean((0, 1))

def avg_aev_by_species(aev, species, n_species):
    bincount = torch.bincount(species, minlength=n_species)
    div = bincount.float()
    div[bincount == 0] = 1.0
    numerator = torch.zeros((n_species, aev.shape[-1]))
    numerator = numerator.index_add(0, species, aev)
    mean_aev = numerator/div.unsqueeze(-1)
    return mean_aev

def get_probe_combo_aev(probe_mols, aev_comp, element_map, n_species, device='cpu'):
    combo_coords = []
    combo_species = []
    for molh in probe_mols:
        mol_nh = AllChem.RemoveAllHs(molh)
        mol_species = map_species(mol_nh, element_map)
        conf_coords = np.array([conf.GetPositions() for conf in mol_nh.GetConformers()])
        center_id = np.abs(
            conf_coords - np.median(conf_coords, axis = 0)
        ).mean((1,2)).argmin()
        center_coords = conf_coords[center_id]
        combo_coords.append(center_coords)
        combo_species.append(mol_species)
    combo_coords = np.concatenate(combo_coords)
    combo_species = np.concatenate(combo_species)
    combo_coords = torch.tensor(combo_coords, device=device, dtype=torch.float).unsqueeze(0)
    combo_species = torch.tensor(combo_species, device=device, dtype=torch.int).unsqueeze(0)
    combo_species, aev = aev_comp((combo_species, combo_coords))
    mean_aev = avg_aev_by_species(aev.squeeze(0), combo_species.squeeze(0), n_species)
    return mean_aev

def query_membership(query_memb, lookup, weight_dict, n_bins, n_neighbors=20):
    nei_ids = {}
    for bin_num in reversed(range(n_bins)):
        aev_dims = torch.where(query_memb == bin_num)[0]
        candidates = []
        for aev_dim in aev_dims:
            candidates.extend(lookup[bin_num][aev_dim.item()])
        nei_ids[bin_num] = torch.bincount(
            torch.tensor(candidates, dtype=torch.int32)
        ).argsort().flip(0)[:n_neighbors]

    weighted_neighbors = {}
    for bin_num, nei_list in nei_ids.items():
        weight = weight_dict[bin_num]
        for nei_idx in nei_list:
            nei_idx = nei_idx.item()
            if nei_idx not in weighted_neighbors:
                weighted_neighbors[nei_idx] = weight
            else:
                weighted_neighbors[nei_idx] += weight
    
    close_neighbors = sorted(weighted_neighbors.keys(), key=lambda x: weighted_neighbors[x], reverse=True)
    close_neighbors = close_neighbors[:n_neighbors]
    return close_neighbors

def query_neighbors(mol, bin_bounds, lookup, aev_comp, weight_dict, n_neighbors=20):
    n_bins = len(bin_bounds)+1
    mean_aev = get_mean_aev(mol, aev_comp)
    mol_memb = torch.bucketize(mean_aev, bin_bounds).to(torch.uint8).squeeze(0)

    return query_membership(mol_memb, lookup, weight_dict, n_neighbors)

