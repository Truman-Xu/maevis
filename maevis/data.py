
# default element mappings
elements = {
    'C': 0, 'N': 1, 'O': 2, 'S': 3, 
    'P': 4, 'F': 5, 'Cl': 5, 'Br': 5, 'I': 5
}
element_map = {0: 'C', 1: 'N', 2: 'O', 3: 'S', 4: 'P', 5: 'Halogen'}
n_species = len(set(elements.values()))
# default AEV parameters
radial_cutoff = 10.4
angular_cutoff = 7.0
EtaR = 16.0
EtaA = 8.0
radial_div = 32
angular_div = 8
zeta = 16.0 #32.0
angular_sec = 8
# default AEV computer Args
AEV_ARGS = (
    radial_cutoff, angular_cutoff, 
    EtaR, EtaA, 
    radial_div, angular_div, 
    zeta, angular_sec, 
    n_species
)