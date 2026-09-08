"""Geometry only: attach individual node sets to the existing PCMI contact mesh."""
from scipy.io import netcdf_file
from generate_b10_cax8t_meshes import ROOT, MESH_ROOT, enrich_quad4_mesh, write_exodus
from generate_b13_pcmi_meshes import names


def main():
    source = ROOT.parent / 'moose' / 'm23_pcmi_coupled_cladding_rz_mesh.e'
    with netcdf_file(str(source), 'r', mmap=False) as file:
        xy = list(zip(file.variables['coordx'][:].tolist(), file.variables['coordy'][:].tolist()))
        blocks = [(name, file.variables['connect%d' % i][:].tolist())
                  for i, name in enumerate(names(file.variables['eb_names']), 1)]
        sides = [(name, file.variables['elem_ss%d' % i][:].tolist(), file.variables['side_ss%d' % i][:].tolist())
                 for i, name in enumerate(names(file.variables['ss_names']), 1)]
    for element in ('cax4t', 'cax4rt', 'cax8t', 'cax8rt'):
        coordinates, topology = enrich_quad4_mesh(xy, blocks) if '8' in element else (xy, blocks)
        node_sets = [('n_%d' % n, [n]) for n in range(1, len(coordinates) + 1)]
        write_exodus(MESH_ROOT / ('b14_nts_' + element + '.e'), 'Native NTS contact identification',
                     coordinates, topology, sides, node_sets)


if __name__ == '__main__':
    main()
