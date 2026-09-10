"""Generate geometry only for the two tracked axisymmetric performance meshes."""
from pathlib import Path
from generate_b10_cax8t_meshes import write_exodus, append_nset, enrich_quad4_mesh

ROOT = Path(__file__).resolve().parent

def write_mesh(size, fuel_nr, clad_nr, element="cax4t"):
    xy, blocks, boundaries = [], [], []
    offset = 0
    for name, r0, r1, height, nr in [('fuel', 0., .004120, .010, fuel_nr),
                                   ('clad', .004122, .004692, .010020, clad_nr)]:
        first = len(xy) + 1
        def node(i, j): return first + j * (nr + 1) + i
        xy += [(r0 + (r1-r0)*i/nr, height*j/64) for j in range(65) for i in range(nr+1)]
        cells = [[node(i,j),node(i+1,j),node(i+1,j+1),node(i,j+1)] for j in range(64) for i in range(nr)]
        blocks.append((name,cells))
        boundaries += [(name+'_bottom',list(range(offset+1,offset+nr+1)),[1]*nr),
                       (name+'_right',[offset+(j+1)*nr for j in range(64)],[2]*64),
                       (name+'_top',list(range(offset+63*nr+1,offset+64*nr+1)),[3]*nr),
                       (name+'_left',[offset+j*nr+1 for j in range(64)],[4]*64)]
        offset += len(cells)
    if element == "cax8t":
        xy, blocks = enrich_quad4_mesh(xy, blocks)
    flat = [e for _,cells in blocks for e in cells]
    sets = [('all',list(range(1,len(xy)+1)))]
    for name,es,ss in boundaries:
        nodes = sorted({flat[e-1][q] for e,s in zip(es,ss) for q in ((s-1,s%4,s+3) if element == "cax8t" else (s-1,s%4))})
        sets.append((name,nodes))
    stem='rz_performance_'+size+('_cax8t' if element == 'cax8t' else '')
    write_exodus(ROOT.parent/'meshes'/(stem+'.e'),stem,xy,blocks,boundaries,sets)
    lines=['** Shared geometry with '+stem+'.e','*Node']
    lines += ['%d, %.17g, %.17g'%(i,r,z) for i,(r,z) in enumerate(xy,1)]
    label=0
    for name,cells in blocks:
        lines.append('*Element, type='+element.upper()+', elset='+name.upper())
        for cell in cells:
            label+=1;lines.append(', '.join(map(str,[label]+cell)))
    for name,nodes in sets:append_nset(lines,name,nodes)
    for name,es,ss in boundaries:
        lines.append('*Surface, type=ELEMENT, name=S_'+name.upper())
        lines += ['%d, S%d'%(e,s) for e,s in zip(es,ss)]
    (ROOT/'rz_performance'/(stem+'_mesh.inc')).write_text('\n'.join(lines)+'\n')

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--element', choices=('cax4t', 'cax8t'), default='cax4t')
    args = parser.parse_args()
    write_mesh('medium',100,16,args.element)
    if args.element == 'cax4t':
        write_mesh('large',200,32)
