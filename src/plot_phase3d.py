#!/usr/local/bin/python3

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Path3DCollection, Line3DCollection
import matplotlib.pyplot as plt
import numpy as np
from polyhedron import Vrep, Hrep
np.set_printoptions(precision=4)
'''
variables in cddlib 
 *   Inc   defines generator-halfspaces relations
 *   Adj   defines generator-generators relations
 *   InInc defines halfspace-generators relations
 *   InAdj defines halfspace-halfspaces relations
'''

def get_mu_dep(A0, b0, mu, indx_dep):
    a_dep = A0[indx_dep]
    A0 = np.delete(A0, indx_dep)
    mu_dep = (b0 - np.dot(A0, mu)) / a_dep
    return mu_dep[0]


def read_cplap_input(file='input.dat'):
    '''
    read cplap input file 
    and return A, b, and lists of secondary phase and elements
    '''
    A = []
    b = []
    # elements
    sec_phases = []

    with open(file) as f:
        lines = [line for line in f.readlines() if not line[0] == '#']
        n_elements = int(lines.pop(0))
        
        # read product
        line = lines.pop(0).split()
        b_i = float(line.pop())
        elements = line[1::2]
        A_i = np.array(line[0::2], dtype=int)
        dep_elem = lines.pop(0).strip()
        indx_dep = elements.index(dep_elem)
        dep_coeff = float(A_i[indx_dep])

        # product
        A.append(-A_i)
        b.append(-b_i)
        sec_phases.append(''.join(['{}{}'.format(i,j) for i,j in zip(elements, A_i)]))

        # elemental each
        for i in np.delete(range(len(elements)), [indx_dep]):
            A_ie = np.zeros(len(elements))
            A_ie[i] = A_i[i]
            A.append(A_ie)
            b.append(0)
            sec_phases.append(''.join(['{}{}'.format(i,j) for i,j in zip(elements, A_ie)]))
            A.append(-A_ie)
            b.append(-b_i)
            sec_phases.append(''.join(['{}{}'.format(i,j) for i,j in zip(elements, A_ie)]))

        # secondary phases
        no_sec = int(lines.pop(0))
        for i_sec in range(no_sec):
            n_elem = int(lines.pop(0)) # ignore
            
            line = lines.pop(0).split()
            b_i = float(line.pop())
            
            elem_i = line[1::2]
            indx_elem = [elements.index(elem) for elem in elem_i]
            A_it = np.array(line[0::2], dtype=int)
            A_i = np.array([A_it[indx_elem.index(i)] 
                if i in indx_elem else 0 
                for i in range(len(elements))])

            A.append(A_i+A[0]/dep_coeff*A_i[indx_dep])
            b.append(b_i+b[0]/dep_coeff*A_i[indx_dep])
            sec_phases.append(''.join(['{}{}'.format(i,j) for i,j in zip(elem_i, A_it)]))

    # remove dependent column
    A0 = A[0] # to calculate mu_dependent 
    A = np.delete(A, indx_dep, axis=1)
    b = np.array(b).reshape((len(b),1))
    elements.pop(indx_dep)

    mu_dep = lambda mu: get_mu_dep(A0, b[0], mu, indx_dep)

    return A, b, elements, sec_phases, mu_dep, dep_elem

def sort_vert(ininc_i, adj):
    ininc_sorted = []
    ininc_i = list(ininc_i)
    while len(ininc_i) > 0:
        v = ininc_i.pop()
        ininc_sorted.append(v)
        # find adj
        adj_i = adj[v]
        ininc_i = sorted(ininc_i, reverse=True, 
            key=lambda x: np.where(np.concatenate([adj_i, np.arange(1000)]) == x)[0][0])
    return ininc_sorted

def draw_plane(ax, verts, ininc, adj, sec_phases=[]):
    # find number of polygons with verts 
    if len([half for half in ininc if len(half) > 3]) > 8:
        cmap = plt.get_cmap("tab10_r")
    else:
        cmap = plt.get_cmap("Set1")

    # draw plane
    color_counter = 0
    for i, ininc_i in enumerate(ininc):
        if len(ininc_i) < 3:
            continue

        ininc_i = sort_vert(ininc_i, adj)
        x = []
        y = []
        z = []
        for v in ininc_i:
            x.append(verts[v][0])
            y.append(verts[v][1])
            z.append(verts[v][2])
        x.append(verts[ininc_i[0]][0])
        y.append(verts[ininc_i[0]][1])
        z.append(verts[ininc_i[0]][2])
        coord = [list(zip(x, y, z))]

        label = sec_phases[i]
        polygon = Poly3DCollection(coord, alpha=0.9, label=label, closed=True)
        polygon.set_facecolor(cmap(color_counter))
        color_counter += 1
        polygon._facecolors2d=polygon._facecolors3d
        polygon._edgecolors2d=polygon._edgecolors3d

        ax.add_collection3d(polygon)
        path = Line3DCollection(coord, lw=2, color='k') 
        ax.add_collection3d(path)
    ax.legend(loc='center left')
    ax.legend(loc='upper center', ncol=5)

def write_vert_output(generators, inc, elements, sec_phases, mu_dep, dep_elem, file='output_vert.dat'):
    with open(file, 'w') as f:
        [f.write('# {0: <7}'.format(elem)) for elem in elements + [dep_elem] + [': secondary phases']]
        f.write('\n')
        for i, vert in enumerate(generators):
            # write chemical potential
            [f.write('{0:0.4f}  '.format(v)) for v in np.hstack((vert, [mu_dep(vert)]))]
            f.write('  :')
            # write secondary phase
            [f.write(' {0: <15}'.format(sec_phases[inc_i])) for inc_i in inc[i]]
            f.write('\n')
        f.write('\n')

def write_half_output(generators, ininc, elements, sec_phases, mu_dep, dep_elem, file='output_half.dat'):
    with open(file, 'w') as f:
        f.write('# {0: <13}:'.format('Sec phase'))
        f.write(' Vertices \n')
        for i, half in enumerate(ininc):
            if len(half) < 3:
                continue
            f.write('{0: <15}:'.format(sec_phases[i]))
            for v_i in half:
                [f.write('  {0:0.4f}'.format(v)) for v in np.hstack((generators[v_i], [mu_dep(generators[v_i])]))]
                f.write(', ')
            f.write('\n')

def draw_pd(ax):
    """ tool_tip_missing
    """
    A, b, elements, sec_phases, mu_dep, dep_elem = read_cplap_input()

    p = Hrep(A, b)

    draw_plane(ax, p.generators, p.ininc, p.adj, sec_phases)
    
    write_vert_output(p.generators, p.inc, elements, sec_phases, mu_dep, dep_elem)
    write_half_output(p.generators, p.ininc, elements, sec_phases, mu_dep, dep_elem)

    set_axis(ax, A, b, elements)

def set_axis(ax, A, b, elements):
    # calc elemental chemical potential and set them 0
    buffer = 0.2 # eV

    mu_elements = []
    for i, ele in enumerate(elements):
        mu_i = []
        for k, a in enumerate(A):
            if b[k] == 0.: continue
            l_elemental = True
            for j, coeff in enumerate(a):
                if i != j and coeff != 0:
                    l_elemental = False
            if l_elemental:
                mu_i.append(b[k] / A[k, i])
        mu_elements.append(np.max(mu_i))
    print('mu_elements', elements, np.ravel(mu_elements))

    ax.set_xlim([-1 -buffer + mu_elements[0], 0 + buffer + mu_elements[0]])
    ax.set_ylim([-2 -buffer + mu_elements[1], -1 + buffer + mu_elements[1]])
    ax.set_zlim([-2 -buffer + mu_elements[2], 0 + buffer + mu_elements[2]])

    ax.set_xticks(-np.arange(2) + mu_elements[0])
    ax.set_yticks(-np.arange(2) + mu_elements[1] -1)
    ax.set_zticks(-np.arange(3) + mu_elements[2])

    ax.set_xticklabels(-1*np.arange(2))  
    ax.set_yticklabels(-1*np.arange(2)-1)
    ax.set_zticklabels(-1*np.arange(3))  

    ax.set_xlabel(r'$\mu_{{{}}}$ $(eV)$'.format(elements[0]))  
    ax.set_ylabel(r'$\mu_{{{}}}$ $(eV)$'.format(elements[1]))  
    ax.set_zlabel(r'$\mu_{{{}}}$ $(eV)$'.format(elements[2]))  

def main(ax):
    draw_pd(ax)

if __name__ == '__main__':
    fig = plt.figure()
    ax=Axes3D(fig)
    ax.set_aspect('equal')

    main(ax)
    ax.view_init(elev=20, azim=165)
    fig.savefig('pd.pdf')
    fig.savefig('pd.png', dpi=700)
    # plt.show()
