#!/usr/local/bin/python3

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Path3DCollection, Line3DCollection
import matplotlib.pyplot as plt
import numpy as np
from polyhedron import Vrep, Hrep
'''
 *   Inc   defines generator-halfspaces relations
 *   Adj   defines generator-generators relations
 *   InInc defines halfspace-generators relations
 *   InAdj defines halfspace-halfspaces relations
'''

def read_dd_data_dd(file='results.txt'):
    '''
    double description method
    see ftp://ftp.ifor.math.ethz.ch/pub/fukuda/cdd/cddlibman/node3.html
    P = { x=(x1, ..., xd)^T :  b - A  x  >= 0 }
    '''
    """
    read 
    4.000 mu_S      1.000 mu_Zn     1.000 mu_Sn  >   -38.6731
    8.000 mu_S      0.000 mu_Zn     0.000 mu_Sn  <   -41.8860
    return 
    A, b, and set of 
    """
    with open(file) as f:
        elements = []
        lines = f.readlines()
        A = np.zeros([len(lines), 3])
        b = np.zeros([len(lines), 1])
        for j, line in enumerate(lines):
            for i in range(3):
                # get coefficient
                coeff_str = line[i*16:i*16+9].strip()
                A[j, i] = 0 if len(coeff_str) <1 else float(coeff_str)
                
                # get inequality
                ele_str = line[9+i*16:9+i*16+6].strip('mu_').strip()
                if len(ele_str) > 0:
                    elements.append(ele_str)
            i+=1
            inequal = line[i*16:i*16+3].strip()
            rhs = float(line[i*16+3:].strip())
            b[j] = rhs
            if inequal == '>':
                A[j, :] *= -1
                b[j] *= -1
            else:
                A[j, :] *= 1
                b[j] *= 1
    return A, b, sorted(set(elements), key=elements.index)

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

def draw_plane(ax, verts, ininc, adj):
    # find number of polygons with verts 
    ininc = [half for half in ininc if len(half) > 2]
    if len(ininc) > 8:
        cmap = plt.get_cmap("tab10_r")
    else:
        cmap = plt.get_cmap("Set2")

    # draw plane
    for i, ininc_i in enumerate(ininc):
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

        polygon = Poly3DCollection(coord, alpha=0.8, closed=True)
        polygon.set_facecolor(cmap(i))

        ax.add_collection3d(polygon)
        path = Line3DCollection(coord, lw=2, color='k') 
        ax.add_collection3d(path)

def draw_pd(ax):
    """ tool_tip_missing
    """
    A, b, elements = read_dd_data_dd()
    
    print(elements)
    hrep = np.append(b, A, axis=1)

    p = Hrep(A, b)
    print(p.generators)
    draw_plane(ax, p.generators, p.ininc, p.adj)

    set_axis(ax, A, b, elements)

def set_axis(ax, A, b, elements):
    # calc elemental chemical potential and set them 0
    buffer = 0.2 # eV
    mu_elements = []
    for i, ele in enumerate(elements):
        for k, a in enumerate(A):
            l_elemental = True
            for j, coeff in enumerate(a):
                if i != j and coeff != 0:
                    l_elemental = False
            if l_elemental:
                mu_elements.append(b[k] / A[k, i])
                break
    print('mu_elements', np.ravel(mu_elements))

    ax.set_xlim([-1 -buffer + mu_elements[0], 0 + buffer + mu_elements[0]])
    ax.set_ylim([-2 -buffer + mu_elements[1], -1 + buffer + mu_elements[1]])
    ax.set_zlim([-2 -buffer + mu_elements[2], 0 + buffer + mu_elements[2]])

    ax.set_xticks(-np.arange(2) + mu_elements[0])
    ax.set_yticks(-np.arange(2) + mu_elements[1] -1)
    ax.set_zticks(-np.arange(3) + mu_elements[2])

    ax.set_xticklabels(-1*np.arange(2))  
    ax.set_yticklabels(-1*np.arange(2)-1)
    ax.set_zticklabels(-1*np.arange(3))  
    # ax.set_aspect('equal', 'datalim')

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
    plt.show()
