#!/usr/bin/env python3
from ete3 import Tree
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-t', dest='intree', type=str, required=True, help="""Input Tree File""")
parser.add_argument('-l', dest='leaf', type=str, required=True, help="""Leaf of Interest""")
args = parser.parse_args()

def get_closest_leaf(tree, leafoi):
    mindist = 1000000
    closestleaf = None
    leaves = tree.get_leaves()
    for leaf in leaves:
        distance = tree.get_distance(leafoi, leaf)
        if distance < mindist and leaf != leafoi:
            mindist = distance
            closestleaf = leaf

    return closestleaf.name, round(mindist, 6)


def loadtree(treefile, leafname):
    tree = Tree(treefile)
    try:
        leaf = tree.get_leaves_by_name(leafname)[0]
    except IndexError:
        raise RuntimeError('Could not find leaf named {} in the tree file.'.format(leafname))
    return tree, leaf


tree, leaf = loadtree(args.intree, args.leaf)
print(get_closest_leaf(tree, leaf))
