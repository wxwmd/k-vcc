import networkx as nx
import operations
from networkx.algorithms.flow import maximum_flow
from networkx.algorithms.flow import minimum_cut
import matplotlib.pyplot as plt
from networkx.algorithms.connectivity import cuts

g=nx.read_edgelist('./graphs/test2.txt')

vcc=operations.kvcc_enum(g,3)