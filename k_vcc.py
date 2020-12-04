import networkx as nx
import operations
from networkx.algorithms.flow import maximum_flow
from networkx.algorithms.flow import minimum_cut
import matplotlib.pyplot as plt
from networkx.algorithms.connectivity import cuts
import time
from networkx import k_core


def k_vcc(g: nx.Graph, k):
    operations.strong_side_vertices(g, k)
    return operations.kvcc_enum(g, k)


g = nx.read_edgelist('./graphs/test7.txt')
nx.draw(g, with_labels=True)
plt.show()

k = 12
print(time.localtime(time.time()))
vcc = k_vcc(g, k)
for x in vcc:
    subg = g.subgraph(x)
    nx.draw(subg, with_labels=True)
    plt.show()
print(time.localtime(time.time()))
