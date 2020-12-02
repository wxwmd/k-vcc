import copy

import networkx as nx
from networkx.algorithms.connectivity import cuts
from networkx.algorithms.traversal.breadth_first_search import bfs_edges
from networkx.algorithms.flow import maximum_flow_value
from networkx.classes.function import common_neighbors
from networkx.algorithms.shortest_paths.generic import shortest_path_length


def directed_flow_graph(g: nx.Graph) -> nx.DiGraph:
    """
    transform the original graph G into a directed flow graph G'.
    :param g:original graph G
    :return:a directed flow graph G'
    """
    flow_graph = nx.DiGraph()
    for node in g.nodes:
        flow_graph.add_edge(node, str(node) + '#', capacity=1)
    for edge in g.edges:
        flow_graph.add_edge(str(edge[0]) + '#', edge[1], capacity=1)
        flow_graph.add_edge(str(edge[1]) + '#', edge[0], capacity=1)
    return flow_graph


def vertex_cut(flow_g, s, t) -> set:
    """
    在流向图计算能分割s和t的最小点集
    :param flow_g: 流向图
    :param s: 源点
    :param t: 终点
    :return: 所有分割的点集
    """
    nodes = cuts.minimum_node_cut(flow_g, s, t)
    s = set()
    for node in nodes:
        node = str(node).split('#')[0]
        s.add(node)
    return s


def overlap_partition(g: nx.Graph, vertex_cut: set) -> list:
    """
    将g用vertex_cut分为一些子图，每个子图用set{nodes}来表示
    :param g: 要切分的图
    :param vertex_cut: 分割g的点集
    :return: 子图的集合 {{1，2}，{3，4，5}}
    """
    s = []
    g1 = copy.deepcopy(g)
    g1.remove_nodes_from(vertex_cut)
    components = list(nx.connected_components(g1))

    for component in components:
        component = component.union(vertex_cut)
        s.append(component)
    return s


def sparse_certificate(g: nx.Graph, k):
    """
    生成g的稀疏等价图，同时生成side-groups
    :param g: 原图
    :return: 稀疏等价图
    """
    g1 = copy.deepcopy(g)
    side_groups = []
    for i in range(k):
        components = list(nx.connected_components(g1))
        # 生成side-groups
        if i == k - 1:
            id = 0
            for component in components:
                if len(component) > 1:
                    side_group = SideGroup(id, 0, component)
                    side_groups.append(side_group)
                    for node in component:
                        g1.nodes[node]['cc'] = id
                    id += 1

        sources = [component.pop() for component in components]
        edges = set()
        for source in sources:
            edges = edges.union(set(bfs_edges(g1, source)))
        g1.remove_edges_from(edges)
    remove_edges = copy.deepcopy(g1.edges)
    g1.add_edges_from(g.edges)
    g1.remove_edges_from(remove_edges)
    return g1, side_groups


def local_cut(u, v, flow_g: nx.Graph, g: nx.Graph, k):
    """
    找到在u,v之间存在的最小点割，如果不存在返回None
    :param u: 起点，存在于g中
    :param v: 终点，存在于g中
    :param flow_g:流图
    :param g: 原图
    :return: vertex-cut
    """
    u1 = str(u) + '#'
    if g.has_edge(u, v) or u == v:
        return None
    if maximum_flow_value(flow_g, u1, v) >= k:
        return None
    else:
        return vertex_cut(flow_g, u1, v)


def global_cut(g: nx.Graph, k):
    """
    输入一个图g和一个k
    输出图g中包含的长度小于k的点切
    如果不存在返回None
    :param g: 图
    :param k: 点切规模上限
    :return: 点切
    """
    sc = sparse_certificate(g, k)[0]
    min_degree = min([x[1] for x in list(sc.degree)])
    u = 0
    for node in g.nodes:
        if sc.degree[node] == min_degree:
            u = node
    flow_sc = directed_flow_graph(sc)

    for v in sc.nodes:
        cut = local_cut(u, v, flow_sc, sc, k)
        if cut is not None:
            return cut

    neighbors = list(sc.adj[u].keys())
    for va in neighbors:
        for vb in neighbors:
            cut = local_cut(va, vb, flow_sc, sc, k)
            if cut is not None:
                return cut
    return None


def kvcc_enum(g: nx.Graph, k) -> list:
    """
    寻找图中所有的k-vcc
    :param g: 输入的无向图
    :param k: k
    :return: set，里面的每一个元素是set，每个set表示一个k-vcc的点集
    """
    vcc = []
    remove_nodes = [node for node in g.nodes if g.degree[node] < k]
    g.remove_nodes_from(remove_nodes)
    components = nx.connected_components(g)
    for component in components:
        component_g = g.subgraph(component).copy()
        cut = global_cut_optimization(component_g, k)
        if cut is None:
            vcc.append(component)
        else:
            sub_components = overlap_partition(component_g, cut)
            for sub_component in sub_components:
                sub_component_g = component_g.subgraph(sub_component).copy()
                vcc = vcc + kvcc_enum(sub_component_g, k)
    return vcc


class SideGroup:
    id = 0
    deposit = 0
    vertex = set()

    def __init__(self, id, deposit, vertex):
        self.id = id
        self.deposit = deposit
        self.vertex = set()
        self.vertex = self.vertex.union(vertex)


def strong_side_vertices(g: nx.Graph, k):
    """
    计算图中的所有strong_side_vertices
    :param g:
    :param k:
    :return:
    """
    side_vettices = set()
    for u in g.nodes:
        flag = True
        neighbors = list(g.adj[u].keys())
        for v in neighbors:
            for v1 in neighbors:
                if v == v1 or g.has_edge(v, v1) or len(list(common_neighbors(g, v, v1))) >= k:
                    continue
                else:
                    flag = False
                    break
            if not flag:
                break
        g.nodes[u]['side_vertex'] = flag
        if flag:
            side_vettices.add(u)
    return side_vettices


def global_cut_optimization(g: nx.Graph, k):
    """
    这个函数实现优化之后的global_cut
    :param g: 图
    :param k:
    :return: 点切
    """
    sc, cs = sparse_certificate(g, k)
    flow_sc = directed_flow_graph(sc)
    side_vertices = strong_side_vertices(sc, k)
    u = 0
    if len(side_vertices) == 0:
        min_degree = min([x[1] for x in list(sc.degree)])
        for node in g.nodes:
            if sc.degree[node] == min_degree:
                u = node
    else:
        u = side_vertices.pop()

    for v in sc.nodes:
        sc.nodes[v]['pru'] = False
        sc.nodes[v]['deposit']=0
    sweep(u, cs, sc, k)

    dist = [v for v in sc.nodes]
    dist.sort(key=lambda v: shortest_path_length(sc, u, v), reverse=True)

    for v in dist:
        if sc.nodes[v]['pru']:
            continue
        else:
            cut = local_cut(u, v, flow_sc, sc, k)
            if cut:
                return cut
            else:
                sweep(v, cs, sc, k)


def sweep(v, cs: list, g: nx.Graph, k):
    """
    跳过对v的check
    :param g:
    :param v: 节点
    :param cs: SideGroup的信息.
    :return:None
    """
    g.nodes[v]['pru'] = True
    for w in list(g.adj[v].keys()):
        if not g.nodes[w]['pru']:
            g.nodes[w]['deposit'] += 1
            if g.nodes[v]['side_vertex'] or g.nodes[w]['deposit'] >= k:
                sweep(w, cs, g, k)
    if not g.nodes[v].get('cc', -1) == -1:
        cc = g.nodes[v]['cc']
        if cs[cc].deposit < k:
            cs[cc].deposit += 1
            if g.nodes[v]['side_vertex'] or cs[cc].deposit >= k:
                cs[cc].deposit = k
                for w in cs[cc].vertex:
                    if not g.nodes[w]['pru']:
                        sweep(w, cs, g, k)
