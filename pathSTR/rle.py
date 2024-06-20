from itertools import groupby
import rustworkx
import logging


def setup_graph(s, k):
    """
    This code was written by Rob Patro
    """
    eps = 1.0 / k
    graph = rustworkx.PyDAG()

    node_list = []

    for o in range(k):
        nl = []
        for i in range(o, len(s) - k, k):
            nl.append(graph.add_node((i, i + k)))
        node_list.append(nl)

    char_nodes = []
    for i in range(len(s)):
        char_nodes.append(graph.add_node(i))

    # add edges between tokens in the same phase
    for nodes in node_list:
        for i in range(1, len(nodes)):
            dat_x = graph.get_node_data(nodes[i - 1])
            dat_y = graph.get_node_data(nodes[i])
            if s[dat_x[0] : dat_x[1]] == s[dat_y[0] : dat_y[1]]:
                graph.add_edge(nodes[i - 1], nodes[i], 0.1)
            else:
                graph.add_edge(nodes[i - 1], nodes[i], k + eps)

    for i in range(1, len(char_nodes)):
        graph.add_edge(char_nodes[i - 1], char_nodes[i], 1)
        m = i % k
        idx = i // k
        # the character can connect to this token
        if len(node_list[m]) < idx:
            graph.add_edge(char_nodes[i - 1], node_list[m][idx], k + eps)

    for o in range(k):
        for idx, node in enumerate(node_list[o]):
            # this token can connect back to the next character
            end_pos = o + (idx + 1) * k
            if end_pos + 1 < len(s):
                graph.add_edge(node_list[o][idx], char_nodes[end_pos + 1], 1)

    start_node = graph.add_node("start")
    sink_node = graph.add_node("end")

    graph.add_edge(start_node, char_nodes[0], 1)
    try:
        graph.add_edge(start_node, node_list[0][0], k + eps)
    except IndexError:
        logging.error(f"IndexError (1) when creating rle-graph for {s} with k={k}")
        return None, None

    graph.add_edge(char_nodes[-1], sink_node, 0)
    # which shift connects to the end
    m = len(s) % k
    try:
        graph.add_edge(node_list[m][-1], sink_node, 0)
    except IndexError:
        logging.error(f"IndexError (2) when creating rle-graph for {s} with k={k}")
        return None, None

    nnodes = sink_node + 1
    best_path_weight = [float("inf")] * nnodes
    best_path_weight[start_node] = 0

    traceback = [-1] * nnodes
    traceback[start_node] = start_node

    for n in rustworkx.topological_sort(graph):
        for pred in graph.predecessor_indices(n):
            weight = graph.get_edge_data(pred, n)
            if best_path_weight[pred] + weight < best_path_weight[n]:
                traceback[n] = pred
                best_path_weight[n] = min(
                    best_path_weight[pred] + weight, best_path_weight[n]
                )

    p = []
    x = sink_node
    while traceback[x] != x:
        x = traceback[x]
        p.append(x)

    return list(reversed(p)), graph


def make_rle(s, p, graph):
    """
    This code was written by Rob Patro
    """
    out = []
    for n in p[1:]:
        ndat = graph.get_node_data(n)
        if isinstance(ndat, tuple):
            out += [s[ndat[0] : ndat[1]]]
        elif isinstance(ndat, int):
            out += [s[ndat]]
        elif isinstance(ndat, str):
            out += [ndat]

    # Use 'groupby' to group consecutive elements and count their occurrences
    os = []
    for key, group in groupby(out):
        l = len(list(group))
        if len(key) > 1 and l > 1:
            os.append(f"({key}){to_subscript(l)}")
        elif l > 1:
            os.append(f"{key*l}")
        else:
            os.append(f"{key}")
    return "".join(os)


def rle(sequence, motif_length):
    """
    This function takes a sequence and a motif length and returns a compact run-length encoded sequence

    kudos to Rob Patro and https://stackoverflow.com/a/78634539/6631639
    """
    if not sequence:
        return None
    elif len(sequence) < motif_length * 2:
        return sequence
    else:
        p, graph = setup_graph(sequence, motif_length)
        if graph:
            return make_rle(sequence, p, graph)
        else:
            return None


def to_subscript(n):
    """
    Convert a number to subscript
    """
    return "".join(chr(0x2080 + int(i)) for i in str(n))
