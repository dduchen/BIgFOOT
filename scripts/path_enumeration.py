# THIS SCRIPT CONTAINS OLD PATH ENUMERATION FUNCTIONS; CURRENTLY NOT IN USE

def check_branches(graph, ncontigs, nvert):
    # keep track of start/end nodes per contig
    contigs_started = [0 for i in range(ncontigs)]
    starts_per_node = {v : [] for v in range(nvert)}
    contig_ends = {}
    for v in graph.vertices():
        for contig in list(graph.vp.contigs[v]):
            contig_ends[contig] = int(v)
            if contigs_started[int(contig)] == 0:
                starts_per_node[int(v)].append(contig)
                contigs_started[int(contig)] = 1
    ends_per_node = {v : [] for v in range(nvert)}
    for c, v in contig_ends.items():
        ends_per_node[v].append(c)
    # process all internal 2-sided branch nodes
    discarded_overlaps = []
    for v in graph.vertices():
        if v.out_degree() <= 1 or v.in_degree() <= 1:
            # not a 2-sided branch
            continue
        contigs_v = list(graph.vp.contigs[v])
        start_end_contigs = set(starts_per_node[int(v)] + ends_per_node[int(v)])
        bridged = (len(contigs_v) != len(start_end_contigs))
        if bridged:
            # discard all overlaps between a start and a stop contig
            for c1 in starts_per_node[int(v)]:
                for c2 in ends_per_node[int(v)]:
                    if c1 != c2:
                        discarded_overlaps.append((int(c1), int(c2)))
    return discarded_overlaps
    #         filtered_contigs = list(set(contigs_v) - start_end_contigs)
    #         graph.vp.contigs[v] = filtered_contigs
    #         contigs_adapted = contigs_adapted.union(start_end_contigs)
    # adapted_count = len(contigs_adapted)
    # print("{} contigs got adapted start/end points".format(adapted_count))
    # return graph, adapted_count


def enumerate_paths(g, abundances, low_ab_nodes, min_ab, reduce_paths, discarded):
    # load graph (.gt) from file
    #g = load_graph(graphfile)
    nvert = len(list(g.vertices()))
    contig_IDs = set()
    for contigs in g.vp.contigs:
        contig_IDs = contig_IDs.union(set(contigs))
    ncontigs = max([int(x) for x in contig_IDs])+1
    # B is matrix with dimension nodes x contigs, where B[i,j] is 1 if contig j
    # is passing through node i
    B = zeros((nvert,ncontigs))
    i = 0
    for v in g.vertices():
        for k in range(ncontigs):
            if str(k) in list(g.vp.contigs[v]):
                B[i,k]=1
        i += 1
    print("{} low abundance nodes (min_abundance = {})".format(
            len(low_ab_nodes), min_ab))
    # for v in low_ab_nodes:
    #     print("node {}, seq = {}".format(v, g.vp.seq[g.vertex(v)]))
    print()

    # now build the paths using BFS
    paths = BFS_paths(g, B, low_ab_nodes, nvert, ncontigs, discarded)

    if reduce_paths and len(low_ab_nodes) > 0:
        # reduce paths by replacing low-abundance nodes with most likely alternative
        paths = correct_paths(g, B, low_ab_nodes, abundances, paths, nvert,
                    ncontigs, min_ab)

    # store strains and paths in files
    seq_lengths = []
    f_haps = open("haps.fasta", 'w')
    f_paths = open("paths.txt", 'w')
    haps = []
    for i in range(len(paths)):
        p = paths[i]
        assert len(p) > 0
        f_paths.write(">path{}\n".format(i))
        hap = ""
        for node in p:
            if hap == "":
                f_paths.write("{}".format(node))
            else:
                f_paths.write(",{}".format(node))
            hap += g.vp.seq[g.vertex(node)]
            # if 'N' in g.vp.seq[g.vertex(node)]:
            #     print(int(node), g.vp.seq[g.vertex(node)])
        f_paths.write("\n")
        # write haplotype sequence to haps file
        f_haps.write(">path{}\n".format(i))
        f_haps.write("{}\n".format(hap))
        seq_lengths.append(len(hap))
        haps.append(hap)
    f_haps.close()
    f_paths.close()
    return [nvert, paths, haps]



def correct_paths(g, B, low_ab_nodes, abundances, paths, nvert, ncontigs, min_ab):
    print("\nCorrecting low-abundance nodes in haplotype paths")
    start_per_node = {v : set() for v in range(nvert)}
    end_per_node = {v : set() for v in range(nvert)}
    for k in range(ncontigs):
        start = nvert
        end = -1
        for i in range(nvert):
            if B[i,k]==1:
                end = max(end, i)
                start = min(start, i)
        start_per_node[start].add(k)
        end_per_node[end].add(k)
    corrected_paths = []
    corrected_haps = []
    hapnum = 0
    for path in tqdm(paths):
        hap = []
        if len(set(low_ab_nodes).intersection(
                set([int(x) for x in path]))) == 0:
            skip = True
        elif len(path) == 1:
            skip = True
        else:
            skip = False
        # trace path through graph, build contig set & list low-abundance nodes
        contigs = []
        active_contigs = set(g.vp.contigs[path[0]])
        assert len(active_contigs) > 0
        idx_v_to_correct = []
        idx = 0
        for node in path:
            hap.append(g.vp.seq[node])
            if skip:
                continue
            end_contigs = end_per_node[node]
            if node in low_ab_nodes:
                # print(int(node), abundances[int(node)])
                idx_v_to_correct.append(idx)
                if idx > 0:
                    for v_alt in g.vertex(path[idx-1]).out_neighbors():
                        active_contigs = active_contigs.union(start_per_node[v_alt])
                else:
                    for v_alt in g.vertex(path[idx+1]).in_neighbors():
                        active_contigs = active_contigs.union(start_per_node[v_alt])
            else:
                start_contigs = start_per_node[node]
                active_contigs = active_contigs.intersection(set(g.vp.contigs[node]))
                active_contigs = active_contigs.union(start_contigs)
            for c in active_contigs:
                if int(c) in end_contigs:
                    contigs.append(c)
            idx += 1
        if not skip:
            assert len(contigs) > 0
        # print(">old_hap_{}".format(hapnum))
        # print(''.join(hap))
        # find best alternative for low-abundance nodes suggested by contigs
        for i in range(len(idx_v_to_correct)):
            idx = idx_v_to_correct[i]
            v = path[idx]
            if idx == 0:
                v_next = path[idx+1]
                candidate_nodes = list(g.vertex(v_next).in_neighbors())
            elif idx == len(path)-1 or ( i+1<len(idx_v_to_correct) and
                                         idx_v_to_correct[i+1] == idx+1 ):
                v_prev = path[idx-1]
                candidate_nodes = list(g.vertex(v_prev).out_neighbors())
            else:
                v_next = path[idx+1]
                v_prev = path[idx-1]
                candidate_nodes = list(
                                    set(g.vertex(v_prev).out_neighbors()) &
                                    set(g.vertex(v_next).in_neighbors())
                )
            max_ab = abundances[int(v)]
            best_node = v
            best_seq = g.vp.seq[v]
            for node in candidate_nodes:
                if abundances[int(node)] <= max_ab:
                    continue
                node_seq = g.vp.seq[node]
                node_ab = abundances[int(node)]
                if (node_seq == best_seq
                        and node_ab + abundances[int(v)] >= min_ab):
                    # equivalent node, allow replacement independent of contigs
                    max_ab = node_ab
                    best_node = node
                    best_seq = node_seq
                    break
                else:
                    for k in g.vp.contigs[node]:
                        if k in contigs:
                            max_ab = abundances[int(node)]
                            best_node = node
                            break
            # replace the node with its best alternative
            path[idx] = best_node
            hap[idx] = best_seq

        corrected_paths.append(path)
        corrected_haps.append(''.join(hap))
        # print(">new_hap_{}".format(hapnum))
        # print(''.join(hap))
        hapnum += 1
    # filter out any duplicate paths
    remove_path_idx = []
    for h1 in range(len(corrected_haps)):
        hap1 = corrected_haps[h1]
        for h2 in range(h1+1, len(corrected_haps)):
            hap2 = corrected_haps[h2]
            if hap1 in hap2:
                remove_path_idx.append(h1)
            elif hap2 in hap1:
                remove_path_idx.append(h2)
    for idx in sorted(list(set(remove_path_idx)), reverse=True):
        corrected_paths.pop(idx)
    print("Corrected #paths from {} to {}.".format(len(paths), len(corrected_paths)))
    return corrected_paths


def get_start_nodes(B, start_end_per_contig, contig_ends_per_node, nvert, ncontigs):
    #print("get_start_nodes")
    global_start_nodes = []
    for k in range(ncontigs):
        # find the contig's start node and end node
        [start_k, end_k] = start_end_per_contig[k]
        if start_k == nvert: # empty contig path: all its nodes were trimmed
            continue
        elif start_k in global_start_nodes:
            continue
        active_contigs = B[start_k] # bitvec of contigs going through start_k
        path = [B[i,k] for i in range(nvert)] # bitvec indicating nodes on path
        current_node = start_k
        is_start = True
        while current_node <= end_k and is_start:
            if path[current_node] == 0:
                # current node not on contig path
                current_node += 1
                continue
            elif sum(active_contigs) <= 1:
                # no other active contigs remaining
                break
            elif current_node > start_k:
                active_contigs = [active_contigs[i] and B[current_node][i]
                                                    for i in range(ncontigs)]
            # check if any of our active contigs ends here
            for i in range(ncontigs):
                if (i != k and active_contigs[i]
                            and contig_ends_per_node[current_node][i]):
                    is_start = False
                    break
            current_node += 1
        if is_start:
            #print("start node {} because of contig {}.".format(start_k, k))
            global_start_nodes.append(start_k)
    return global_start_nodes


def get_end_nodes(B, start_end_per_contig, contig_starts_per_node, nvert, ncontigs):
    #print("get_end_nodes")
    global_end_nodes = []
    for k in range(ncontigs):
        # find the contig's start node and end node
        [start_k, end_k] = start_end_per_contig[k]
        if end_k == -1: # empty contig path: all its nodes were trimmed
            continue
        elif end_k in global_end_nodes:
            continue
        active_contigs = B[end_k] # bitvec of contigs going through end_k
        path = [B[i,k] for i in range(nvert)] # bitvec indicating nodes on path
        current_node = end_k
        is_end = True
        while current_node >= start_k and is_end:
            if path[current_node] == 0:
                # current node not on contig path
                current_node -= 1
                continue
            elif sum(active_contigs) <= 1:
                # no other active contigs remaining
                break
            elif current_node < end_k:
                active_contigs = [active_contigs[i] and B[current_node][i]
                                                    for i in range(ncontigs)]
            # check if any of our active contigs ends here
            for i in range(ncontigs):
                if (i != k and active_contigs[i]
                            and contig_starts_per_node[current_node][i]):
                    is_end = False
                    break
            current_node -= 1 # walk the contig path backwards!
        if is_end:
            #print("end node {} because of contig {}.".format(end_k, k))
            global_end_nodes.append(end_k)
    return global_end_nodes


def BFS_paths(g, B, low_ab_nodes, nvert, ncontigs, discard_overlaps):
    #print("BFS_paths...")
    #print(B)
    path_start_nodes = []
    path_end_nodes = []
    final_paths = []

    # build adj_out
    adj_out = []
    for vertex in g.vertices():
        v = int(vertex)
        adj_out.append([int(w) for w in vertex.out_neighbors()])

    # store start and end nodes per contig
    contig_starts_per_node = [[0 for k in range(ncontigs)] for i in range(nvert)]
    contig_ends_per_node = [[0 for k in range(ncontigs)] for i in range(nvert)]
    start_end_per_contig = []
    for k in range(ncontigs):
        start = nvert
        end = -1
        for i in range(nvert):
            if B[i,k]==1:
                end = max(end, i)
                start = min(start, i)
        start_end_per_contig.append([start, end])
        if start == nvert or end == -1:
            # contig has no nodes in trimmed graph
            continue
        assert start >= 0 and start < nvert
        assert end >= start and end < nvert
        contig_starts_per_node[start][k] = 1
        contig_ends_per_node[end][k] = 1

    path_start_nodes = get_start_nodes(B, start_end_per_contig,
                                    contig_ends_per_node, nvert, ncontigs)
    path_end_nodes = get_end_nodes(B, start_end_per_contig,
                                    contig_starts_per_node, nvert, ncontigs)
    print("path start nodes:", path_start_nodes)
    print("path end nodes:", path_end_nodes)
    #print("contig ends per node:")

    paths_to_process = []
    for v in path_start_nodes:
        # keep track of the contigs corresponding to our current path
        active_contigs = [int(x) for x in contig_starts_per_node[v]]
        #print(v, active_contigs)
        current_path = [v]
        reserve_contigs = [0 for k in range(ncontigs)]
        paths_to_process.append([current_path, active_contigs, reserve_contigs])

    progress = 0
    while len(paths_to_process) > 0:
        [path, active_contigs, reserve_contigs] = paths_to_process.pop(0) # BFS
        v = path[-1]
        #if g.vertex(v).out_degree() == 0: # check if we are in an end-node
        if v in path_end_nodes and sum([contig_ends_per_node[v][i] and
                active_contigs[i] for i in range(ncontigs)]) == sum(active_contigs):
            #print("path found:", path)
            final_paths.append(path)
            continue
        elif sum(active_contigs) == 0: # paths ends here
            # print("path ends here")
            # print(path)
            # print(active_contigs)
            final_paths.append(path)
            continue
        # elif progress%1000==0: # show progress
        #     print("current stack size:", len(paths_to_process))
        #     print("path length:", len(path))
        #     print("end node:", path[-1])
        #     print(active_contigs)
        progress += 1
        paths_to_process += extend_path(v, active_contigs, reserve_contigs,
                        adj_out, B, low_ab_nodes, contig_starts_per_node,
                        contig_ends_per_node, path, discard_overlaps)
    print("Total number of paths: {}".format(len(final_paths)))
    print("Final progress:", progress)
    return final_paths


def extend_path(node, active_contigs, reserve_contigs, adj_out, B, low_ab_nodes,
        contig_starts_per_node, contig_ends_per_node, path, discard_overlaps):
    # print("extend_path...",)
    ncontigs = len(active_contigs)
    nvert = len(B)
    path_list = []
    for v in adj_out[node]:
        if v == node:
            print("cycle skipped during path enumeration")
            continue
        # keep track of contigs starting and ending in this node
        start_contigs = contig_starts_per_node[v]
        end_contigs = [contig_ends_per_node[v][k] if active_contigs[k] else 0
                                                    for k in range(ncontigs)]
        # update active and reserve contigs
        if v in low_ab_nodes:
            # low-abundance node so we don't require a match at this position
            new_active_contigs = active_contigs
            new_reserve_contigs = reserve_contigs
            for v_alt in adj_out[node]:
                alt_start_contigs = contig_starts_per_node[v_alt]
                for k in range(ncontigs):
                    new_reserve_contigs[k] += alt_start_contigs[k]
        else:
            # update active contigs matching our current path
            new_active_contigs = [1 if (active_contigs[k] and int(B[v][k]))
                                                else 0 for k in range(ncontigs)]
            tmp_reserve_contigs = [1 if (reserve_contigs[k] and int(B[v][k]))
                                                else 0 for k in range(ncontigs)]
            new_reserve_contigs = [tmp_reserve_contigs[k] + start_contigs[k]
                                                for k in range(ncontigs)]
        if sum(end_contigs) > 0:
            # print("Contig ends here, allowing a switch to reserve contigs!")
            new_active_contigs = [new_active_contigs[k] + new_reserve_contigs[k]
                                                    for k in range(ncontigs)]
            new_reserve_contigs = [0 for k in range(ncontigs)]
            # discard any active contigs from discarded overlaps
            for k in range(ncontigs):
                if end_contigs[k] == 1:
                    for (c1, c2) in discard_overlaps:
                        if c1 == k:
                            new_active_contigs[c2] = 0
        # if sum(selected_contigs) > 0 or sum(start_contigs) > 0:
        if sum(new_active_contigs) > 0:
            # extend path
            new_path = path + [v]
            # remove end contigs
            for k in range(ncontigs):
                if end_contigs[k] == 1:
                    new_active_contigs[k] = 0
            active_count = sum(new_active_contigs)
            if active_count > 0:
                for i in range(v+1, nvert):
                    support = sum([B[i,k] if new_active_contigs[k]==1 else 0
                                                        for k in range(ncontigs)])
                    assert support >= 0 and support <= active_count
                    if support == 0:
                        continue
                    elif i in low_ab_nodes: # low abundance node
                        break
                    elif support < active_count: # path branches
                        break
                    elif sum(contig_starts_per_node[i]) > 0: # start node
                        break
                    elif sum(contig_ends_per_node[i]) > 0: # end node
                        break
                    else:
                        new_path += [i]
                        new_reserve_contigs = [1 if (
                                    new_reserve_contigs[k] and int(B[i,k])
                                ) else 0 for k in range(ncontigs)]
            path_list.append([new_path, new_active_contigs, new_reserve_contigs])
    return path_list
