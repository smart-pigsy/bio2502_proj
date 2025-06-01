import heapq
import numpy as np
from concurrent.futures import ProcessPoolExecutor

def dijkstra_single_source_adjlist(args):
    adjacency_list, node_list, source_idx = args
    num_nodes = len(node_list)
    source = node_list[source_idx]
    distances = {node: float('inf') for node in node_list}
    distances[source] = 0
    heap = [(0, source)]
    while heap:
        current_dist, u = heapq.heappop(heap)
        if current_dist > distances[u]:
            continue
        for v, weight in adjacency_list.get(u, []):
            alt = current_dist + weight
            if alt < distances[v]:
                distances[v] = alt
                heapq.heappush(heap, (alt, v))
    # 返回与 node_list 顺序一致的距离数组
    return [distances[node] for node in node_list]

class ParallelDijkstraAdjList:
    def __init__(self, adjacency_list):
        """
        :param adjacency_list: dict，邻接表，{node: [(neighbor, weight), ...]}
        """
        self.adjacency_list = adjacency_list
        self.node_list = list(adjacency_list.keys())
        self.num_nodes = len(self.node_list)

    def all_pairs_shortest_paths(self, num_processes=None):
        all_distances = np.full((self.num_nodes, self.num_nodes), float('inf'))
        args = [
            (self.adjacency_list, self.node_list, idx)
            for idx in range(self.num_nodes)
        ]
        with ProcessPoolExecutor(max_workers=num_processes) as executor:
            results = executor.map(dijkstra_single_source_adjlist, args)
            for idx, distances in enumerate(results):
                all_distances[idx] = distances
        return all_distances, self.node_list
