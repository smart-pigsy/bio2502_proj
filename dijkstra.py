import heapq

class DijkstraHeap:
    def __init__(self, graph):
        """
        初始化图对象
        :param graph: 字典表示的图，形如 {u: [(v1, w1), (v2, w2), ...]}
        """
        self.graph = graph
        self.distances = {}
        self.previous = {}

    def shortest_paths_from(self, source):
        """
        从源点 source 计算最短路径
        :param source: 源点编号
        :return: 返回距离字典
        """
        self.distances = {node: float('inf') for node in self.graph}
        self.previous = {node: None for node in self.graph}
        self.distances[source] = 0

        heap = [(0, source)]  # (distance, node)

        while heap:
            current_dist, u = heapq.heappop(heap)

            # 若该路径不是当前最短路径，则跳过（过期节点）
            if current_dist > self.distances[u]:
                continue

            for v, weight in self.graph.get(u, []):
                alt = current_dist + weight
                if alt < self.distances[v]:
                    self.distances[v] = alt
                    self.previous[v] = u
                    heapq.heappush(heap, (alt, v))

        return self.distances

    def get_shortest_path(self, target):
        """
        获取从源点到目标点的最短路径（调用 shortest_paths_from 之后）
        :param target: 目标节点
        :return: 节点路径列表
        """
        path = []
        while target is not None:
            path.append(target)
            target = self.previous[target]
        path.reverse()
        return path
