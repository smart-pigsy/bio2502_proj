class BellmanFord:
    def __init__(self, graph):
        """
        初始化图对象
        :param graph: 邻接表表示的图，格式为字典:
                     {节点: [(邻居1, 权重1), (邻居2, 权重2), ...]}
        """
        self.graph = graph
        self.vertices = list(graph.keys())
        self.distances = {}
        self.predecessors = {}
        self.has_negative_cycle = False

    def shortest_paths_from(self, source):
        """
        计算从源点 source 到所有点的最短路径
        :param source: 源点
        :return: 如果无负环则返回 True，否则返回 False
        """
        # 初始化距离和前驱
        self.distances = {v: float('inf') for v in self.vertices}
        self.predecessors = {v: None for v in self.vertices}
        self.distances[source] = 0
        self.has_negative_cycle = False

        # 转换邻接表为边列表
        edges = []
        for u in self.graph:
            for v, w in self.graph[u]:
                edges.append((u, v, w))

        # 进行 V-1 次松弛操作
        for _ in range(len(self.vertices) - 1):
            updated = False
            for u, v, w in edges:
                if self.distances[u] + w < self.distances[v]:
                    self.distances[v] = self.distances[u] + w
                    self.predecessors[v] = u
                    updated = True
            if not updated:
                break  # 提前退出优化

        # 检查是否存在负权环
        for u, v, w in edges:
            if self.distances[u] + w < self.distances[v]:
                self.has_negative_cycle = True
                return False  # 存在负权环

        return True

    def get_shortest_path(self, target):
        """
        获取从源点到目标点的最短路径（调用 shortest_paths_from 后）
        :param target: 目标点
        :return: 路径列表，如果无法到达则为空列表
        """
        if self.distances[target] == float('inf'):
            return []

        path = []
        current = target
        while current is not None:
            path.append(current)
            current = self.predecessors.get(current)
        path.reverse()
        return path

    def get_all_shortest_paths(self):
        """
        获取所有节点的最短路径
        :return: 字典 {节点: (距离, 路径)}
        """
        result = {}
        for v in self.vertices:
            result[v] = (self.distances[v], self.get_shortest_path(v))
        return result