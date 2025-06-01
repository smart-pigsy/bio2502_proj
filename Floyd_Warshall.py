class FloydWarshall:
    def __init__(self, graph):
        """
        使用邻接表初始化图对象
        :param graph: 邻接表表示的图，格式为字典:
                     {节点: [(邻居1, 权重1), (邻居2, 权重2), ...]}
        """
        # 获取所有节点
        self.vertices = list(graph.keys())
        self.num_vertices = len(self.vertices)
        
        # 创建节点到索引的映射
        self.index_map = {v: i for i, v in enumerate(self.vertices)}
        self.reverse_map = {i: v for i, v in enumerate(self.vertices)}
        
        # 初始化距离矩阵和路径矩阵
        self.dist = [[float('inf')] * self.num_vertices for _ in range(self.num_vertices)]
        self.next = [[None] * self.num_vertices for _ in range(self.num_vertices)]
        
        # 初始化自身到自身的距离为 0
        for i in range(self.num_vertices):
            self.dist[i][i] = 0
            self.next[i][i] = i
        
        # 根据邻接表添加边
        for u, neighbors in graph.items():
            u_idx = self.index_map[u]
            for v, weight in neighbors:
                v_idx = self.index_map[v]
                # 如果存在多条边，保留权重最小的一条
                if weight < self.dist[u_idx][v_idx]:
                    self.dist[u_idx][v_idx] = weight
                    self.next[u_idx][v_idx] = v_idx

    def compute(self):
        """
        执行 Floyd-Warshall 算法计算所有点对的最短路径
        """
        for k in range(self.num_vertices):
            for i in range(self.num_vertices):
                # 如果 i->k 不可达，跳过
                if self.dist[i][k] == float('inf'):
                    continue
                for j in range(self.num_vertices):
                    # 检查是否存在更短路径
                    new_dist = self.dist[i][k] + self.dist[k][j]
                    if new_dist < self.dist[i][j]:
                        self.dist[i][j] = new_dist
                        self.next[i][j] = self.next[i][k]

    def get_path(self, u, v):
        """
        获取从 u 到 v 的最短路径列表
        :param u: 起点
        :param v: 终点
        :return: 路径列表，如 ['A', 'B', 'C']
        """
        u_idx, v_idx = self.index_map[u], self.index_map[v]
        
        # 检查路径是否存在
        if self.next[u_idx][v_idx] is None:
            return []
        
        # 重建路径
        path = []
        current = u_idx
        while current != v_idx:
            path.append(self.reverse_map[current])
            current = self.next[current][v_idx]
        path.append(self.reverse_map[v_idx])
        
        return path

    def get_distance(self, u, v):
        """
        获取从 u 到 v 的最短距离
        """
        u_idx, v_idx = self.index_map[u], self.index_map[v]
        return self.dist[u_idx][v_idx]

    def get_all_pairs_shortest_paths(self):
        """
        获取所有点对的最短距离和路径
        :return: 字典 { (起点, 终点): (距离, 路径) }
        """
        result = {}
        for u in self.vertices:
            for v in self.vertices:
                dist = self.get_distance(u, v)
                path = self.get_path(u, v) if dist != float('inf') else []
                result[(u, v)] = (dist, path)
        return result