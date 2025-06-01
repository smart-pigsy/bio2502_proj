from collections import deque

class SPFA:
    def __init__(self, graph):
        """
        使用邻接表初始化图对象
        :param graph: 邻接表表示的图，格式为字典:
                     {节点: [(邻居1, 权重1), (邻居2, 权重2), ...]}
        """
        self.graph = graph
        self.vertices = list(graph.keys())
        
        # 初始化数据结构
        self.dist = {v: float('inf') for v in self.vertices}
        self.in_queue = {v: False for v in self.vertices}
        self.count = {v: 0 for v in self.vertices}
        self.prev = {v: None for v in self.vertices}
        self.has_negative_cycle = False

    def shortest_paths_from(self, source):
        """
        计算从源点到所有点的最短路径
        :param source: 源点
        :return: 如果没有负权环则返回 True，否则返回 False
        """
        # 重置状态
        self.dist = {v: float('inf') for v in self.vertices}
        self.in_queue = {v: False for v in self.vertices}
        self.count = {v: 0 for v in self.vertices}
        self.prev = {v: None for v in self.vertices}
        self.has_negative_cycle = False
        
        # 设置源点
        self.dist[source] = 0
        queue = deque([source])
        self.in_queue[source] = True
        self.count[source] = 1

        while queue:
            u = queue.popleft()
            self.in_queue[u] = False

            # 遍历u的所有邻居
            for neighbor, weight in self.graph.get(u, []):
                # 松弛操作
                if self.dist[u] + weight < self.dist[neighbor]:
                    self.dist[neighbor] = self.dist[u] + weight
                    self.prev[neighbor] = u
                    
                    # 如果邻居不在队列中，加入队列
                    if not self.in_queue[neighbor]:
                        queue.append(neighbor)
                        self.in_queue[neighbor] = True
                        self.count[neighbor] += 1
                        
                        # 检查负环
                        if self.count[neighbor] > len(self.vertices):
                            self.has_negative_cycle = True
                            return False
        
        return True

    def get_shortest_path(self, target):
        """
        获取从源点到目标点的路径
        :param target: 目标点
        :return: 路径列表，如果无法到达则返回空列表
        """
        if self.dist[target] == float('inf'):
            return []

        path = []
        current = target
        while current is not None:
            path.append(current)
            current = self.prev.get(current)
        return path[::-1]

    def get_all_shortest_paths(self):
        """
        获取所有节点的最短路径
        :return: 字典 {节点: (距离, 路径)}
        """
        result = {}
        for v in self.vertices:
            result[v] = (self.dist[v], self.get_shortest_path(v))
        return result