import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor
from Floyd_Warshall import FloydWarshall
from Bellman_Ford import BellmanFord
from dijkstra import DijkstraHeap
from paralleldijkstra import ParallelDijkstraAdjList
from SPFA import SPFA
class PPINetwork:
    def __init__(self, file_path=None, min_score=400, weight_type="inverse"):
        """
        初始化PPI网络对象
        
        :param file_path: STRING数据库文件路径
        :param min_score: 最小置信度阈值
        :param weight_type: 权重类型 ("inverse" 或 "direct")
        """
        self.graph = nx.Graph()
        self.df = None
        self.node_mapping = {}
        self.adjacency_list = None
        if file_path:
            self.load_data(file_path, min_score, weight_type)
    
    def load_data(self, file_path, min_score=400, weight_type="inverse"):
        """
        加载并预处理蛋白质互作数据
        
        :param file_path: 数据文件路径
        :param min_score: 最小置信度阈值
        :param weight_type: 权重类型 ("inverse" 或 "direct")
        """
        # 读取原始数据
        df = pd.read_csv(file_path, sep=' ')
        
        # 过滤低置信度互作
        df = df[df['combined_score'] >= min_score]
        
        # 提取关键列并重命名
        df = df[['protein1', 'protein2', 'combined_score']]
        df.columns = ['protein1', 'protein2', 'combined_score']
        
        # 去除物种前缀 (如"10090.")
        df['protein1'] = df['protein1'].str.split('.').str[-1]
        df['protein2'] = df['protein2'].str.split('.').str[-1]
        
        # 创建节点映射字典
        all_nodes = pd.concat([df['protein1'], df['protein2']]).unique()
        self.node_mapping = {node: idx for idx, node in enumerate(all_nodes)}
        
        # 添加反向映射
        self.id_to_node = {idx: node for node, idx in self.node_mapping.items()}
        
        # 保存处理后的数据
        self.df = df
        self.weight_type = weight_type
        
        # 创建图
        self.create_graph()
    
    def create_graph(self):
        """创建加权无向图"""
        if self.df is None:
            raise ValueError("请先加载数据")
        
        # 重置图
        self.graph = nx.Graph()
        
        # 添加所有节点
        for node in self.node_mapping.keys():
            self.graph.add_node(node)
        
        # 添加带权重的边
        for _, row in self.df.iterrows():
            protein1 = row['protein1']
            protein2 = row['protein2']
            
            # 计算权重
            if self.weight_type == "inverse":
                # 高置信度 → 低权重 (便于最短路径计算)
                weight = 1 - (row['combined_score'] / 1000.0)
            else:
                # 直接使用combined_score (0-1000)
                weight = row['combined_score'] / 1000.0
            
            # 添加边
            self.graph.add_edge(protein1, protein2, weight=weight)
    
    def get_graph(self):
        """获取NetworkX图对象"""
        return self.graph
    
    def get_node_mapping(self):
        """获取节点映射字典"""
        return self.node_mapping
    
    def get_id_to_node(self):
        """获取ID到节点的映射"""
        return self.id_to_node
    
    def create_adjacency_list(self):
        """创建邻接表表示"""
        if self.graph.number_of_nodes() == 0:
            raise ValueError("图尚未创建")
        
        self.adjacency_list = {}
        
        # 遍历所有节点
        for node in self.graph.nodes():
            self.adjacency_list[node] = []
            
            # 获取所有邻居
            for neighbor in self.graph.neighbors(node):
                weight = self.graph[node][neighbor]['weight']
                self.adjacency_list[node].append((neighbor, weight))
    
    def get_graph(self):
        """获取NetworkX图对象"""
        return self.graph
    
    def get_adjacency_list(self):
        """获取邻接表表示"""
        if self.adjacency_list is None:
            self.create_adjacency_list()
        return self.adjacency_list
    
    def run_bellman_ford(self, source):
        """
        运行Bellman-Ford算法计算从源点到所有点的最短路径
        :param source: 源点
        :return: BellmanFord对象
        """
        if source not in self.graph:
            raise ValueError(f"源点 {source} 不在图中")
        
        adjacency_list = self.get_adjacency_list()
        bf = BellmanFord(adjacency_list)
        success = bf.shortest_paths_from(source)
        
        if not success:
            print("警告: 检测到负权环，结果可能不准确")
        
        return bf
    
    def run_floyd_warshall(self):
        """
        运行Floyd-Warshall算法计算所有点对的最短路径
        :return: FloydWarshall对象
        """
        adjacency_list = self.get_adjacency_list()
        fw = FloydWarshall(adjacency_list)
        fw.compute()
        return fw
    
    def run_spfa(self, source):
        """
        运行SPFA算法计算从源点到所有点的最短路径
        :param source: 源点
        :return: SPFA对象
        """
        if source not in self.graph:
            raise ValueError(f"源点 {source} 不在图中")
        
        adjacency_list = self.get_adjacency_list()
        spfa = SPFA(adjacency_list)
        success = spfa.shortest_paths_from(source)
        
        if not success:
            print("警告: 检测到负权环，结果可能不准确")
        
        return spfa
    
    def run_dijkstra(self, source):
        """
        运行Dijkstra算法计算从源点到所有点的最短路径
        :param source: 源点
        :return: DijkstraHeap对象
        """
        if source not in self.graph:
            raise ValueError(f"源点 {source} 不在图中")
        adjacency_list = self.get_adjacency_list()
        dijkstra = DijkstraHeap(adjacency_list)
        dijkstra.shortest_paths_from(source)
        return dijkstra
    
    def run_parallel_dijkstra_adjlist(self, num_processes=None):
        """
        基于邻接表的并行Dijkstra全源最短路径
        """
        adjacency_list = self.get_adjacency_list()
        pd = ParallelDijkstraAdjList(adjacency_list)
        all_distances, node_list = pd.all_pairs_shortest_paths(num_processes=num_processes)
        return all_distances, node_list
    def analyze_graph(self):
        """分析图的基本属性"""
        if self.graph.number_of_nodes() == 0:
            raise ValueError("图尚未创建")
        
        print(f"图的基本属性:")
        print(f"- 节点数: {self.graph.number_of_nodes()}")
        print(f"- 边数: {self.graph.number_of_edges()}")
        
        # 计算平均度
        degrees = [d for n, d in self.graph.degree()]
        avg_degree = sum(degrees) / len(degrees) if degrees else 0
        print(f"- 平均度: {avg_degree:.2f}")
        
        # 连通分量分析
        components = list(nx.connected_components(self.graph))
        print(f"- 连通分量数: {len(components)}")
        if components:
            largest_cc = max(components, key=len)
            print(f"- 最大连通分量大小: {len(largest_cc)} 节点")
        
        # 权重统计
        if self.graph.number_of_edges() > 0:
            weights = [d['weight'] for _, _, d in self.graph.edges(data=True)]
            print(f"- 权重范围: {min(weights):.4f} - {max(weights):.4f}")
            print(f"- 平均权重: {sum(weights) / len(weights):.4f}")
        else:
            print("- 图中没有边")
        
        return {
            "node_count": self.graph.number_of_nodes(),
            "edge_count": self.graph.number_of_edges(),
            "avg_degree": avg_degree,
            "components": len(components),
            "largest_component": len(largest_cc) if components else 0
        }
    
    def visualize(self, max_nodes=100, output_file=None):
        """
        可视化图结构
        
        :param max_nodes: 最大可视化节点数
        :param output_file: 输出文件路径 (可选)
        """
        if self.graph.number_of_nodes() == 0:
            raise ValueError("图尚未创建")
        
        # 如果图过大，使用子图
        if self.graph.number_of_nodes() > max_nodes:
            # 提取最大连通分量中的前max_nodes个节点
            largest_cc = max(nx.connected_components(self.graph), key=len)
            subgraph_nodes = list(largest_cc)[:max_nodes]
            subgraph = self.graph.subgraph(subgraph_nodes)
            print(f"图过大，仅可视化前{max_nodes}个节点")
        else:
            subgraph = self.graph
        
        # 设置可视化参数
        plt.figure(figsize=(12, 10))
        pos = nx.spring_layout(subgraph, seed=42)
        
        # 绘制节点和边
        nx.draw_networkx_nodes(subgraph, pos, node_size=50, node_color='skyblue')
        
        # 绘制边 - 根据权重设置透明度
        edges = subgraph.edges(data=True)
        edge_alphas = [d['weight'] for _, _, d in edges]
        nx.draw_networkx_edges(subgraph, pos, alpha=0.1)
        
        # 添加标签（仅对度高的节点）
        high_degree_nodes = [n for n, d in subgraph.degree() if d > 5]
        labels = {n: n for n in high_degree_nodes}
        nx.draw_networkx_labels(subgraph, pos, labels, font_size=8)
        
        plt.title('Protein-Protein Interaction Network')
        plt.axis('off')
        plt.tight_layout()
        
        # 保存或显示图像
        if output_file:
            plt.savefig(output_file, dpi=300)
            print(f"图像已保存至 {output_file}")
        else:
            plt.show()
        
        plt.close()
    
    def save_graph(self, output_file):
        """保存图数据为GraphML格式"""
        if self.graph.number_of_nodes() == 0:
            raise ValueError("图尚未创建")
        
        nx.write_graphml(self.graph, output_file)
        print(f"图已保存为 {output_file}")

    def calculate_betweenness_centrality(self, k=None):
        """
        计算介数中心性
        
        :param k: 用于近似的节点数量 (None表示精确计算)
        :return: 介数中心性字典
        """
        if self.graph.number_of_nodes() == 0:
            raise ValueError("图尚未创建")
        
        print(f"计算介数中心性 (k={k if k else '精确'})...")
        
        if k and k < self.graph.number_of_nodes():
            return nx.betweenness_centrality(self.graph, k=k, weight='weight')
        else:
            return nx.betweenness_centrality(self.graph, weight='weight')

    def export_for_visualization(self, output_file, source=None, target=None, 
                            highlight_path=True, highlight_central=True, 
                            centrality_threshold=0.8):
        """
        导出用于 Cytoscape 或 Gephi 可视化的网络数据
        :param output_file: 输出文件路径
        :param source: 最短路径的源点
        :param target: 最短路径的目标点
        :param highlight_path: 是否高亮最短路径
        :param highlight_central: 是否高亮中心节点
        :param centrality_threshold: 中心节点阈值
        """
    # 创建图的副本
        G = self.graph.copy()
    
        # 计算介数中心性
        if highlight_central:
            betweenness = nx.betweenness_centrality(G, weight='weight')
            # 设置中心节点属性
            max_bc = max(betweenness.values())
            for node, bc in betweenness.items():
                G.nodes[node]['betweenness'] = bc
                # 标记高中心性节点
                G.nodes[node]['central'] = 1 if bc > centrality_threshold * max_bc else 0
        
        # 标记最短路径上的边
        if highlight_path and source and target:
            try:
                # 计算最短路径
                path = nx.shortest_path(G, source, target, weight='weight')
                # 标记路径上的节点
                for node in path:
                    G.nodes[node]['on_path'] = 1
                
                # 标记路径上的边
                path_edges = list(zip(path[:-1], path[1:]))
                for u, v in G.edges:
                    if (u, v) in path_edges or (v, u) in path_edges:
                        G[u][v]['on_path'] = 1
                    else:
                        G[u][v]['on_path'] = 0
            except nx.NetworkXNoPath:
                print(f"警告: {source} 和 {target} 之间没有路径")
        
        # 保存为 GraphML 格式
        nx.write_graphml(G, output_file)
        print(f"可视化数据已导出至 {output_file}")

    
# 使用示例
if __name__ == "__main__":
    # 1. 初始化并加载数据
    ppi_network = PPINetwork(
        file_path="7209.protein.physical.links.v12.0.txt",
        min_score=400,
        weight_type="inverse"  # 用于最短路径计算
    )
    
    # 2. 分析图结构
   # stats = ppi_network.analyze_graph()
    source_node = "A0A1I7V567"
    target_node = "A0A1I7W541"
    print("\n运行 Bellman-Ford 算法:")
    bf = ppi_network.run_bellman_ford(source_node)
    path = bf.get_shortest_path(target_node)
    distance = bf.distances[target_node]
    print(f"从 {source_node} 到 {target_node} 的最短路径: {' → '.join(path)} (距离: {distance:.4f})")

    print("\n运行 Dijkstra 算法:")
    dijkstra = ppi_network.run_dijkstra(source_node)
    dijkstra_path = dijkstra.get_shortest_path(target_node)
    dijkstra_distance = dijkstra.distances[target_node]
    print(f"从 {source_node} 到 {target_node} 的最短路径: {' → '.join(dijkstra_path)} (距离: {dijkstra_distance:.4f})")

    print("\n运行 Parallel Dijkstra 算法:")
    all_distances, node_list = ppi_network.run_parallel_dijkstra_adjlist(num_processes=4)

    try:
        source_idx = node_list.index(source_node)
        target_idx = node_list.index(target_node)
        distance = all_distances[source_idx][target_idx]
        # 用networkx查找路径（只查一对，速度快）
        path = nx.shortest_path(ppi_network.get_graph(), source=source_node, target=target_node, weight='weight')
        print(f"从 {source_node} 到 {target_node} 的最短路径: {' → '.join(path)} (距离: {distance:.4f})")
    except ValueError:
        print("源点或目标点不在节点列表中")
    except nx.NetworkXNoPath:
        print(f"{source_node} 到 {target_node} 不连通，无路径。")
'''
    ppi_network.export_for_visualization(
        output_file="ppi_visualization.graphml",
        source=source_node,
        target=target_node,
        highlight_path=True,
        highlight_central=True,
        centrality_threshold=0.85  # 仅显示前15%的中心节点
    )


    print("\n运行 Floyd-Warshall 算法:")
    fw = ppi_network.run_floyd_warshall()
    path = fw.get_path(source_node, target_node)
    distance = fw.get_distance(source_node, target_node)
    print(f"从 {source_node} 到 {target_node} 的最短路径: {' → '.join(path)} (距离: {distance:.4f})")
'''