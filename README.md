这个包实现了蛋白质互作网络（PPI Network）的构建与分析，并集成了多种经典最短路径算法（Bellman-Ford、Dijkstra、SPFA、Floyd-Warshall）。下面是详细的用户手册和一个测试程序示例。

---

# 用户手册

## 1. 安装方法和主要功能
**安装方法**
在dist的文件夹中，有一个文件：ppinetwork-0.1.0-py3-none-any.whl

可以通过以下方法安装
` pip install path/to/ppinetwork-0.1.0-py3-none-any.whl`
- **PPI网络构建**：从STRING数据库文件加载蛋白质互作数据，自动过滤低置信度边，支持权重自定义。
- **最短路径算法**：支持 Bellman-Ford、Dijkstra、SPFA、Floyd-Warshall 四种算法。
- **网络分析**：节点数、边数、平均度、连通分量、权重统计等。
- **可视化与导出**：支持网络可视化和GraphML格式导出，便于Cytoscape/Gephi等工具进一步分析。

## 2. 依赖

- pandas
- networkx
- matplotlib

请先安装依赖：

```bash
pip install pandas networkx matplotlib
```

## 3. 主要类与方法

### PPINetwork

#### 初始化

```python
PPINetwork(file_path=None, min_score=400, weight_type="inverse")
```

- `file_path`：STRING数据库的物理互作文件路径
- `min_score`：最小置信度阈值（默认400）
- `weight_type`：边权重类型，"inverse"（1-score）或"direct"（直接用score）

#### 主要方法

- `load_data(file_path, min_score=400, weight_type="inverse")`：加载数据
- `create_graph()`：构建NetworkX图
- `get_graph()`：获取NetworkX图对象
- `get_adjacency_list()`：获取邻接表
- `run_bellman_ford(source)`：运行Bellman-Ford算法
- `run_floyd_warshall()`：运行Floyd-Warshall算法
- `run_spfa(source)`：运行SPFA算法
- `run_dijkstra(source)`：运行dijkstra算法
- `run_parallel_dijkstra_adjlist()`：运行并行dijkstra算法
- `analyze_graph()`：输出网络基本属性
- `visualize(max_nodes=100, output_file=None)`：可视化网络
- `save_graph(output_file)`：保存为GraphML
- `export_for_visualization(output_file, source=None, target=None, ...)`：导出用于Cytoscape/Gephi的GraphML

### 最短路径算法类

- 普通方法共分为四类`BellmanFord(graph)`、`DijkstraHeap(graph)`、`SPFA(graph)`、`FloydWarshall(graph)`，执行PPINetwork的run_xx方法则返回一个该类的对象。
- 并行化的dijkstra算法单独作为一类`paralleldijkstra(graph)`，返回全对节点的最短路径矩阵（n*n）
- 统一接口：`shortest_paths_from(source)`、`get_shortest_path(target)`、`get_all_shortest_paths()`
- 对于FloydWarshall，接口为`get_path(u, v)`（获取两节点的最短路径），`get_distance(u, v)`（获取两节点的最短距离），`get_all_pairs_shortest_path(u, v)`）（获取所有节点对的最短路径）

---

# 测试程序示例

假设你有一个小型的蛋白质互作数据文件`test_ppi.txt`，内容如下（或用真实数据）：

```
protein1 protein2 combined_score
A B 800
B C 700
A C 900
C D 500
```

## test_ppi.py

```python
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from main import PPINetwork

# 1. 初始化并加载数据
ppi = PPINetwork(
    file_path="test_ppi.txt",
    min_score=400,
    weight_type="inverse"
)

# 2. 分析网络结构
stats = ppi.analyze_graph()
print("网络统计信息：", stats)

# 3. 运行 Bellman-Ford 算法
source, target = "A", "D"
bf = ppi.run_bellman_ford(source)
path = bf.get_shortest_path(target)
distance = bf.distances[target]
print(f"Bellman-Ford: {source} 到 {target} 的最短路径: {' → '.join(path)} (距离: {distance:.4f})")

# 4. 运行 Floyd-Warshall 算法
fw = ppi.run_floyd_warshall()
fw_path = fw.get_path(source, target)
fw_dist = fw.get_distance(source, target)
print(f"Floyd-Warshall: {source} 到 {target} 的最短路径: {' → '.join(fw_path)} (距离: {fw_dist:.4f})")

# 5. 运行 SPFA 算法
spfa = ppi.run_spfa(source)
spfa_path = spfa.get_shortest_path(target)
spfa_dist = spfa.dist[target]
print(f"SPFA: {source} 到 {target} 的最短路径: {' → '.join(spfa_path)} (距离: {spfa_dist:.4f})")
# 6. 运行dijkstra算法
print("\n运行 Dijkstra 算法:")
dijkstra = ppi_network.run_dijkstra(source_node)
dijkstra_path = dijkstra.get_shortest_path(target_node)
dijkstra_distance = dijkstra.distances[target_node]
print(f"从 {source_node} 到 {target_node} 的最短路径: {' → '.join(dijkstra_path)} (距离: {dijkstra_distance:.4f})")
# 7. 运行并行的dijkstra算法
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
#8.可视化网络
ppi.visualize(max_nodes=10, output_file="test_ppi.png")
print("网络可视化已保存为 test_ppi.png")

# 8. 导出GraphML
ppi.save_graph("test_ppi.graphml")
print("GraphML已导出为 test_ppi.graphml")
```

---

如需进一步自定义或扩展功能，请参考各类的源码注释。
