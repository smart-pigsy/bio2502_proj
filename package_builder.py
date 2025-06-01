import os
import shutil

# 创建包目录结构
def setup_package():
    base = os.path.dirname(os.path.abspath(__file__))
    pkg = os.path.join(base, 'ppinetwork')
    os.makedirs(pkg, exist_ok=True)
    os.makedirs(os.path.join(pkg, 'test_data'), exist_ok=True)
    # 拷贝源码
    for fname in ['Bellman_Ford.py', 'dijkstra.py', 'Floyd_Warshall.py', 'SPFA.py', 'paralleldijkstra.py', 'main.py']:
        shutil.copy(os.path.join(base, fname), os.path.join(pkg, fname.lower()))
    # 拷贝测试数据
    shutil.copy(os.path.join(base, 'test_data', '7209.protein.physical.links.v12.0.txt'), os.path.join(pkg, 'test_data', '7209.protein.physical.links.v12.0.txt'))
    # 拷贝文档
    shutil.copy(os.path.join(base, 'ppinetwork', 'manual.md'), os.path.join(pkg, 'manual.md'))
    shutil.copy(os.path.join(base, 'ppinetwork', 'README.md'), os.path.join(pkg, 'README.md'))
    # 生成__init__.py
    with open(os.path.join(pkg, '__init__.py'), 'w', encoding='utf-8') as f:
        f.write('# PPI Network package init\n')
    print('Package structure created at', pkg)

if __name__ == '__main__':
    setup_package()
