# EdMot_C
论文复现："EdMot: An Edge Enhancement Approach for Motif-aware Community Detection" (KDD 2019) [论文链接](https://arxiv.org/abs/1906.04560)
论文作者 MATLAB 版本实现 [here](https://github.com/lipzh5/EdMot_pro)
其他复现：基于 NetworkX 的 Python 版本实现 [here](https://github.com/benedekrozemberczki/EdMot)

个人觉得论文的重要工作在于提出在基本社区发现算法（如 Louvain, kmeans-Cluster）的结果图上作边缘强化。项目对论文进行部分复现，手撸了一些基本算法如 Louvain， Motif 矩阵的计算等等，并未充分考虑效率问题，好处是没有依赖可以直接下载运行。最后，在数据集 polbooks 上进行了简单验证。







