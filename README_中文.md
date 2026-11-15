# 巴西和墨西哥高温暴露与尿石症住院

**两国逐日病例交叉研究：急性医疗服务负担及胆道—胰腺疾病比较**

代码版本：**v1.0.0**  
结果／图件版本：**SRR-20261114-C3-display-3**  
作者：Shuo Feng、Chao Shi、Tianming Liu、Haoran Ding、Fenghe Hang、Zhenghao Huang、Xianzhi Meng。

本仓库保存最终候选稿的封存汇总结果，并据此重建五幅图和三张主要结果表。默认构建过程不下载原始数据、不重建分析面板，也不重新拟合统计模型。请先阅读[重建说明](docs/reproduction.md)和[范围与局限](docs/limitations.md)，再判断本仓库能支持何种复现。

## 重建候选稿图表

绘图环境测试使用 Python 3.12.3。激活 Python 环境后，安装图表依赖并运行：

```sh
python -m pip install -r requirements-figures.txt
python reproduce.py --verify-only
python reproduce.py
```

`--verify-only` 按文件字节数和 SHA-256 校验清单中声明的输入文件。默认命令先校验这些文件，再重建图和表；不会重新拟合模型。生成文件写入 `artifacts/`，该目录由 Git 忽略。输出及操作细节见[重建说明](docs/reproduction.md)。

## 仓库内容

- `data/aggregate_results/`：补充汇总结果，以及 `main_tables/table1.csv` 至 `table3.csv` 三张主要结果表。
- `data/figure_inputs/`：五幅图所用的冻结汇总结果和修订阶段输入。
- `manifest.csv`：仓库相对路径、SHA-256、字节数和文件角色。
- `src/`：由 `reproduce.py` 调用的绘图与主要结果表构建脚本。
- `data/processed_inputs/`：六份地区或州逐日汇总输入，含无损拆分的巴西面板；见[输入说明](docs/processed-inputs.md)。
- `analysis/work/upgrade_20261112/` 和 `analysis/work/revision_20261113/`：归档的模型与分析脚本；它们不是受支持的汇总结果重建流程。
- `reference/figures/`：如有收录，则为可选参考 PDF；构建脚本不会读取这些文件。

整理包包含 121 个补充汇总文件、额外一份历史冻结参考表、19 份绘图输入副本；三张主表 CSV 已包含在汇总文件中。[可选模型入口](docs/model_runner.md)可在隔离目录准备服务负担与非线性分析，已检查输入与路径转换，但本次未重新拟合模型。这不代表完整的原始数据处理流程。

## 数据、代码与使用

数据来源类别包括巴西 SIH/DATASUS 和 CNES、墨西哥 DGIS/SAEH SSA、ERA5-Land、NASA POWER、IPCA，以及 Hu 等人（2025）描述的 GFED 驱动火灾 PM2.5 浓度产品（DOI：[10.5194/essd-17-3741-2025](https://doi.org/10.5194/essd-17-3741-2025)）。[数据来源说明](docs/data-sources.md)区分了已命名的数据来源与经过核实的下载记录。

作者编写的代码及配套代码文档采用[MIT许可证](LICENSE)。第三方数据、来源派生输入和软件依赖仍遵循其适用条款，详见[LICENSE_STATUS.md](LICENSE_STATUS.md)。

## 引用

请使用 [CITATION.cff](CITATION.cff) 中的引文元数据。仓库地址为https://github.com/FENG1567/brazil-mexico-temperature-urolithiasis，代码发布版本为`v1.0.0`，未分配DOI。

## 图表展示修订

本次仅更新五幅图和三张主表的展示，未重新拟合模型。修订前三张主表的完整内容保留在 `data/aggregate_results/full_model_tables/`，五组次要结果CSV见 `display_tables/`。

## 代码可用性

代码及配套汇总输入公开于https://github.com/FENG1567/brazil-mexico-temperature-urolithiasis（版本v1.0.0）。仓库包含分析脚本、处理后的地区–日输入，以及从已发布汇总结果重建五张主图和三张主表的操作说明。作者编写的代码采用MIT许可证；第三方数据仍遵循各提供方的使用条款。模型重新运行的范围与前提条件见仓库说明。
