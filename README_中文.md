# 巴西和墨西哥短期气温与泌尿系结石住院的关联：巴西胆石症和胆囊炎的探索性高温—烟霾比较

v1.1.1仅恢复论文题目和图1中已确认的A区纵向位置。科学结果与暴露输入沿用v1.1.0。



代码版本：**v1.1.1**  
修正结果标识：**TEMP-20261115-C1**  
历史版本：**v1.0.0** 保持原样并继续保留。  
作者：Shuo Feng、Chao Shi、Tianming Liu、Haoran Ding、Fenghe Hang、Zhenghao Huang、Xianzhi Meng。

v1.1.0 发布包包含修正后的巴西气温面板、固定网格权重和结果标识 TEMP-20261115-C1。五幅图和三张主表从本包提供的汇总结果重建；该构建流程不下载原始数据、不重建天气输入，也不拟合模型。请阅读[气温输入修正说明](docs/temperature_input_correction.md)、[重建说明](docs/reproduction.md)和[范围与局限](docs/limitations.md)，了解支持的范围和前提。

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
- `data/processed_inputs/`：地区或州逐日汇总输入，包括保留的历史巴西面板和紧凑格式的修正面板；见[输入说明](docs/processed-inputs.md)及[气温输入修正说明](docs/temperature_input_correction.md)。
- `analysis/work/upgrade_20261112/` 和 `analysis/work/revision_20261113/`：归档的模型与分析脚本；它们不是受支持的汇总结果重建流程。
- `reference/figures/`：如有收录，则为可选参考 PDF；构建脚本不会读取这些文件。

[可选模型入口](docs/model_runner.md)可在隔离目录准备服务负担与非线性分析。v1.1.0 需通过 [气温输入修正说明](docs/temperature_input_correction.md)中的 --brazil-panel 参数显式选择修正面板。仅准备阶段不会拟合模型；完整的上游原始天气重建不属于本包范围。

## 数据、代码与使用

数据来源类别包括巴西 SIH/DATASUS 和 CNES、墨西哥 DGIS/SAEH SSA、ERA5-Land 与 ERA5 单层再分析（用于修正巴西气温面板）、NASA POWER、IPCA，以及 Hu 等人（2025）描述的 GFED 驱动火灾 PM2.5 浓度产品（DOI：[10.5194/essd-17-3741-2025](https://doi.org/10.5194/essd-17-3741-2025)）。见[数据来源说明](docs/data-sources.md)和[气温输入修正说明](docs/temperature_input_correction.md)，了解下载限制与修正范围。

作者编写的代码及配套代码文档采用[MIT许可证](LICENSE)。第三方数据、来源派生输入和软件依赖仍遵循其适用条款，详见[LICENSE_STATUS.md](LICENSE_STATUS.md)。

## 引用

引文元数据见 [CITATION.cff](CITATION.cff)。本次准备的修正代码版本为 v1.1.0，结果标识为 TEMP-20261115-C1；历史 v1.0.0 保持原样。尚未分配 DOI。

## 图表展示修订

v1.1.0 的五幅图和主表从本次发布包内已修正的汇总结果生成。构建步骤只负责呈现这些结果；模型运行是独立步骤，且需要满足仓库说明的输入与软件前提。

## 代码可用性

本次准备的 v1.1.0 发布包（含修正后的巴西地区–日面板和固定网格权重）对应仓库地址为 https://github.com/FENG1567/brazil-mexico-temperature-urolithiasis（结果标识 TEMP-20261115-C1）。仓库说明如何从汇总结果重建五幅图和三张主表，并介绍单独选择的模型重运行流程。作者编写的代码采用 MIT 许可证；第三方数据仍遵循各提供方条款。历史版本 v1.0.0 继续作为原始版本保留。
