# Issue 128 六小时 HPC 冲刺运行记录

记录日期：2026-07-31（Asia/Shanghai）

## 主张演进与边界

本轮没有覆盖既有十文件冻结包；该包仍独立复现
`11791/2911 = 4.050498110614909...`，对应 `r=393 → r=97`。

冲刺期间，已有的 exact D5 sidecar 被正式接入一个新的主证书。新证书
`issue128-d5-integrated-certificate.json` 接受 `r=95`、拒绝 `r=94`，并已
同时通过 fast 与 deep verifier。因此本轮新增且可报告的倍率是
`11791/2851 = 4.135741844966678...`。两份证书各自保留，哈希与 provenance
不混写。D8 分项仍不得单独宣传为五倍整体结果。

## 可复现源状态

- 本地分支：`codex/issue128-hpc-six-hour`
- 生产源提交：`e65c95f4ae76dc1a78e81715cac21cfd7a81dbf8`
- 远端入口：`xh5-acamtw`
- 远端仓库：`/work/home/acamtw70yu/quantum-harness-issue128-hpc`
- Slurm account/QOS/partition：`giggleliu` / `user_acamtw70yu` / `xhacnormalb`
- 计算类型：CPU 上的 Python 精确整数/有理数与稀疏 Pauli 字典；GPU 不适用。

## 快速门槛

- Python 文件通过 `py_compile`。
- 坐标无关制品、确定性 gzip、manifest 与损坏拒绝测试：3 passed。
- 本地 stage 30 / order 8 smoke 成功，重复制品哈希稳定。
- 本地常规测试：93 passed，11 deselected，41.12 秒。
- 本地 deep verifier：`valid=true`，`deep_proof_regenerated=true`。
- 远端 smoke：作业 `23044170`，stages 28--30 完成且 manifest/hash 正常；stage 27 用作较重梯度检查。

## D5-integrated 本地闭环

- D5 sidecar：605,832 个规范 Pauli 项，123,106 个同支撑两两反对易组，最大组大小 10。
- 新边界：`r=95` 全局上界 `9.90449140283245e-7`，`r=94` 上界 `1.04680611656037e-6`。
- 新资源：2,851 个 merged groups、205,272 个 bond propagators、615,816 个 CNOT 上界。
- 精确倍率：`11791/2851 = 4.135741844966678...`。
- fast verifier：通过。
- deep verifier：`deep_proof_regenerated=true`，167.62 秒，最大 RSS 1,533,149,184 bytes。
- 常规测试：100 passed、11 deselected；新增最小性与篡改拒绝测试 3 passed。

这一改善不依赖正在运行的 D8 任务，因此 D8 的成功、超时或 OOM 都不会撤销
4.1357418449× 的新证书。

## 生产图

`23044178` 的 25 个重型阵列单元在 2026-07-31 06:51 CST 仍处于 RUNNING，
并最终分别在 `4:30:11`--`4:30:13` 达到墙钟上限；因此本轮“生产计算至少
4.25 小时”的时间要求已由 `sacct` 实际记录满足，而不是用申请的 time limit
代替运行时长。其余 6 个单元正常完成；父 array job 结束后，要求全阵列成功的
`afterok` reducer `23044212` 被 Slurm 取消，没有从部分和生成全局结果。

| 角色 | Slurm job | 资源/上限 | 依赖 |
|---|---:|---|---|
| 31-stage exact D8 array | `23044178` | 每单元 17 CPU、64 GiB、4:30:00 | 无 |
| deep verifier + 全测试 | `23044196` | 9 CPU、32 GiB、4:30:00 | 无 |
| exact reducer | `23044212` | 26 CPU、96 GiB、1:15:00 | `afterok:23044178` |
| cancellation-aware exact D8 fallback | `23044280` | 63 CPU、240 GiB、4:15:00 | 无 |
| exact D6 early-result lane | `23044494` | 9 CPU、32 GiB、3:30:00 | 无；08:15 deadline |
| exact D6 extractor | `23044520` | 9 CPU、32 GiB、0:30:00 | `afterok:23044494`；08:15 deadline |
| exact D6 stage array | `23044609` | 每单元 9 CPU、32 GiB、2:30:00 | 无；31 路；08:15 deadline |
| exact D6 array reducer | `23044651` | 26 CPU、96 GiB、1:00:00 | `afterok:23044609`；08:32 deadline |
| D6/D8 cross-order check | `23044659` | 17 CPU、64 GiB、0:30:00 | stages 25--30；已完成 |
| direct exact D6 identity | `23044709` | 17 CPU、64 GiB、1:45:00 | 08:25 deadline |
| exact E7 word array | `23044736` | 每单元 17 CPU、64 GiB、1:00:00 | 16 路；08:20 deadline |
| E7-to-D6 exact reducer | `23044752` | 26 CPU、96 GiB、0:35:00 | `afterok:23044736`；08:34 deadline |
| exact E7 CPU fanout | `23044754` | 16 单元 × 17 进程；每单元 17 CPU、64 GiB、0:45:00 | 272 路；08:10 deadline |
| fanout E7-to-D6 reducer | `23044907` | 26 CPU、96 GiB、0:40:00 | `afterok:23044754`；替代未运行的 `23044755` |
| exact E7 fine fanout | `23044779` | 32 单元 × 17 进程；每单元 17 CPU、64 GiB、0:40:00 | 544 路；08:08 deadline |
| fine-fanout E7-to-D6 reducer | `23044906` | 26 CPU、96 GiB、0:45:00 | 替代 fail-closed 的 `23044780`；08:34 deadline |
| exact E7 ultra-fine fanout | `23044815` | 64 单元 × 17 进程；每单元 17 CPU、64 GiB、0:30:00 | 1,088 路；08:05 deadline |
| ultra-fine E7-to-D6 reducer | `23044816` | 26 CPU、96 GiB、0:24:00 | `afterok:23044815`；08:34 deadline |

数组最初以 16 路并发启动；确认账号/QOS允许后，用 Slurm 正常调度接口把 `ArrayTaskThrottle` 提升到 31。所有剩余单元随后进入 RUNNING；调度和授权仍由集群控制器执行。

另一路单遍精确递推在每个 stage 后先合并所有贡献的相同 Pauli 项，再继续共轭，用于降低 31 个独立分片在早期 stage 的重复工作。其本地 order-4 结果为 74,448 个精确项；冻结区间侧车的 75,324 项包含依赖信息丢失后保留的 ghost terms，因此不是矛盾。更强的 order-5 检查得到 605,832 项，与冻结精确 D5 侧车完全一致；本地耗时 1385.7 秒，峰值 RSS 约 1.47 GiB。

24 个 matching-order discovery 单元 `23044214` 运行 34 分钟后仍为 0 complete。为避免它们在一小时上限整体 timeout，并把调度窗口留给精确 D8 fallback，该非可信数组被取消。它们的 Slurm 日志和 running manifest 被保留；未把部分状态当作排序结果。

在 D5-integrated 证书闭合后，又提交了一条 order-6 单体路径。它复用同一
exact-cubic recurrence，但只保留到 D6，因此预期早于 D8 返回，并把 D0--D6
全系列写入 canonical gzip。该作业通过 `sbatch --test-only`，正式 job
`23044494` 随即进入 RUNNING；设置 2026-07-31 08:15 硬截止，避免在六小时
冲刺窗口之后留下计算。若其 exact D6 Pauli-l1 或严格分组界优于现用通用
D6 majorant，它只能在独立 artifact 校验与 fast/deep 接入后形成下一次倍率更新。

随后在本地利用 matching 的平移/旋转重标号对称性做了降维复核。24 个排列的 order-4 指标严格形成三个 8 元轨道，代表 permutation indices 为 0、2、3。冻结顺序 index 0 的 D4/D5 coefficient-l1 分别为约 20.160966/95.679103；另外两轨道为约 20.664900/100.103 与 20.664900/100.125。冻结顺序在两个已测阶数均占优，因此停止 order-6 permutation 路线，不把更多算力投向较差候选。

### 运行中审计更新（2026-07-31 05:20 CST）

- D8 数组 stage 25 已成功完成：D8 精确项数 17,211,714，D0--D8 项数为
  `[6, 72, 582, 3744, 22212, 128616, 710298, 3650880, 17211714]`；墙钟
  7,349.722 秒，峰值 RSS 52,152,229,888 bytes，输出 SHA-256 为
  `c64c0729a1eff25b771d2ccb94170648aa4235a3ac7f40366201a95a4eb3cf7b`。
- `23044494` 的 exact D6 路径已完成 stage 6 前的递推并持续运行；最近一次
  manifest 进度显示 D6 项数从 152,934、448,302、1,183,386 增长到
  3,244,176，36 分 35 秒时最大 RSS 约 5.79 GiB。
- `23044520` 将在 D6 主作业成功后提取独立、canonical-gzip 的 D6-only
  sidecar；只有该 sidecar 的压缩哈希、父制品摘要、精确 l1 和本地 verifier
  全部通过，才允许更新主倍率。
- 单体 D6 到 stage 8 时累计墙钟约 65.8 分钟、D6 项数 4,654,758；按后半程
  增长估计，3:30 上限存在风险。因此在不取消单体与 D8 生产的前提下，补交
  31 路 order-6 stage contribution 阵列 `23044609`。它复用已通过 smoke 的
  hash-bound shard 契约，只把截断阶数从 8 改为 6；所有 31 个单元立即进入
  RUNNING。`23044651` 仅在全阵列成功后做正序/逆序一致的精确归并。原占位
  `23044642` 因 08:15 deadline 无法在线延后而在 PENDING 状态取消；它从未运行、
  未产生输出。替代作业沿用同一 `afterok:23044609` 依赖，并把硬截止设为 08:32。
- 远端 `23044196` 的 deep verifier 本身通过，但整作业以 pytest exit 1
  结束：五个 delivery-package 测试因远端同步目录没有可用 Git 元数据而失败，
  另一个既有慢测试的 symplectic 数值与当前发现值不一致。该失败被保留为
  环境/测试证据，不会冒充全绿；本地主证书的 fast/deep 与干净克隆复核仍独立通过。
- `23044659` 在 15 分 49 秒内逐项比较了两条独立生产路径。D6 阵列与 D8
  阵列在 stages 25--30 的 degree-6 规范 Pauli 系数图完全相同；比较覆盖的
  单 stage D6 项数依次为 710,298、446,094、150,858、2,028、48、0。
  比较对象是精确三次域系数映射，不只是最终 l1 小数。
- 为规避早编号 stage contribution 的长共轭链，又增加了直接恒等式路径
  `23044709`。它从同一 31-stage 公式的 16,380 个七阶自由词构造 exact E7，
  然后使用 `D6 = 7 E7 + (2/3) ad_A^2(E5)` 形成 exact D6。该路线与 stage
  分解算法不同；输出仍采用同一 canonical coordinate/cubic-field sidecar 契约，
  只有与阵列归并结果逐项一致后才可作为独立复核。
- 为消除直接 E7 枚举中的串行瓶颈，在 07:13 CST 提交 16 路、按自由词索引
  精确分片的生产阵列 `23044736`，源提交为
  `23dbe9532929d81c8cbb366e12a3d86c03bd4e9e`。提交前远端
  `cubic_local.py` 与本地 SHA-256 均为
  `dbe33aec6ca5835c4db46f38fe403817c635ce867c4de38f6dbcc04d47e82d1b`；
  16 个单元均立即进入 RUNNING。依赖归并作业 `23044752` 要求 16 份 manifest
  完整、词区间无缝覆盖、payload hash 正确、正反向精确合并一致，随后才用
  `D6 = 7 E7 + (2/3) ad_A^2(E5)` 写出 D6 sidecar；任何单元失败都会阻止归并。
- 16 路车道首分钟暴露出负载不均：3 个连续索引块约 43--45 秒完成且 E7
  贡献为零，其余单元停留在高成本词；同时 `sstat` 显示单进程 RSS 仅约
  165--203 MiB，远低于 64 GiB 分配。于是增加 CPU 饱和车道 `23044754`：
  仍为 16 个 Slurm array cells，但每个 cell 通过 `srun` 启动 17 个单核
  精确 evaluator，共把 16,380 个词细分为 272 个不重叠区间。16 个 cells
  均在提交后 11 秒内进入 RUNNING，272 份 running manifest 全部出现；抽样
  `sstat` 显示一个 cell 的 step 有 17 tasks 且 CPU 时间随墙钟增长。依赖
  reducer `23044755` 使用同一严谨归并契约，不能读取 16 路车道的目录。
- 272 路在约 5 分钟时达到 155 complete / 117 running / 0 failed，说明细分
  有效，但未完成集合仍集中于高成本词。窗口允许的情况下再提交 544 路独立
  车道 `23044779`，把每片缩到约 30 个词；32 个 cells 同样在 11 秒内全部
  RUNNING，并建立 544 份 manifest。其 reducer `23044780` 与 272 路目录隔离；
  如果两条车道均闭合，最终 D6 必须逐项相同，否则任何一方都不得进入证书。
- 归并预检发现初版 Slurm wrapper 把 `--shard-count` 写死为 16。所有 fanout
  父作业当时仍在 RUNNING、两个 reducer 仍为 dependency PENDING，因此在任何
  reduction 开始前把 wrapper 改为读取 `ISSUE128_SHARD_COUNT`（缺省仍为 16），
  通过 `bash -n` 后同步到远端并核对实际脚本文本。272/544 reducer 的提交环境
  分别已导出 272/544；该修复提交为
  `08643023ac62b98e2d2d97ef82cb3fb616710b76`。
- Slurm 会在提交时固化 batch script，而不是在作业启动时重新读取远端文件。
  因此修复同步前已经提交的 544 reducer `23044780` 仍以旧参数启动，并在 2 秒内
  因 shard 0 的 `shard_count=544` 与旧期望 16 不符而 fail-closed；输出目录为空，
  没有形成部分结果。修复后重提的 `23044906` 使用 45 分钟上限并立即 RUNNING。
  同理，尚未运行的旧 272 reducer `23044755` 被取消并以 `23044907` 替代，仍
  保留 `afterok:23044754`。1088 reducer `23044816` 是修复后提交，已越过元数据
  检查并持续归并。
- 544 路在约 4 分钟达到 383 complete / 161 running / 0 failed，但已完成
  高成本片的单进程峰值升至约 3 GiB。于是提交第三条 1,088 路车道
  `23044815`，每片约 15 个词，64 个 cells 全部进入 RUNNING；其独立 reducer
  为 `23044816`。三种宽度共同提供“最快返回者”和跨分片边界逐项复核，且
  任一 reducer 都只能消费自己的完整目录。
- 07:28:33 CST 的一次集中快照显示当前用户共有 1,748 个 Slurm CPU 处于
  RUNNING allocation。三条 fine-grained E7 车道分别为 242/272、482/544、
  567/1,088 complete，均为 0 failed；对应已完成分片的最大 RSS 约为
  6.41、4.15、1.33 GiB。该证据说明加密分片同时提高了真实 CPU 并发并降低
  单进程内存与最长尾，而不是只扩大申请数字。

## 制品契约

每个 `stage-XX` 目录包含：

- `shard.json.gz`：坐标寻址、registry 无关、canonical JSON、确定性 gzip；
- `manifest.json`：stage/order、开始结束时间、wall time、峰值 RSS、term counts、源提交、dirty 标志和 SHA-256；
- Slurm stdout/stderr：启动与完成状态。

归并器要求 31 个 manifest 全部成功、公式 ID/order/source commit 一致、每个 shard 哈希匹配，并要求正向与反向精确归并得到相同字典。任何失败均不得生成倍率主张。

## 终局运行审计（2026-07-31 08:22 CST）

- D8 stage array `23044178`：6 complete、25 timeout。25 个重型单元实测
  `4:30:11`--`4:30:13`；依赖 reducer 未运行。D8 monolithic `23044280`
  在 `4:15:29` timeout。两条路径都没有全局 D8 制品。
- D6 monolithic `23044494` 在 `3:30:05` timeout，extractor 未运行。D6 stage
  array `23044609` 为 13 complete、18 timeout，依赖 reducer 未运行；已完成的
  stages 25--30 与 D8 路径的 exact D6 map 逐项一致证据仍有效，但部分数组不能
  代替完整 D6。
- E7 原始 16 路 `23044736` 为 4 complete、12 timeout。CPU 饱和的
  `23044754`、`23044779`、`23044815` 三条父阵列则分别 16/16、32/32、64/64
  cells 全部完成，对应 272、544、1,088 份 word manifests 全覆盖、0 failed。
- 三条完整 E7 目录的 D6 reducer 分别在 `24:21`、`45:06`、`40:07` 达到
  墙钟上限；它们都遵守原子写，`reduced/` 中没有可误用的半成品。旧 544
  reducer 的 2 秒 metadata fail-closed 和 replacement 过程已在上文单独记录。
  因没有两份完整 D6，三个 exact-map compare 依赖均未运行。
- 串行 direct identity `23044709` 在 `1:39:06` timeout，最后可见进度为
  16,380 个七阶词中的 1,792 个；没有输出完整 D6。08:24:40 CST 复核时当前
  用户队列为空，没有作业越过本次冲刺边界。
- 本地最终常规回归为 112 passed、11 deselected、175.21 秒。新增比较器的
  正例与系数篡改拒绝测试均通过。

终局判定是严格的负操作结果：本轮显著提高了并发和精确分片覆盖，也证明了
4.25 小时生产下限，但没有在六小时内得到完整 D6/D8 reduction。因此不更新
步数、倍率或五倍主张；当前可报告结果仍是已经 fast/deep 验证的
`11791/2851 = 4.135741844966678...`。

## 判定规则

1. 若数组和 reducer 完成：读取 exact D8 以及同一单体制品中的 D6/D7 系数，评估下一轮分组空间并保留 hash-bound 证据。
2. 若任一 stage OOM/timeout：保留成功制品和 sacct 证据，按失败类别只重提缺失单元，不把部分和当全局界。
3. 五倍仍需新的 D4 紧化：当前 D4 单项在 `r=78` 已约为 `1.259084e-6`，所以 D8 再小也不能在当前非负误差账本中闭合五倍。
4. D8 结果只用于诊断和设计下一份证书；当前主倍率保持新近验证的 4.1357418449×。
