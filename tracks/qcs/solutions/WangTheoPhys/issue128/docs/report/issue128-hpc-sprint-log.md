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

| 角色 | Slurm job | 资源/上限 | 依赖 |
|---|---:|---|---|
| 31-stage exact D8 array | `23044178` | 每单元 17 CPU、64 GiB、4:30:00 | 无 |
| deep verifier + 全测试 | `23044196` | 9 CPU、32 GiB、4:30:00 | 无 |
| exact reducer | `23044212` | 26 CPU、96 GiB、1:15:00 | `afterok:23044178` |
| cancellation-aware exact D8 fallback | `23044280` | 63 CPU、240 GiB、4:15:00 | 无 |
| exact D6 early-result lane | `23044494` | 9 CPU、32 GiB、3:30:00 | 无；08:15 deadline |
| exact D6 extractor | `23044520` | 9 CPU、32 GiB、0:30:00 | `afterok:23044494`；08:15 deadline |
| exact D6 stage array | `23044609` | 每单元 9 CPU、32 GiB、2:30:00 | 无；31 路；08:15 deadline |
| exact D6 array reducer | `23044642` | 26 CPU、96 GiB、1:00:00 | `afterok:23044609`；08:15 deadline |

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
  RUNNING。`23044642` 仅在全阵列成功后做正序/逆序一致的精确归并。
- 远端 `23044196` 的 deep verifier 本身通过，但整作业以 pytest exit 1
  结束：五个 delivery-package 测试因远端同步目录没有可用 Git 元数据而失败，
  另一个既有慢测试的 symplectic 数值与当前发现值不一致。该失败被保留为
  环境/测试证据，不会冒充全绿；本地主证书的 fast/deep 与干净克隆复核仍独立通过。

## 制品契约

每个 `stage-XX` 目录包含：

- `shard.json.gz`：坐标寻址、registry 无关、canonical JSON、确定性 gzip；
- `manifest.json`：stage/order、开始结束时间、wall time、峰值 RSS、term counts、源提交、dirty 标志和 SHA-256；
- Slurm stdout/stderr：启动与完成状态。

归并器要求 31 个 manifest 全部成功、公式 ID/order/source commit 一致、每个 shard 哈希匹配，并要求正向与反向精确归并得到相同字典。任何失败均不得生成倍率主张。

## 终局判定

1. 若数组和 reducer 完成：读取 exact D8 以及同一单体制品中的 D6/D7 系数，评估下一轮分组空间并保留 hash-bound 证据。
2. 若任一 stage OOM/timeout：保留成功制品和 sacct 证据，按失败类别只重提缺失单元，不把部分和当全局界。
3. 五倍仍需新的 D4 紧化：当前 D4 单项在 `r=78` 已约为 `1.259084e-6`，所以 D8 再小也不能在当前非负误差账本中闭合五倍。
4. D8 结果只用于诊断和设计下一份证书；当前主倍率保持新近验证的 4.1357418449×。
