# Issue 128：六小时超算部署与收件手册

这份手册把本轮实测有效的流程固化为可重复操作。目标不是让 Slurm 里出现
最大的申请数字，而是在六小时内最大化**通过验证的精确制品**。当前可信
headline 仍是 `11791/2851 = 4.135741844966678...`；任何 D6/D8 新倍率都必须
等完整 ledger、相邻步边界、fast/deep verifier 和制品哈希全部通过后再发布。

## 1. 六小时墙钟切片

| 时段 | 上限 | 退出条件 |
|---|---:|---|
| 本地准备 | 10 分钟 | 语法、focused tests、源哈希通过 |
| 远端 smoke | 10 分钟 | 1 个小制品、manifest、SHA-256 通过 |
| 主生产 | 至少 4 小时 15 分 | 生产任务实测墙钟达到 4:15，或完整结果提前返回 |
| 归并/验证/收件 | 最多 1 小时 25 分 | reducer、比较、证书、报告全部落盘 |

不要把 `#SBATCH --time` 当成实际生产时长。生产下限必须由 `sacct` 的
`Elapsed` 或仍处于 `RUNNING` 的调度快照证明。本轮 D8 重型单元实测达到
`4:30:11`--`4:30:13`，因此 4.25 小时要求已由运行记录满足。

## 2. 十分钟内完成本地准备

在本地仓库根目录执行：

```bash
ISSUE=tracks/qcs/solutions/WangTheoPhys/issue128
git status --short
git rev-parse HEAD
bash -n hpc/issue128_e7_word_fanout.sbatch
bash -n hpc/issue128_e7_word_reduce.sbatch
PYTHONPATH="$ISSUE/src:$ISSUE" python -m py_compile \
  "$ISSUE/scripts/build_e7_word_shard.py" \
  "$ISSUE/scripts/reduce_e7_word_shards.py"
PYTHONPATH="$ISSUE/src:$ISSUE" pytest -q \
  "$ISSUE/tests/test_hpc_artifacts.py" \
  "$ISSUE/tests/test_compare_exact_degree_payloads.py"
```

只同步任务需要的文件，不使用 `--delete`，也不覆盖远端结果目录：

```bash
REMOTE=xh5-acamtw:/work/home/acamtw70yu/quantum-harness-issue128-hpc
rsync -a hpc/ "$REMOTE/hpc/"
rsync -a "$ISSUE/src/" "$REMOTE/$ISSUE/src/"
rsync -a "$ISSUE/scripts/" "$REMOTE/$ISSUE/scripts/"
```

同步后对关键源码做本地/远端 SHA-256 对照。source commit 与关键文件摘要要写入
运行日志；不要用“目录大致相同”代替 provenance。

## 3. 把 smoke 限制在十分钟内

先运行 `bash -n` 和 `sbatch --test-only`。随后只在计算节点上计算一个七阶词：

```bash
ssh xh5-acamtw 'cd /work/home/acamtw70yu/quantum-harness-issue128-hpc && \
  sbatch --partition=xhacnormalb --account=giggleliu --qos=user_acamtw70yu \
  --cpus-per-task=5 --mem=16G --time=00:10:00 \
  --output=results/slurm-issue128-e7-smoke-%j.out \
  --wrap="ISSUE=tracks/qcs/solutions/WangTheoPhys/issue128; \
  PYTHONPATH=\$ISSUE/src .venv/bin/python -u \
  \$ISSUE/scripts/build_e7_word_shard.py --shard-index 0 --shard-count 16380 \
  --output results/e7-smoke/shard.json.gz \
  --manifest results/e7-smoke/manifest.json"'
```

只检查：作业 exit 0、manifest 为 complete、payload SHA 与 manifest 相符、重复写
哈希稳定。smoke 失败就按认证、资源门槛、代码异常三类记录；本轮不在登录节点
直接运行计算，也不无上限反复试错。

## 4. 生产提交：先 544 路，必要时升到 1,088 路

此算法是稀疏 Python 精确整数/有理数运算，GPU 不适用。集群按约
3931 MiB/CPU 约束内存，因此每个 Slurm cell 申请 17 CPU/64 GiB，并用
`srun` 真正启动 17 个单核 evaluator，避免为内存申请的 CPU 闲置。

推荐从 32 cells × 17 workers = 544 shards 开始。以下整段在
`ssh xh5-acamtw` 登录后的远端 shell 中执行：

```bash
export ISSUE128_RUN_ID=issue128-prod-YYYYMMDD-HHMM-e7-fanout544
# 粘贴准备阶段在本地记录的完整 40 位 commit；不要依赖远端同步目录的 .git。
export ISSUE128_SOURCE_COMMIT=LOCAL_FULL_COMMIT
export ISSUE128_WORKERS_PER_TASK=17
export ISSUE128_ARRAY_CELL_COUNT=32
export ISSUE128_SHARD_COUNT=544
cd /work/home/acamtw70yu/quantum-harness-issue128-hpc

ARRAY_JOB=$(sbatch --parsable --array=0-31%32 --time=00:40:00 \
  --export=ALL,ISSUE128_RUN_ID,ISSUE128_SOURCE_COMMIT,ISSUE128_WORKERS_PER_TASK,ISSUE128_ARRAY_CELL_COUNT,ISSUE128_SHARD_COUNT \
  hpc/issue128_e7_word_fanout.sbatch)
```

提交 reducer 前必须确认远端 wrapper 使用
`--shard-count "$ISSUE128_SHARD_COUNT"`，不能写死 16。Slurm 在提交时固化
batch script；修改远端文件不会改变已经 pending 的旧作业。

```bash
REDUCE_JOB=$(sbatch --parsable --dependency=afterok:"$ARRAY_JOB" --time=00:45:00 \
  --export=ALL,ISSUE128_RUN_ID,ISSUE128_SHARD_COUNT \
  hpc/issue128_e7_word_reduce.sbatch)
```

如果 5 分钟快照显示单片 RSS 接近 3--4 GiB、完成数停滞或仍有明显长尾，再提交
64 cells × 17 workers = 1,088 shards 的独立 run。不要写入同一 `RUN_ID`；不同
宽度的 reducer 必须使用目录隔离，完成后逐项比较 exact D6 map。

## 5. 低扰动监控

每 60--120 秒取一次紧凑快照，避免高频 `squeue/sstat` 给控制器施压：

```bash
ssh xh5-acamtw 'squeue -h -j JOB_IDS -o "%i|%T|%M|%L|%R"'
ssh xh5-acamtw 'sacct -n -X -j JOB_IDS \
  --format=JobIDRaw,State,Elapsed,ExitCode -P'
ssh xh5-acamtw 'sstat -n -P -j REDUCER.batch \
  --format=JobID,AveCPU,AveRSS,MaxRSS,MaxVMSize'
```

数组成功条件是所有 manifest 完整且 `failed=0`。Reducer 还必须验证：词区间从
0 到 16,380 无缝覆盖、source commit 单一、payload SHA 全部匹配、正序/逆序
精确求和相同。任何一个分片缺失、timeout、OOM 或 hash 错误，都不能把部分和
写成全局 D6。

## 6. 收件与证书门槛

只收 reducer 的三件套：`degree-6.json.gz`、`word-parent.json`、`summary.json`，
以及 Slurm 日志。先在远端确认摘要与 SHA，再 `rsync` 到本地。两条 lane 都完成
时执行：

```bash
PYTHONPATH="$ISSUE/src:$ISSUE" python \
  "$ISSUE/scripts/compare_exact_degree_payloads.py" \
  --left PATH_A/degree-6.json.gz \
  --right PATH_B/degree-6.json.gz \
  --degree 6 --output PATH_COMPARISON/comparison.json
```

只有 `exact_coefficient_map_equal=true` 才进入证书构建：

```bash
cp VERIFIED/degree-6.json.gz "$ISSUE/certificates/issue128-d6-exact.json.gz"
cp VERIFIED/word-parent.json "$ISSUE/certificates/issue128-d6-parent.json"
PYTHONPATH="$ISSUE/src:$ISSUE" python \
  "$ISSUE/scripts/build_d6_integrated_certificate.py"
PYTHONPATH="$ISSUE/src:$ISSUE" python \
  "$ISSUE/scripts/verify.py" \
  "$ISSUE/certificates/issue128-d6-integrated-certificate.json"
```

最后再运行完整常规回归、deep regeneration、SHA256SUMS 和干净 worktree 复核。
如果 deep 未在窗口内完成，就把 D6 明确标为 provisional，并继续使用已认证的
4.1357418449 倍；绝不能用“算出了一个更小数字”替代完整可验证证书。
