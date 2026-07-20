# Full 3-D PO/IPO versus C8 modal PO/IPO

固定コミット
[`d37311da872406c1a439f1c184f12d80fb8269d7/I1_オクタゴン`](https://github.com/kouske2222/concentrator_2/tree/d37311da872406c1a439f1c184f12d80fb8269d7/I1_%E3%82%AA%E3%82%AF%E3%82%BF%E3%82%B4%E3%83%B3)
のFortranコードを基礎に、開放端を持つ先すぼまりPEC八角形について、full 3-DとC8 DFTモード分解を同じパネル・核・観測点で比較する独立プロジェクトです。元リポジトリは変更していません。

## 実行

GNU Fortran、OpenMP、Python 3、NumPy、Matplotlibを使います。

```bash
make debug
./build/compare_full_vs_c8 config_validation.nml
python3 plot_compare.py --result results/validation
```

通常の最適化ビルドは次です。

```bash
make clean
make -j
OMP_NUM_THREADS=8 ./build/compare_full_vs_c8 config_validation.nml
make plot
```

出力は `results/validation` に保存されます。`make debug` は境界検査、浮動小数点例外、警告を有効にします。

Full参照は`use_dense_full`でdense／matrix-freeを切り替えます。両者は同じ`pair_map_local`を使用します。94 GHzの135条件掃引は、まず次の安全な事前見積もりを実行します。

```bash
./build/run_94ghz_sweep config_94ghz.nml
```

同じ周波数では形状依存の伝搬演算子を偏波角・軸ずれ・入射角の間で再利用できます。周波数が変わればGreen関数が変わるため再構築が必要です。

## 元コードから継承したもの

- `mod_types`、`mod_config`、`mod_geometry`、`mod_incident`、電場積分、IPO反復というモジュール分割
- paraxial Gaussian beamの局所波数方向と、x偏波を横波へ射影するベクトル計算
- \(\mathbf J^{(0)}=2\hat{\mathbf n}\times\mathbf H_{inc}\)
- 前次数の表面電流が作る \(\mathbf E,\mathbf H\) と \(\mathbf J^{(r+1)}=2\hat{\mathbf n}\times\mathbf H^{(r)}\)
- パネル中心、面積、代表長、壁面・観測点除外の考え方
- OpenMPによるパネル対ループの並列化

依頼に合わせ、閉じた放物面キャップは除去し、基準45°面を厳密に8回回転複製した開放八角形テーパー＋直管へ変更しました。

## 時間調和規約の修正

指定規約は \(e^{-j\omega t}\) です。したがって、前進波と外向きGreen関数を

\[
e^{+jkz},\qquad G(R)=\frac{e^{+jkR}}{4\pi R}
\]

としました。固定コミットの `mod_incident.f90` と `mod_field_integrals.f90` は `exp(-I_C*phase)`、`exp(-I_C*k0*R)` なので、実質的に逆の時間規約です。本比較版では、Green関数だけでなく近接項の複素符号も一貫して共役側へ変更しています。

sourceパネル電流を一定とした磁界核は

\[
\mathbf H_t=rac{A_s}{4\pi}e^{jkR}
\left(\frac{jk}{R}-\frac{1}{R^2}\right)
\hat{\mathbf R}\times\mathbf J_s
\]

です。これをtarget内向き法線に対して \(2\hat{\mathbf n}_t\times\mathbf H_t\) とし、3×3のベクトル伝搬写像を作ります。

## 自己項、可視判定、線形性

自己相互作用はPO境界条件に含まれるものとして伝搬核から除外します。距離softeningや代表長による近接パネル削除は行いません。隣接・近接パネルは重心則で評価するため、メッシュ細分化への感度確認が必要です。観測電場だけはPEC壁面直近を共通マスクで除外します。

sourceからtargetへの線分をz依存正八角形の解析的半空間式でサンプリングし、内部にある場合だけ採用します。さらに、光線がsource内面から出てtarget内面へ入る固定の表裏判定を使います。基準セクタで用いる判定は回転複製形状から得られるためC8対称です。

固定コミットのIPOは、各次数で合成場のPoyntingベクトルを計算し、targetごとに電流をゼロ化します。この判定は入力電流に依存する非線形演算なので、モードを独立に更新できません。厳密なmodal分離のため、本比較では固定された線形ペア判定へ変更しました。元の動的マスクを保持する場合は、各次数で「逆DFT→物理セクタでマスク→DFT」を行う必要があり、モード間結合が生じます。

## C8定式化

電流3成分はセクタ角 \(\phi_p=2\pi p/8\) だけ逆回転した局所座標で保存します。testingセクタ \(p\)、sourceセクタ \(q\) の伝搬ブロックは

\[
\mathbf K_{pq}=\mathbf K_{(q-p)\bmod8}
\]

です。`mod_operator.f90` はfull行列を全パネル対から独立に構築し、同時にtestingセクタ0と8個のsourceセクタから \(K_d\) を別経路で構築します。

DFT規約は

\[
\widetilde{\mathbf J}_m=\sum_{p=0}^{7}\mathbf J_p e^{-j2\pi mp/8},\qquad
\mathbf J_p=\frac18\sum_{m=0}^{7}\widetilde{\mathbf J}_m e^{+j2\pi mp/8},
\]

\[
\widetilde{\mathbf K}_m=\sum_{d=0}^{7}\mathbf K_d e^{+j2\pi md/8}
\]

です。全8モード版と \(m=1,7\) 限定版を実行します。

x偏波をセクタ局所座標へ戻すと \(\cos\phi_p\)、\(-\sin\phi_p\) が現れます。これらは \(e^{+j\phi_p}\) と \(e^{-j\phi_p}\) の和・差なので、軸上C8対称包絡の励振は \(m=1,7\) だけになります。各セクタ電流を同一とは仮定していません。

## 軽量検証設定

| 項目 | 値 |
|---|---:|
| 周波数 | 3 GHz |
| 入口／出口直径（対向頂点） | 180／65 mm |
| テーパー／直管長 | 160／100 mm |
| 1セクタの三角形 | 60 |
| 全パネル | 480 |
| ベクトル自由度 | 1440 |
| 最大反射次数 | 4 |
| Gaussian waist／位置 | 80 mm／入口上流250 mm |

これはfull–modalの代数検証用です。3 GHzでの物理的メッシュ収束を保証する設定ではありません。

## 実行済み参照結果

この作業環境にはFortranコンパイラが無いため、FortranソースはF2008パーサで全ファイルの構文を確認し、同じ形状・核・DFT規約を持つPython検証核で軽量ケースを実行しました。数値結果は `results/python_reference_validation` に保存されています。Fortran実行時の合格基準は同じです。

| 比較量 | 結果 | 合格基準 |
|---|---:|---:|
| block-circulant誤差 | \(9.63\times10^{-16}\) | \(<10^{-10}\) |
| 累積電流 full–全モード誤差 | \(6.34\times10^{-16}\) | \(<10^{-10}\) |
| 累積電流 full–\(m=1,7\)誤差 | \(6.40\times10^{-16}\) | \(<10^{-10}\) |
| XZ電場誤差 | \(7.75\times10^{-16}\) | \(<10^{-10}\) |
| 軸上電場誤差 | \(4.47\times10^{-16}\) | \(<10^{-10}\) |
| XY電場誤差最大 | \(6.53\times10^{-16}\) | \(<10^{-10}\) |
| 出口電場誤差 | \(4.88\times10^{-16}\) | \(<10^{-10}\) |

初期モードノルムは \(m=1:0.0738545\)、\(m=7:0.0738532\)、他モードは \(10^{-17}\) 以下でした。保存量はfull 2,073,600複素要素（31.64 MiB）、C8ブロック259,200要素（3.96 MiB）で8分の1です。Python参照での構築時間はfull 0.622 s、C8ブロック0.131 sでした。

参照実行で反射電流は減少しましたが停止条件 \(10^{-4}\) には未到達したため、最大4次で打ち切られた比較値です。更新後のFortranコードは \(\|J^{(r)}\|/\|J_{cum}^{(r)}\|\) を停止指標にします。full–C8同値性には合格していますが、IPO物理収束完了とは解釈しません。

## 出力

- `order_comparison.csv`: 各反射次数の電流比、電流誤差、更新時間
- `mode_norms.csv`: \(m=0,\ldots,7\) の初期・最終ノルム
- `surface_current.csv`: full/modal累積表面電流
- `xz_*.csv`, `xy*_*.csv`, `axis_*.csv`, `exit_*.csv`: 複素 \(E_x,E_y,E_z\)、\(|E|\)、位相、共通マスク
- `exit_order_*`, `exit_cumulative_*`, `exit_power_by_order.csv`: 反射次数別・累積出口電場と電力
- `metrics.txt`: 電流・電場誤差、出口電力、構築時間、保存量
- `plot_compare.py`: 共通カラーバーのXZ・XY比較図、各XY面の相対誤差図、`y=0`横断プロファイル

## 94 GHz実寸見積もり

入口180 mm、出口65 mm、テーパー160 mm、直管500 mmを \(\lambda/6\)（約0.532 mm）で三角形化すると、概算は次です。

| 項目 | full | C8全モード | \(m=1,7\) |
|---|---:|---:|---:|
| パネル対相当／反射次数 | \(6.6845\times10^{12}\) | \(8.3556\times10^{11}\) | \(2.0889\times10^{11}\) |
| dense複素保存量 | 896,461 GiB | 112,058 GiB | 28,014 GiB |
| 理想速度比 | 1 | 8 | 32 |

概算パネル数は1セクタ323,180、全体2,585,440です。C8だけではdense保存不能です。

適用箇所は `K_d @ J_m` です。

- GPU：Green関数、可視判定、パネル対MVPをバッチ化する。
- ACA/H-matrix：各 \(K_d\) の遠方ブロックを低ランク化し、近接だけdenseにする。
- FMM/MLFMA：\(K_d\) を保存せず、モード別MVPを階層計算する。
- \(m=1,7\) 専用：8個の \(K_d\) を保存せず、生成時に \(\widetilde K_1,\widetilde K_7\) またはそのmatrix-free作用へ直接蓄積する。

IPO次数削減は、電流ノルムと出口電力が収束した場合だけ有効です。未収束で次数を切ると高速化ではなくモデル誤差になります。
