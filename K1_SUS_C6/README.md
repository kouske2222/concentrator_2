# K1_SUS_C6 — 170 GHz SUS C6 IPO

`J2_hexiagonical` を基にした新しい独立ディレクトリです。既存ディレクトリは変更しません。
Fortran/OpenMP が1次数の計算を担当し、Python が履歴・収束・保存・再開を管理します。
本番周波数・スモークテストとも **170 GHz 専用**です。

**検証状態:** Python の保存・再開管理の単体テストは実施済み。
こちらの環境には Fortran コンパイラがないため、Fortran のコンパイル、C6/full 比較、
物理的精度の確認は未実施です。まず WSL で `make test` を実行してください。
現時点では研究用実装であり、検証済みの電磁界計算結果ではありません。

## WSL で開始

このディレクトリに移動後:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
make
OMP_NUM_THREADS=2 make test
make preflight
```

`gfortran` と `make` がなければ、WSL 側で `sudo apt install gfortran make python3-venv` を実行してください。
`make test` は、TE/TM 平面波応答、C6 と full の1次数比較、連続実行と中断再開の一致、PEC 極限を検査します。
粗い `config_smoke.nml` はソフトウェア検査用で、物理的な結果として使わないでください。

```bash
# 小規模テスト（0次から2次まで）
OMP_NUM_THREADS=2 python run.py --config config_smoke.nml --output results/smoke --max-order 2

# 本計算。まず preflight の相互作用数・ディスク量を確認してください。
OMP_NUM_THREADS=8 python run.py --config config_170ghz.nml --output results/170ghz --max-order 100
```

**再開は同じコマンドを実行するだけです。** 完了した最後の次数の次から計算します。
`--max-order` は増やせます。上限到達は収束扱いにしません。
`--steps 3` は今回新たに3次数を保存して一時停止する指定です（初回は0次も数えます）。

別ターミナルから `touch results/170ghz/STOP` とすると、現在の次数の保存後に停止します。
再開前に `mv results/170ghz/STOP results/170ghz/STOP.used` で停止要求を解除してください。
Ctrl+C、SIGTERM、プロセスの強制終了では、未完了の次数だけを再計算します。
メモリ内の途中の二重ループからの再開ではありません。

## 条件と、確認が必要な仮定

|項目|初期値・出所|
|---|---|
|周波数|170 GHz（他周波数は拒否）|
|形状|既存 J2 の六角テーパ + 六角直管、開口端、端面キャップなし|
|寸法|対頂点径 180 → 65 mm、テーパ160 mm、直管200 mm（既存 J2 の設定を継承）|
|入射|軸上 +z 方向、x偏波 Gaussian、振幅1 V/m|
|ウエスト|170 GHz の過去の値を参考に20.4 mm。位置は既存 J2 の −2.2 m を仮置き|
|材料|室温 SUS304 相当の非磁性近似、σ=1.37×10^6 S/m、μr=1|
|メッシュ|各面180分割、テーパ340分割、管400分割。精度確定値ではない|

**実験の170 GHz光学系の正確なウエスト位置・形状が確認できた場合は設定を合わせてください。**
SUS の鋼種、温度、冷間加工による磁性、表面粗さ、酸化膜は未同定です。
σ と μr は `/material/` で変更できます。`pec=.true.` は検証用 PEC 極限です。
170 GHz での実測複素導電率ではなく、室温抵抗率を周波数非依存とした近似です。

## 物理モデルとコード対応

時間依存は `exp(-iωt)`、伝搬は `exp(+ikR)`。法線 n は金属から空洞内向きです。
SI 単位で J は A/m、等価磁流 M は V/m とします。

```text
Zs = (1-i) sqrt(pi f mu0 mur / sigma)
J = n × H_total
E_t = Zs J
M = -n × E_total = -Zs n × J
```

局所平面波として c=cos(入射角)、s を入射面に垂直な接線、t を入射面内の接線とすると、
PEC 電流 `Jpec=2 n×Hinc` を次で置換します。

```text
J_s = eta/(eta+Zs*c) * Jpec_s
J_t = eta*c/(eta*c+Zs) * Jpec_t
```

初回も相互作用先でも同じ処理を適用します。次数番号に応じた任意の減衰率は掛けません。
発生源の J と M の双方を放射積分に含めます。
電流だけを小さくして PEC の電流放射をそのまま使うモデルではありません。

|計算要素|コード|
|---|---|
|境界応答・等価磁流・電磁界核|`common/mod_material.f90`|
|初回 Gaussian 入射による電流|`common/mod_incident.f90`|
|パネル間作用と電場観測|`common/mod_operator.f90`|
|6セクタの回転、三角形重心・面積|`common/mod_geometry.f90`|
|DFT、モード別作用、逆DFT|`c6_modal/mod_c6_modal.f90`|
|1次数の計算とバイナリ受渡し|`worker.f90`|
|収束とトランザクション保存|`run.py`|

C6 はベクトルをセクタ局所座標に回転してから DFT します。
全6モードを保持し、励振に対して m=±1 だけと決めつけません。
行列は構築しないため RAM は次数に比例しませんが、作用の計算時間は O(6Q²) です。
メッシュ Q=266,400/セクタでは約4.26×10^11ペア/次数となり、依然として非常に重い計算です。
ACA/MLFMA の導入や演算子キャッシュはこの版にはありません。

## 収束条件

0次を最初の PO 壁電流とし、J^(n)=K J^(n-1)、S_n=Σ_(r=0)^n J^(r) とします。
電流ノルムは面積重み付き `||J||A = sqrt(Σ_panel area*|J|²)` です。
観測電場は複素振幅のまま加算し、`E_n = E_inc + Σ E[J^(r),M^(r)]` とします。

次の **3条件すべてを5次数連続**で満たし、かつ n≥10 のとき打ち切ります。

1. `||J^(n)||A / ||S_n||A < 1e-4`
2. `||J^(n)||A / ||J^(0)||A < 1e-4`（累積電流の増大による見かけの収束も抑制）
3. `||ΔE_n||2 / ||E_n||2 < 1e-3`

分母のキャンセル対策として初期ノルムの10^-12を下限にします。
電場観測は21軸方向位置 ×（軸上1点 + 局所内接半径70%の周上8点）=189点です。
`--tol-j`, `--tol-e`, `--consecutive`, `--min-order` で変更できます。
閾値は**本実装で提案する工学的な打切り基準**であり、参考文献で保証された万能値ではありません。
連続した小さい増分でも高Q構造の遅い尾部を完全には保証しません。
より厳しい閾値・追加次数・観測点密度・メッシュ細分化で観測量の安定性を別に確認してください。
`CONVERGED` は全空間精度、Maxwell 方程式の残差、MoM との一致を意味しません。

## 保存ファイルと安全性

`results/170ghz/order_000010/` の例:

- `state.npz`: 次数電流 `j`、累積 `jsum`、次数電場 `e`、入射込み累積 `esum`、
  入射電場 `ei`、パネル面積 `areas`、観測座標 `xyz`、初期電流ノルム `jref`。
- `meta.json`: 次数、SHA-256、条件・コードの指紋、0次からの収束比・計算時間履歴。

配列は Fortran の連続データです。Python では `j.reshape(6,Q,3)`、
`e.reshape(189,3)`、`xyz.reshape(189,3)` とするとセクタ・パネル・ベクトル成分の順になります。
電流は各セクタの局所座標です。全体座標へは z 軸まわりにセクタ角を回転してください。
`esum` は189点での電場です。全面2D電場マップの自動出力は含みません。
累積電流を保存するので、後から追加の観測点・断面へ放射積分を適用できます。

```python
import numpy as np
d = np.load('results/170ghz/order_000010/state.npz')
xyz = d['xyz'].reshape(-1, 3)
E = d['esum'].reshape(-1, 3)
Eabs = np.linalg.norm(E, axis=1)
```

保存は一時ディレクトリへ書き込み、ファイルとディレクトリを fsync した後に
同一ファイルシステム上で atomic rename します。最後の保存が破損していれば、
そのディレクトリを `.corrupt-...` に移して保全し、直前の正常な保存から再開します。
`.pending-*` は未完了なので無視します。容量不足時にも直前の保存を上書きしません。
`flock` で同一出力先への二重起動を拒否します。
WSL の Linux ファイルシステム内への保存を推奨します。ディスク自体の故障の保証やバックアップの代替ではありません。

設定ファイル、実行ファイル、Fortran ソース、管理コード、収束基準が変わると再開を拒否します。
変更時は別の `--output` を指定してください。`--max-order`、`--steps`、スレッド数は変更可能です。
スレッド数・コンパイラ最適化の違いによる丸め誤差に注意してください。
本番メッシュでは保存量は約0.143 GiB/次数、0〜100次で約14.4 GiBです。

## 適用限界と残る検証

局所平面波・表面インピーダンス PO の近似であり、SUS の厳密な MoM ではありません。
相互作用方向はパネル中心間方向です。自己項除外、パネル重心一点積分、
幾何学的可視性と表裏判定は既存 J2 を継承します。
角部回折、特異・近接積分の高精度補正、粗面損失、プラズマは含みません。
界面近接場の一意な入射方向という仮定にも限界があります。
表面インピーダンスは厚さ・曲率半径に対して表皮深さが十分小さい良導体を想定します。
全体の電力収支、メッシュ収束、実験・MoM との比較は**未検証**です。
反射次数の収束だけでこれらの検証を代替しないでください。

## 参照

- 元コード: [J2_hexiagonical](https://github.com/kouske2222/concentrator_2/tree/bdc32f79a7e1ef57cd49f8932c9f8c8a08d132b3/J2_hexiagonical)
- 室温304の抵抗率0.73 Ω·mm²/m: [Outokumpu Core datasheet](https://www.outokumpu.com/-/media/files/products/core/outokumpu-core-range-datasheet.pdf)
- 表面インピーダンスの適用の背景: [COMSOL AC/DC introduction](https://doc.comsol.com/6.4/doc/com.comsol.help.acdc/IntroductionToACDCModule.pdf)

上記資料は材料値・境界条件の背景資料です。このコードのIPO近似や収束閾値の妥当性を直接検証するものではありません。
