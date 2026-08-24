# J2_hexiagonical

J1_octagonical を基にした、開放端を持つ先細り PEC 正六角形コンセントレーターの Full 3-D / C6 モード PO・IPO 解析です。ディレクトリ名は依頼時の表記 J2_hexiagonical を維持していますが、形状は regular hexagonal concentrator（正六角形）です。

## 形状

入口から l_cone の区間で正六角形の外接円半径を線形に絞り、その後 l_pipe の直管へ接続します。入口・出口はともに開放で、解析対象は6枚の PEC 内壁です。

- M=6: 60°回転で一致する6セクタ
- d_in, d_out: diameter_across_vertices=.true. のとき対向頂点間直径
- n_face: 各面の横方向分割数
- n_z_cone, n_z_pipe: テーパー部・直管部の軸方向分割数
- 総三角形パネル数: 6 * 2*n_face*(n_z_cone+n_z_pipe)

common/mod_geometry.f90 は基準面（角度 -30° から +30° の辺）を作り、60°ずつ厳密に回転複製します。内外判定にも正六角形の解析的半空間式を使います。

## C6 モードと物理的意味

各セクタの電流3成分は、セクタ角だけ逆回転した局所座標で保存されます。セクタ番号を p=0,...,5 とすると、

[
widehat{mathbf J}_m=sum_{p=0}^{5}mathbf J_p e^{-i2pi mp/6},qquad
mathbf J_p=rac{1}{6}sum_{m=0}^{5}widehat{mathbf J}_m e^{+i2pi mp/6}.
]

純粋なモード m は、隣接面へ移るたび局所電流の位相が 360*m/6 度進む離散角方向調和波です。

| m | 符号付き角次数 | 隣接面の位相差 | 物理的意味 |
|---:|---:|---:|---|
| 0 | 0 | 0° | 6面が同位相の C6 共通モード。回転に対して完全に共変 |
| 1 | +1 | +60° | 正方向一次角モード。m=5 と共役対を作る |
| 2 | +2 | +120° | 正方向二次角モード。四極子的な角変化を持つ |
| 3 | 3 | 180° | 隣接面で符号が反転する交互モード。自己共役 |
| 4 | -2 | -120° | m=2 の逆回転（共役）モード |
| 5 | -1 | -60° | m=1 の逆回転（共役）モード |

中心軸上から入る一様な横直線偏波を局所座標へ戻すと cos(phi_p) と -sin(phi_p) が現れるため、理想的な正入射平面波（または伝搬方向を z に固定した paraxial Gaussian）では一次共役対 m=1,5 だけが励起されます。本実装の Gaussian は有限な波面曲率から局所伝搬方向を計算し、その方向へ偏波を再直交化するため、小さな m=3 成分が現れることがあります。中心対称なら m=0,2,4 は丸め誤差レベルです。ビームの偏心、斜入射、形状誤差があればさらに他モードも励起され得ます。

## モード別に誘起される電場

write_mode_fields=.true. のとき、C6 IPO で得た累積表面電流を6モードへ再分解し、各モードだけを逆DFTして出口直前の観測面で電場積分を行います。

- 入射場は加えず、各表面電流モードが誘起する散乱電場だけを出力
- mode_m0_exit.csv ... mode_m5_exit.csv: 複素 Ex,Ey,Ez、|E|、マスク
- mode_field_summary.csv: 符号付き角次数、面間位相差、累積DFT電流ノルム、最大誘起電場
- 未励起モードのCSVはゼロ場となり、選択則をそのまま確認可能

plot_modes.py は次の2図を作ります。

- c6_mode_physical_meaning.png: 6面間の局所電流位相と各モードの意味
- c6_mode_induced_exit_fields.png: モード別誘起電場強度（全モード共通基準のdB）と位相を合わせた横電場方向

## ビルドと実行

    cd J2_hexiagonical
    make
    make run-modal
    make plot-modes

軽量な Full / C6 同値性比較を一式実行する場合は次のとおりです。

    make run
    make plot-all

主なターゲット:

- make run-modal: C6 モード解析とモード別電場CSV
- make run: Full、C6 modal、保存電流の後処理比較
- make plot-modes: モードの物理図と誘起電場図
- make plot: Full / C6 比較図
- make preflight-94: 94 GHz 設定の規模見積り（大規模計算は実行しない）

94 GHz 設定ではモード別の全観測点積分が非常に高価なため、既定で write_mode_fields=.false. にしています。必要な場合だけ有効化し、メッシュ・観測点数と計算資源を確認してください。

## 実装上の注意

- フェーザ規約、パネル対核、可視性判定、IPO 停止比は J1_octagonical と同じです。
- C6 モード演算子は相対セクタだけを使う行列フリー実装です。
- mode_policy=0 は全6モード、1 は一次対 m=1,5、2 は初期ノルムから自動選択します。
- 可視化に必要な Python パッケージは requirements.txt に記載されています。
