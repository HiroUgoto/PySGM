# PySGM

Python code to parse and analyze Strong Ground Motion records

強震動データを読み込んで解析するPythonツール群


## Requirements

* Python3

以下のパッケージが必要です
* NumPy
* SciPy
* Matplotlib

## Usage

K-NET強震動データを読み込んで，最大加速度（PGA）を求めてみよう

```py
import PySGM

file_name = "~~~.EW"  # .EWファイルのみ指定すれば良い
acc = PySGM.parse(file_name,fmt="nied")   # 自動的に .NS .UD ファイルも読み込んで，3成分記録をvectorオブジェクト（独自）に変換

acc.trend_removal()  # vectorオブジェクトでは様々な解析ができます．まずは基線補正．
acc.peak_ground_3d() # 3成分合成の最大加速度が出力されます

acc.output("~~~~~~~.acc")  # ファイルに出力できます
```

基本的な使い方は，以下の手順です
1. データ形式を指定して読み込む
2. 解析する

### Data Format

PySGM.parseで指定可能なformatです．増える可能性もあります．

* fmt="vector" (default)

PySGMで出力したテキストデータを再読み込みする

* fmt="nied"

防災科学技術研究所のK-NET, KiK-netデータ（ASCII形式）を読み込む．拡張子（.EW，.EW1，.EW2）で自動的に判断するので，K-NET，KiK-net地表／地中の区別は不要

* fmt="jma"

気象庁の強震記録を読み込む．拡張子が.csvであれば，気象庁HPで公開されているデータを読み込める．その他の拡張子であれば，DVD記録を読みこめる．



### Wave Analysis

読み込み後に利用可能な解析ツール群．以下のようにvectorオブジェクトを生成したとする．

```py
vector = PySGM.parse(file_name)
```

* 基線補正
```py
vector.trend_removal()
```

* 最大値
```py
vector.peak_ground_2d()  # 水平2成分合成値
vector.peak_ground_3d()  # 3成分合成値
PGA = vector.peak_ground_3d(print_result=False)  # コンソールに結果を出力せずに，値を変数に代入することもできる
```

* 気象庁計測震度
```py
vector.jma_seismic_intensity()
JSI = vector.jma_seismic_intensity(print_result=False)  # コンソールに結果を出力せずに，値を変数に代入することもできる
```

* 積分（加速度->速度）
```py
vel = vector.integration()  # 数値積分された結果が得られる．
vel = vector.integration(low=0.2,high=10)  # デフォルトは0.2-50Hzの帯域通過フィルタ済みだが，変更することもできる
PGV = vel.peak_ground_2d(print_result=False)  # 組み合わせれば，PGVも求まります
```

* 波形の描画
```py
velocity.plot_all()   # 3成分の時刻歴波形が表示される
velocity.plot_all(to_end=True)   # defaultでは0-60秒間のみの表示．記録全体を表示するにはto_end=Trueとする
```

* テキストファイルの出力
```py
velocity.output(output_file_name)  # 指定したファイルにテキスト形式で出力される，
```

* 応答スペクトルの計算
```py
velocity.response_spectrum()   # 成分毎に応答スペクトル（h=5%）が計算される
velocity.plot_rs_all()  # 応答スペクトルが表示される
velocity.output_rs(output_rs_file_name) # 応答スペクトルが出力される
```

この他にも色々ありますが，省略します
