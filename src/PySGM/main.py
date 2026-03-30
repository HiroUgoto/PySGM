from . import vector
from . import nied,jma
from . import win
from . import gns
from . import sac
import os
import pickle

def parse(file_name,file_name2=None,file_name3=None,fmt="vector",
    gain=None,bit=None,noheader=False):
    """
    強震動データを読み込み、vectorsオブジェクトとしてパースします。

    fmtに指定したフォーマットに応じて適切に変換します。例えば，防災科研K-NET/KiK-netの場合、fmt="nied"と指定すれば拡張子に応じて自動的に変換される。
    
    Args:
        file_name (str): 読み込むファイルのパス（K-NET/KiK-netの場合はEWファイル名を指定すること。"~~~~.EW"，"~~~~~.EW2"など）
        file_name2 (str, optional): 追加のファイルパス（SAC形式用など）。ほとんどの場合は不要。
        file_name3 (str, optional): 追加のファイルパス（SAC形式用など）。ほとんどの場合は不要。
        fmt (str, optional): データのフォーマット ("vector", "nied", "jma", "win", "sac")
        gain (float, optional): 振幅倍率（win形式のみ）
        bit (int, optional): データビット数（win形式のみ）
        noheader (bool, optional): ヘッダーがあるかどうか（vectors形式のみ）

    Returns:
        (vectors): 読み込まれたデータの vectorsインスタンス
    """
    if fmt == "vector":
        v = vector.parse(file_name,noheader)
        return v

    elif fmt == "nied":
        file_basename,ext = os.path.splitext(file_name)
        v = nied.parse(file_basename,ext)
        return v

    elif fmt == "jma":
        file_basename,ext = os.path.splitext(file_name)
        if ext == ".csv":
            v = jma.csv_parse(file_name)
        else:
            v = jma.parse(file_name)
        return v

    elif fmt == "win":
        v = win.parse(file_name,gain,bit)
        return v

    elif fmt == "gns":
        v = gns.parse(file_name)
        return v

    elif fmt == "sac":
        v = sac.parse(file_name,file_name2,file_name3)
        return v

def input_pickle(input_file):
    with open(input_file,'rb') as f:
        v = pickle.load(f)
    return v