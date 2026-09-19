# figstyle.py —— 全部图的统一字体/主题（Times New Roman优先，缺失时用公制兼容的Liberation Serif）
# 把 times.ttf/timesbd.ttf/timesi.ttf/timesbi.ttf 放入 温州市科技局结题/fonts/ 即自动启用正版TNR
import os, glob
import matplotlib
from matplotlib import font_manager
from plotnine import theme_classic, theme, element_text

def setup_font(base):
    fdir = os.path.join(os.path.dirname(base), "fonts")
    tnr = False
    for p in ([] if not os.path.isdir(fdir) else glob.glob(os.path.join(fdir, "*.ttf"))):
        try:
            font_manager.fontManager.addfont(p); tnr = True
        except Exception: pass
    # 用户也可能直接把字体放在根目录
    for p in glob.glob(os.path.join(os.path.dirname(base), "times*.ttf")):
        try:
            font_manager.fontManager.addfont(p); tnr = True
        except Exception: pass
    names = {f.name for f in font_manager.fontManager.ttflist}
    fam = "Times New Roman" if ("Times New Roman" in names and tnr) or "Times New Roman" in names else "Liberation Serif"
    matplotlib.rcParams["font.family"] = fam
    matplotlib.rcParams["pdf.fonttype"] = 42
    return fam

def TH(fam, base_size=13):
    return (theme_classic(base_size=base_size) + theme(
        text=element_text(family=fam),
        axis_text=element_text(color="black", size=10, family=fam),
        axis_title=element_text(size=12, weight="bold", family=fam),
        plot_title=element_text(size=12, weight="bold", ha="center", family=fam),
        legend_text=element_text(size=9, family=fam),
        legend_title=element_text(size=10, family=fam)))

NPG = {"red":"#E64B35","blue":"#4DBBD5","green":"#00A087","navy":"#3C5488",
       "salmon":"#F39B7F","purple":"#8491B4","teal":"#91D1C2","darkred":"#DC0000","brown":"#7E6148"}
