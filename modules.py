# coding: utf-8
def makedir(dirname):
    """ディレクトリ作成関数"""
    import os
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    del(os)

def makeCmap(c_range={"red": 30, "pink": 5, "white": 10, "light_green": 5, "green": 13, "blue": 17},
            c_num = {"red": "#ff0000", "pink": "#ffc0cb", "white": "#ffffff", "light_green": "#90ee90", "green": "#008000", "blue": "#0000ff"},
            c_list = ["red", "pink", "white", "light_green", "green", "blue"]): # eをかませるためのカラーマップを作る関数
    """自分で定義したカラーマップを返す(固定)"""

    import matplotlib

    c_result = []
    for color in c_list:
        for i in range(c_range[color]):
            c_result.append(c_num[color])
    cmap = matplotlib.colors.ListedColormap(c_result)
    del(matplotlib)
    return cmap

def generate_cmap(colors):
    """自分で定義したカラーマップを返す(線形補完)"""
    from matplotlib.colors import LinearSegmentedColormap
    values = range(len(colors))

    vmax = np.ceil(np.max(values))
    color_list = []
    for v, c in zip(values, colors):
        color_list.append( ( v/ vmax, c) )
    return LinearSegmentedColormap.from_list('custom_cmap', color_list)


def epsilon(x, y, val):
    """
    Nature genetics 2006のepsilonを計算する関数
    e = (Wxy - WxWy)/|min(Wx, Wy) - WxWy|
    """
    result = (val - x * y) / abs(min(x, y) - x * y)
    return result
