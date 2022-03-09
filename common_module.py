import pandas as pd


def console_wide_view():
    desired_width = 320
    pd.set_option('display.width', desired_width)
    # np.set_printoptions(linewidth=desired_width)
    pd.set_option('display.max_columns', 10)