import matplotlib.pyplot
import seaborn as sns
import pandas as pd


def draw_first():
    sns.set_theme(style="darkgrid")
    df = pd.read_csv('A_plot_data.csv')
    sns.lineplot(data=df)

    matplotlib.pyplot.savefig("A_plot.jpeg", dpi=2000)
    matplotlib.pyplot.close()


def draw_second():
    sns.set_theme(style="darkgrid")
    df = pd.read_csv('A2_plot_data.csv')
    sns.lineplot(data=df)

    matplotlib.pyplot.savefig("A2_plot.jpeg", dpi=2000)
    matplotlib.pyplot.close()


def main():
    draw_second()
    draw_first()


if __name__ == "__main__":
    main()
