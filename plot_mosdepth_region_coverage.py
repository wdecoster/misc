import pandas as pd
import plotly.graph_objects as go
import sys
import plotly


def main():
    df = pd.read_csv(sys.argv[1],
                     sep="\t",
                     header=None,
                     names=["chrom", "begin", "end", "coverage"]) \
        .rolling(3, center=True).mean() \
        .dropna()
    median = df["coverage"].median()
    fig = go.Figure(
        [go.Scatter(x=df["begin"], y=df["coverage"], mode='markers', marker=dict(size=3), name='coverage'),
         go.Scatter(x=[85e6, 88e6], y=[median] * 2, mode='lines',
                    name="median", line=dict(dash='dot')),
         go.Scatter(x=[85e6, 88e6], y=[median * 1.5] * 2, mode='lines',
                    name="duplication", line=dict(dash='dot')),
         go.Scatter(x=[85e6, 88e6], y=[median / 2] * 2, mode='lines',
                    name="deletion", line=dict(dash='dot'))
         ]

    )
    fig["layout"]["xaxis"].update(range=[85.5e6, 88e6])
    with open("region-plot.html", 'w') as output:
        output.write(plotly.offline.plot(fig, output_type="div",
                                         show_link=False, include_plotlyjs='cdn'))


if __name__ == '__main__':
    main()
