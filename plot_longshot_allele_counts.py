import pandas as pd
import plotly.graph_objects as go
import sys
import plotly
from cyvcf2 import VCF


def main():
    df = pd.DataFrame(data=[(v.start, v.FILTER, *v.INFO["AC"])
                            for v in VCF(sys.argv[1])
                            if v.INFO["AQ"] > 15],
                      columns=['begin', 'filter', 'ref_count', 'alt_count'])
    df["ratio"] = df['alt_count'] / (df['alt_count'] + df['ref_count'])
    colordict = {'dn': 'red', 'dp': 'orange', 'sb': 'yellow', None: 'blue'}
    df["colors"] = df['filter'].apply(lambda x: colordict[x])
    fig = go.Figure(
        [go.Scatter(x=df["begin"],
                    y=df["ratio"],
                    mode='markers',
                    marker=dict(size=3, color=df['colors']),
                    name='alt-allele-ratio'),
         go.Scatter(x=[85e6, 88e6], y=[0.5] * 2, mode='lines',
                    name="heterozygous", line=dict(dash='dot')),
         go.Scatter(x=[85e6, 88e6], y=[0.66] * 2, mode='lines',
                    name="2/3", line=dict(dash='dot')),
         go.Scatter(x=[85e6, 88e6], y=[0.33] * 2, mode='lines',
                    name="1/3", line=dict(dash='dot')),
         go.Scatter(x=[85e6, 88e6], y=[1] * 2, mode='lines',
                    name="homozygous", line=dict(dash='dot')),
         ]
    )
    fig["layout"]["xaxis"].update(range=[85.5e6, 88e6])
    with open("allelic-ratio-plot.html", 'w') as output:
        output.write(plotly.offline.plot(fig, output_type="div",
                                         show_link=False, include_plotlyjs='cdn'))


if __name__ == '__main__':
    main()
