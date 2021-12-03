import pandas as pd
from upsetplotly import UpSetPlotly
import plotly.express as px

def main():
    df = pd.read_csv('counts_matrix.tsv', sep="\t").set_index('ids')
    datasets = df.columns.to_list()
    annot = pd.read_csv("d4905.flair.isoforms.sqantiannot.txt", sep="\t", usecols=["isoform", "gene", "structural_category"])
    annot["identifier"] = annot["isoform"].str.cat(annot['gene'], sep="_")
    df = df.join(annot.set_index("identifier"))
    
    with open("comparison_plots.html", 'w') as output:
        fig = px.scatter(df,
                         x="d4905_dircDNArepeat_batch1",
                         y="d4905_pcrcDNA_batch1",
                         color="structural_category",
                         title="direct-cDNA vs PCR-cDNA",
                         hover_data=["gene"])
        html = fig.to_html(full_html=True, include_plotlyjs='cdn')
        output.write(html)

        for min_count in [0,5,10, 20, 50]:
            usp = UpSetPlotly(samples=[get_transcripts(df, c, min_count) for c in datasets],
                              sample_names=datasets
                              )
            html = usp.plot(order_by='decreasing', show_fig=False, return_fig=True) \
                 .update_layout(title=f"Number of transcripts detected by more than {min_count} supporting reads") \
                 .to_html(full_html=True, include_plotlyjs='cdn')
            output.write(html)

        df2 = df[df["structural_category"] != "intergenic"]
        for min_count in [0,5,10, 20, 50]:
            usp = UpSetPlotly(samples=[get_transcripts(df2, c, min_count) for c in datasets],
                              sample_names=datasets
                              )
            html = usp.plot(order_by='decreasing', show_fig=False, return_fig=True) \
                 .update_layout(title=f"Number of non-intergenic transcripts detected by more than {min_count} supporting reads") \
                 .to_html(full_html=True, include_plotlyjs='cdn')
            output.write(html)


def get_transcripts(df, column, min_count):
    return df[column][df[column] > min_count].index.to_list()

if __name__ == '__main__':
    main()
