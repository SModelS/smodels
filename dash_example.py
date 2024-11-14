#!/usr/bin/env python3

"""
.. module:: runSModelS
   :synopsis: Main code for running SModelS.

"""

from smodels.matching.modelTester import getParameters,loadDatabase,loadDatabaseResults,testPoint
from dash import Dash, dcc, html, Input, Output, callback, dash_table
import numpy as np


def setLayout(app,theoryPredictions,fname):
    
    tp_dict = {}
    for tp in theoryPredictions:
        tp_dict.setdefault(tp.analysisId(),[tp.getRValue(),tp.getRValue(expected=True)])
    tp = theoryPredictions[0]
    tp_table = [{"ID" : tp.analysisId(), 
                "DataSet" : tp.dataset.getID(),
                "r(obs)" : tp.getRValue(),
                "r(exp)" : tp.getRValue(expected=True),
                "TxNames" : str(tp.txnames),
               }]

    app.layout = html.Div([
        html.Div(
            children=[
            html.H1(children="SModelS Result"),
            html.P(
                children=(
                    f"Summary of SModelS results for {fname}"
                ),
            ),
            dcc.Graph(
                id='rvaluesGraph',
                figure={
                    "data": [
                        {
                            "x": list(tp_dict.keys()),
                            "y": np.array(list(tp_dict.values()))[:,0],
                            "type": "lines",
                            "name" : "r(obs)",
                        },
                        {
                            "x": list(tp_dict.keys()),
                            "y": np.array(list(tp_dict.values()))[:,1],
                            "type": "lines",
                            "line" : {"dash": "dash", "color" : "red"},
                            "name" : "r(exp)",
                        },
                    ],
                    "layout": {"title": "r-values"},
                },
            ),        
            ],
        style={'width': '49%', 'display': 'inline-block'}
        ),
        html.Div(
                children=[
                    html.P(children=(f"Info for Analysis")),    
                    dash_table.DataTable(
                        id='analysis-data',
                        columns=[{"name" : i, "id" :  i, "deletable": True, "selectable": True} for i in tp_table[0].keys()],
                        data=tp_table,
                        editable=True,
                        filter_action="native",
                        sort_action="native",
                        sort_mode="multi",
                        column_selectable="single",
                        row_selectable="multi",
                        row_deletable=True,
                        selected_columns=[],
                        selected_rows=[],
                        page_action="native",
                        page_current= 0,
                        page_size= 10,),
                    ],
        style={'width': '49%', 'display': 'inline-block'}
        ),
        
    ],
    style={'padding': '10px 25px'})


    return app


parameterFile = './parameters.ini'
inputFile = './inputFiles/slha/simplyGluino.slha'

# Get parameters
parser = getParameters(parameterFile)

# Load database and results
database = loadDatabase(parser, None)
loadDatabaseResults(parser, database)

# Run SModelS for a single file and get the printer
outputDir = './'
output = testPoint(inputFile, outputDir, parser,
                    database)
fname = list(output.keys())[0]
data = output[fname]
theoryPredictions = data.Printers['summary'].toPrint[1]
theoryPredictions = sorted(theoryPredictions, key = lambda tp: tp.getRValue(), reverse=True)        

app = Dash(__name__)
app = setLayout(app,theoryPredictions,fname)

@callback(
    Output('analysis-data', 'data'),
    Input('rvaluesGraph', 'clickData'))
def update_table(clickData):
    if clickData is None:
        tp = theoryPredictions[0]
    else:
        anaID = clickData['points'][0]['x']
        for tp in theoryPredictions:
            if tp.analysisId() == anaID:
                break

    tp_table = [{"ID" : tp.analysisId(), 
                "DataSet" : tp.dataset.getID(),
                "r(obs)" : tp.getRValue(),
                "r(exp)" : tp.getRValue(expected=True),
                "TxNames" : str(tp.txnames),
               }]

    return tp_table



if __name__ == '__main__':
    app.run(debug=True)