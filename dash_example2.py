#!/usr/bin/env python3

"""
.. module:: runSModelS
   :synopsis: Main code for running SModelS.

"""

from smodels.matching.modelTester import getParameters,loadDatabase,loadDatabaseResults,testPoint
from dash import Dash, dcc, html, Input, Output, callback, dash_table
import numpy as np
import json


def setLayout(app,theoryPredictions,fname):
    
    styles = {
        'pre': {
            'border': 'thin lightgrey solid',
            'overflowX': 'scroll'
        }
    }


    app.layout = html.Div([
                    dcc.Graph(
                        id='basic-interactions',
                        figure={
                                    "data": [
                                        {
                                            "x": ['A','B','C'],
                                            "y": [2.,3.,4.],
                                            "type": "lines",
                                            "name" : "r(obs)",
                                        },
                                        {
                                            "x": ['A','B','C'],
                                            "y": [1.5,1.6,1.7],
                                            "type": "lines",
                                            "line" : {"dash": "dash", "color" : "red"},
                                            "name" : "r(exp)",
                                        },
                                    ],
                                    "layout": {"title": "r-values"},
                                },
                    ),

                    html.Div(className='row', children=[
                        html.Div([
                            dcc.Markdown("""
                                **Hover Data**

                                Mouse over values in the graph.
                            """),
                            html.Pre(id='hover-data', style=styles['pre'])
                        ], className='three columns'),

                        html.Div([
                            dcc.Markdown("""
                                **Click Data**

                                Click on points in the graph.
                            """),
                            html.Pre(id='click-data', style=styles['pre']),
                        ], className='three columns'),

                        html.Div([
                            dcc.Markdown("""
                                **Selection Data**

                                Choose the lasso or rectangle tool in the graph's menu
                                bar and then select points in the graph.

                                Note that if `layout.clickmode = 'event+select'`, selection data also
                                accumulates (or un-accumulates) selected data if you hold down the shift
                                button while clicking.
                            """),
                            html.Pre(id='selected-data', style=styles['pre']),
                        ], className='three columns'),

                        html.Div([
                            dcc.Markdown("""
                                **Zoom and Relayout Data**

                                Click and drag on the graph to zoom or click on the zoom
                                buttons in the graph's menu bar.
                                Clicking on legend items will also fire
                                this event.
                            """),
                            html.Pre(id='relayout-data', style=styles['pre']),
                        ], className='three columns')
                    ])
                ])


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
app.run_server(debug=True)

@callback(
    Output('click-data', 'children'),
    Input('basic-interactions', 'clickData'))
def display_click_data(clickData):
    return json.dumps(clickData, indent=2)
