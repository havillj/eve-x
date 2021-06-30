import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
from lxml import etree
from matplotlib.backends.backend_pdf import PdfPages
import io
import base64
from pipeline_utilities import *
from dash.dependencies import Input, Output, State
import matplotlib
matplotlib.use('Agg')
import math
#import dash_design_kit as ddk
#import dash_daq as daq



import dash_bootstrap_components as dbc

all_viruses = ['Aedes aegypti anphevirus', 'Aedes anphevirus', 'American plum line pattern virus', 'Australian Anopheles totivirus',
 'Cell fusing agent virus', 'Culex Flavi-like virus', 'Culex ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Culiseta flavivirus',
  'Dengue virus', 'Kamiti River virus', 'Lassa mammarenavirus', 'Liao ning virus', 'Mercadeo virus', 'Merida virus',
   'Ochlerotatus caspius flavivirus', 'Ochlerotatus scapularis flavivirus', 'Ohlsdorf ohlsrhavirus', 'Pestivirus A', 'Radi vesiculovirus',
   'Tongilchon ohlsrhavirus', 'West Nile virus', 'Wuhan Mosquito Virus 6', 'Xishuangbanna aedes flavivirus']

vFam = {'Arenaviridae': ['Lassa mammarenavirus'], 'Rhabdoviridae': ['Culex ohlsrhavirus', 'Ohlsdorf ohlsrhavirus', 'Culex tritaeniorhynchus rhabdovirus', 'Tongilchon ohlsrhavirus',
        'Merida virus', 'Radi vesiculovirus'], 'Xinmoviridae': ['Aedes anphevirus', 'Aedes aegypti anphevirus'], 'Orthomyxoviridae': ['Wuhan Mosquito Virus 6'],
        'Sedoreovirinae': ['Liao ning virus'], 'Flaviviridae': ['Xishuangbanna aedes flavivirus', 'Pestivirus A', 'Culiseta flavivirus', 'Mercadeo virus', 'Ochlerotatus scapularis flavivirus',
        'Cell fusing agent virus', 'West Nile virus', 'Ochlerotatus caspius flavivirus', 'Kamiti River virus', 'Dengue virus', 'Culex Flavi-like virus'],
        'Totiviridae': ['Australian Anopheles totivirus'], 'Bromoviridae': ['American plum line pattern virus']}


all_regions = {
    'Angola': [1,2,3,4,5,6,8,9,10,11,12,14,16,17,18,19,20],
    'Gabon': [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,23,24],
    'South-Africa': [2,3,4,5,6,7,8,9,10,11,12,13,14,15],
    'Mexico': [1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24],
    'USA': [1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,18,19,20,21,22,23,24],
    'Australia': [1,2,5,7,8,10,11,12,13,15,16,17,18,19,20,21,22,23,24],
    'French-Polynesia': [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,10],
    'Philippines': [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21,22,23,24],
    'Argentina': [2,3,8,11,31],
    'Brazil': [1,2,3,4,6,7,8,10,11,12,13,14,15,18,19,23,24],
    'Thailand': [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],
    'Vietnam': [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],
}
label2Specimen = {('Mexico', 1): 'Amacuzac-Mexico-1.LIN210A1618', ('Mexico', 10): 'Amacuzac-Mexico-10.LIN210A1625', ('Mexico', 11): 'Amacuzac-Mexico-11.LIN210A1626', ('Mexico', 13): 'Amacuzac-Mexico-13.LIN210A1627', ('Mexico', 14): 'Amacuzac-Mexico-14.LIN210A1628', ('Mexico', 16): 'Amacuzac-Mexico-16.LIN210A1629', ('Mexico', 17): 'Amacuzac-Mexico-17.LIN210A1630', ('Mexico', 18): 'Amacuzac-Mexico-18.LIN210A1631', ('Mexico', 19): 'Amacuzac-Mexico-19.LIN210A1632', ('Mexico', 20): 'Amacuzac-Mexico-20.LIN210A1633', ('Mexico', 22): 'Amacuzac-Mexico-22.LIN210A1634', ('Mexico', 23): 'Amacuzac-Mexico-23.LIN210A1635', ('Mexico', 24): 'Amacuzac-Mexico-24.LIN210A1636', ('Mexico', 3): 'Amacuzac-Mexico-3.LIN210A1619', ('Mexico', 4): 'Amacuzac-Mexico-4.LIN210A1620', ('Mexico', 5): 'Amacuzac-Mexico-5.LIN210A1621', ('Mexico', 6): 'Amacuzac-Mexico-6.LIN210A1622', ('Mexico', 7): 'Amacuzac-Mexico-7.LIN210A1623', ('Mexico', 8): 'Amacuzac-Mexico-8.LIN210A1624', ('Mexico', 15): 'Amacuzac_Mexico-_15.2-116629', ('Mexico', 2): 'Amacuzac_Mexico-_2.2-116616', ('Mexico', 21): 'Amacuzac_Mexico-_21.2-116635', ('Thailand', 1): 'Bangkok_Thailand_01.LIN210A1679', ('Thailand', 2): 'Bangkok_Thailand_02.LIN210A1680', ('Thailand', 3): 'Bangkok_Thailand_03.LIN210A1681', ('Thailand', 4): 'Bangkok_Thailand_04.LIN210A1682', ('Thailand', 5): 'Bangkok_Thailand_05.LIN210A1683', ('Thailand', 6): 'Bangkok_Thailand_06.LIN210A1684', ('Thailand', 7): 'Bangkok_Thailand_07.LIN210A1685', ('Thailand', 8): 'Bangkok_Thailand_08.LIN210A1686', ('Thailand', 9): 'Bangkok_Thailand_09.LIN210A1687', ('Thailand', 10): 'Bangkok_Thailand_10.LIN210A1688', ('Thailand', 11): 'Bangkok_Thailand_11.LIN210A1689', ('Thailand', 12): 'Bangkok_Thailand_12.LIN210A1690', ('Thailand', 13): 'Bangkok_Thailand_13.LIN210A1691', ('Thailand', 14): 'Bangkok_Thailand_14.LIN210A1692', ('Thailand', 15): 'Bangkok_Thailand_15.LIN210A1693', ('Thailand', 16): 'Bangkok_Thailand_16.LIN210A1694', ('Thailand', 17): 'Bangkok_Thailand_17.LIN210A1695', ('Thailand', 18): 'Bangkok_Thailand_18.LIN210A1696', ('Thailand', 19): 'Bangkok_Thailand_19.LIN210A1697', ('Thailand', 20): 'Bangkok_Thailand_20.LIN210A1698', ('Australia', 1): 'Cairns-Australia-1.LIN210A2151', ('Australia', 10): 'Cairns-Australia-10.LIN210A2160', ('Australia', 11): 'Cairns-Australia-11.LIN210A2161', ('Australia', 12): 'Cairns-Australia-12.LIN210A2162', ('Australia', 13): 'Cairns-Australia-13.LIN210A2163', ('Australia', 15): 'Cairns-Australia-15.LIN210A2165', ('Australia', 16): 'Cairns-Australia-16.LIN210A2166', ('Australia', 17): 'Cairns-Australia-17.LIN210A2167', ('Australia', 18): 'Cairns-Australia-18.LIN210A2168', ('Australia', 19): 'Cairns-Australia-19.LIN210A2169', ('Australia', 2): 'Cairns-Australia-2.LIN210A2152', ('Australia', 20): 'Cairns-Australia-20.LIN210A2170', ('Australia', 21): 'Cairns-Australia-21.LIN210A2171', ('Australia', 22): 'Cairns-Australia-22.LIN210A2172', ('Australia', 23): 'Cairns-Australia-23.LIN210A2173', ('Australia', 24): 'Cairns-Australia-24.LIN210A2174', ('Australia', 5): 'Cairns-Australia-5.LIN210A2155', ('Australia', 7): 'Cairns-Australia-7.LIN210A2157', ('Australia', 8): 'Cairns-Australia-8.LIN210A2158', ('Philippines', 1): 'Cebu_City_Philippines-_1.2-116713', ('Philippines', 10): 'Cebu_City_Philippines-_10.2-116722', ('Philippines', 11): 'Cebu_City_Philippines-_11.2-116723', ('Philippines', 12): 'Cebu_City_Philippines-_12.2-116724', ('Philippines', 13): 'Cebu_City_Philippines-_13.2-116725', ('Philippines', 14): 'Cebu_City_Philippines-_14.2-116726', ('Philippines', 15): 'Cebu_City_Philippines-_15.2-116727', ('Philippines', 16): 'Cebu_City_Philippines-_16.2-116728', ('Philippines', 18): 'Cebu_City_Philippines-_18.2-116730', ('Philippines', 19): 'Cebu_City_Philippines-_19.2-116731', ('Philippines', 2): 'Cebu_City_Philippines-_2.2-116714', ('Philippines', 20): 'Cebu_City_Philippines-_20.2-116732', ('Philippines', 21): 'Cebu_City_Philippines-_21.2-116733', ('Philippines', 22): 'Cebu_City_Philippines-_22.2-116734', ('Philippines', 23): 'Cebu_City_Philippines-_23.2-116735', ('Philippines', 24): 'Cebu_City_Philippines-_24.2-116736', ('Philippines', 3): 'Cebu_City_Philippines-_3.2-116715', ('Philippines', 4): 'Cebu_City_Philippines-_4.2-116716', ('Philippines', 5): 'Cebu_City_Philippines-_5.2-116717', ('Philippines', 6): 'Cebu_City_Philippines-_6.2-116718', ('Philippines', 7): 'Cebu_City_Philippines-_7.2-116719', ('Philippines', 8): 'Cebu_City_Philippines-_8.2-116720', ('Philippines', 9): 'Cebu_City_Philippines-_9.2-116721', ('Angola', 1): 'Cuanda_Angola_01.LIN210A1719', ('Angola', 2): 'Cuanda_Angola_02.LIN210A1720', ('Angola', 3): 'Cuanda_Angola_03.LIN210A1721', ('Angola', 4): 'Cuanda_Angola_04.LIN210A1722', ('Angola', 5): 'Cuanda_Angola_05.LIN210A1723', ('Angola', 6): 'Cuanda_Angola_06.LIN210A1724', ('Angola', 8): 'Cuanda_Angola_08.LIN210A1726', ('Angola', 9): 'Cuanda_Angola_09.LIN210A1727', ('Angola', 10): 'Cuanda_Angola_10.LIN210A1728', ('Angola', 11): 'Cuanda_Angola_11.LIN210A1729', ('Angola', 12): 'Cuanda_Angola_12.LIN210A1730', ('Angola', 14): 'Cuanda_Angola_14.LIN210A1732', ('Angola', 16): 'Cuanda_Angola_16.LIN210A1734', ('Angola', 17): 'Cuanda_Angola_17.LIN210A1735', ('Angola', 18): 'Cuanda_Angola_18.LIN210A1736', ('Angola', 19): 'Cuanda_Angola_19.LIN210A1737', ('Angola', 20): 'Cuanda_Angola_20.LIN210A1738', ('Argentina', 11): 'El_Dorado_Argentina_U_11.LIN210A2770', ('Argentina', 2): 'El_Dorado_Argentina_U_2.LIN210A2698', ('Argentina', 3): 'El_Dorado_Argentina_U_3.LIN210A2706', ('Argentina', 8): 'El_Dorado_Argentina_U_8.LIN210A2746', ('Argentina', 31): 'El_Dorado_US_U_31', ('Vietnam', 1): 'HoChiMin_Vietnam_01.LIN210A1659', ('Vietnam', 2): 'HoChiMin_Vietnam_02.LIN210A1660', ('Vietnam', 3): 'HoChiMin_Vietnam_03.LIN210A1661', ('Vietnam', 4): 'HoChiMin_Vietnam_04.LIN210A1662', ('Vietnam', 5): 'HoChiMin_Vietnam_05.LIN210A1663', ('Vietnam', 6): 'HoChiMin_Vietnam_06.LIN210A1664', ('Vietnam', 7): 'HoChiMin_Vietnam_07.LIN210A1665', ('Vietnam', 8): 'HoChiMin_Vietnam_08.LIN210A1666', ('Vietnam', 9): 'HoChiMin_Vietnam_09.LIN210A1667', ('Vietnam', 10): 'HoChiMin_Vietnam_10.LIN210A1668', ('Vietnam', 11): 'HoChiMin_Vietnam_11.LIN210A1669', ('Vietnam', 12): 'HoChiMin_Vietnam_12.LIN210A1670', ('Vietnam', 13): 'HoChiMin_Vietnam_13.LIN210A1671', ('Vietnam', 14): 'HoChiMin_Vietnam_14.LIN210A1672', ('Vietnam', 15): 'HoChiMin_Vietnam_15.LIN210A1673', ('Vietnam', 16): 'HoChiMin_Vietnam_16.LIN210A1674', ('Vietnam', 17): 'HoChiMin_Vietnam_17.LIN210A1675', ('Vietnam', 18): 'HoChiMin_Vietnam_18.LIN210A1676', ('Vietnam', 19): 'HoChiMin_Vietnam_19.LIN210A1677', ('Vietnam', 20): 'HoChiMin_Vietnam_20.LIN210A1678', ('Gabon', 1): 'La_Lope-Gabon-1.LIN210A1637', ('Gabon', 10): 'La_Lope-Gabon-10.LIN210A1646', ('Gabon', 11): 'La_Lope-Gabon-11.LIN210A1647', ('Gabon', 12): 'La_Lope-Gabon-12.LIN210A1648', ('Gabon', 13): 'La_Lope-Gabon-13.LIN210A1649', ('Gabon', 14): 'La_Lope-Gabon-14.LIN210A1650', ('Gabon', 15): 'La_Lope-Gabon-15.LIN210A1651', ('Gabon', 16): 'La_Lope-Gabon-16.LIN210A1652', ('Gabon', 17): 'La_Lope-Gabon-17.LIN210A1653', ('Gabon', 18): 'La_Lope-Gabon-18.LIN210A1654', ('Gabon', 2): 'La_Lope-Gabon-2.LIN210A1638', ('Gabon', 20): 'La_Lope-Gabon-20.LIN210A1655', ('Gabon', 23): 'La_Lope-Gabon-23.LIN210A1656', ('Gabon', 24): 'La_Lope-Gabon-24.LIN210A1658', ('Gabon', 3): 'La_Lope-Gabon-3.LIN210A1639', ('Gabon', 4): 'La_Lope-Gabon-4.LIN210A1640', ('Gabon', 5): 'La_Lope-Gabon-5.LIN210A1641', ('Gabon', 6): 'La_Lope-Gabon-6.LIN210A1642', ('Gabon', 7): 'La_Lope-Gabon-7.LIN210A1643', ('Gabon', 8): 'La_Lope-Gabon-8.LIN210A1644', ('Gabon', 9): 'La_Lope-Gabon-9.LIN210A1645', ('Gabon', 19): 'La_Lope_Gabon-_19.2-116707', ('Gabon', 22): 'La_Lope_Gabon-_22.2-116710', ('USA', 1): 'Maricopa_County_AZ-_1.2-116541', ('USA', 12): 'Maricopa_County_AZ-_12.2-116602', ('USA', 13): 'Maricopa_County_AZ-_13.2-116603', ('USA', 14): 'Maricopa_County_AZ-_14.2-116604', ('USA', 15): 'Maricopa_County_AZ-_15.2-116605', ('USA', 16): 'Maricopa_County_AZ-_16.2-116606', ('USA', 18): 'Maricopa_County_AZ-_18.2-116608', ('USA', 19): 'Maricopa_County_AZ-_19.2-116609', ('USA', 21): 'Maricopa_County_AZ-_21.2-116611', ('USA', 24): 'Maricopa_County_AZ-_24.2-116614', ('USA', 3): 'Maricopa_County_AZ-_3.2-116543', ('USA', 4): 'Maricopa_County_AZ-_4.2-116544', ('USA', 6): 'Maricopa_County_AZ-_6.2-116546', ('USA', 8): 'Maricopa_County_AZ-_8.2-116548', ('USA', 9): 'Maricopa_County_AZ-_9.2-116549', ('USA', 10): 'Maricopa_County_Arizona-USA-10.LIN210A1614', ('USA', 2): 'Maricopa_County_Arizona-USA-2.LIN210A1611', ('USA', 20): 'Maricopa_County_Arizona-USA-20.LIN210A1615', ('USA', 22): 'Maricopa_County_Arizona-USA-22.LIN210A1616', ('USA', 23): 'Maricopa_County_Arizona-USA-23.LIN210A1617', ('USA', 5): 'Maricopa_County_Arizona-USA-5.LIN210A1612', ('USA', 7): 'Maricopa_County_Arizona-USA-7.LIN210A1613', ('Brazil', 1): 'Paqueta-Brazil-1.LIN210A2127', ('Brazil', 10): 'Paqueta-Brazil-10.LIN210A2136', ('Brazil', 11): 'Paqueta-Brazil-11.LIN210A2137', ('Brazil', 12): 'Paqueta-Brazil-12.LIN210A2138', ('Brazil', 13): 'Paqueta-Brazil-13.LIN210A2139', ('Brazil', 14): 'Paqueta-Brazil-14.LIN210A2140', ('Brazil', 15): 'Paqueta-Brazil-15.LIN210A2141', ('Brazil', 18): 'Paqueta-Brazil-18.LIN210A2144', ('Brazil', 19): 'Paqueta-Brazil-19.LIN210A2145', ('Brazil', 2): 'Paqueta-Brazil-2.LIN210A2128', ('Brazil', 23): 'Paqueta-Brazil-23.LIN210A2149', ('Brazil', 24): 'Paqueta-Brazil-24.LIN210A2150', ('Brazil', 3): 'Paqueta-Brazil-3.LIN210A2129', ('Brazil', 4): 'Paqueta-Brazil-4.LIN210A2130', ('Brazil', 6): 'Paqueta-Brazil-6.LIN210A2132', ('Brazil', 7): 'Paqueta-Brazil-7.LIN210A2133', ('Brazil', 8): 'Paqueta-Brazil-8.LIN210A2134', ('South-Africa', 2): 'Skukusa_South_Africa_02.LIN210A1740', ('South-Africa', 3): 'Skukusa_South_Africa_03.LIN210A1741', ('South-Africa', 4): 'Skukusa_South_Africa_04.LIN210A1742', ('South-Africa', 5): 'Skukusa_South_Africa_05.LIN210A1743', ('South-Africa', 6): 'Skukusa_South_Africa_06.LIN210A1744', ('South-Africa', 7): 'Skukusa_South_Africa_07.LIN210A1745', ('South-Africa', 8): 'Skukusa_South_Africa_08.LIN210A1746', ('South-Africa', 9): 'Skukusa_South_Africa_09.LIN210A1747', ('South-Africa', 10): 'Skukusa_South_Africa_10.LIN210A1748', ('South-Africa', 11): 'Skukusa_South_Africa_11.LIN210A1749', ('South-Africa', 12): 'Skukusa_South_Africa_12.LIN210A1750', ('South-Africa', 13): 'Skukusa_South_Africa_13.LIN210A1751', ('South-Africa', 14): 'Skukusa_South_Africa_14.LIN210A1752', ('South-Africa', 15): 'Skukusa_South_Africa_15.LIN210A1753', ('French-Polynesia', 1): 'Tahiti_FrenchPolynesia_01.LIN210A1699', ('French-Polynesia', 2): 'Tahiti_FrenchPolynesia_02.LIN210A1700', ('French-Polynesia', 3): 'Tahiti_FrenchPolynesia_03.LIN210A1701', ('French-Polynesia', 4): 'Tahiti_FrenchPolynesia_04.LIN210A1702', ('French-Polynesia', 5): 'Tahiti_FrenchPolynesia_05.LIN210A1703', ('French-Polynesia', 6): 'Tahiti_FrenchPolynesia_06.LIN210A1704', ('French-Polynesia', 7): 'Tahiti_FrenchPolynesia_07.LIN210A1705', ('French-Polynesia', 8): 'Tahiti_FrenchPolynesia_08.LIN210A1706', ('French-Polynesia', 9): 'Tahiti_FrenchPolynesia_09.LIN210A1707', ('French-Polynesia', 10): 'Tahiti_FrenchPolynesia_10.LIN210A1708', ('French-Polynesia', 11): 'Tahiti_FrenchPolynesia_11.LIN210A1709', ('French-Polynesia', 12): 'Tahiti_FrenchPolynesia_12.LIN210A1710', ('French-Polynesia', 13): 'Tahiti_FrenchPolynesia_13.LIN210A1711', ('French-Polynesia', 14): 'Tahiti_FrenchPolynesia_14.LIN210A1712', ('French-Polynesia', 15): 'Tahiti_FrenchPolynesia_15.LIN210A1713', ('French-Polynesia', 16): 'Tahiti_FrenchPolynesia_16.LIN210A1714', ('French-Polynesia', 17): 'Tahiti_FrenchPolynesia_17.LIN210A1715', ('French-Polynesia', 18): 'Tahiti_FrenchPolynesia_18.LIN210A1716', ('French-Polynesia', 19): 'Tahiti_FrenchPolynesia_19.LIN210A1717', ('French-Polynesia', 20): 'Tahiti_FrenchPolynesia_20.LIN210A1718'}


checklist = dbc.FormGroup(
    [
        dbc.Label("Matches or Inversions"),
        dbc.Checklist(
            options=[
                {"label": "Matches", "value": 1},
                {"label": "Inversions", "value": 2},

            ],
            value=[1],
            id="match-inv-input",
        ),
    ]
)

app = dash.Dash()
colors = {
    'background': '#FFFFFF',
    'text': '#111111'
}
app.layout = html.Div(style={'backgroundColor': colors['background']}, children=[
    html.H1(
        children='EVEs In Aedes aegypti Mosquitoes',
        style={
            'textAlign': 'center',
            'color': colors['text']
        }
    ),
    html.Div(style={
        'textAlign': 'center',
        'color': colors['text']},
        children=['Contig Diagrams for EVE in AAa genome',
        #html.Figure(figure1)
        ]
    ),

    html.P("To enable more convenient appraisal of individual contigs, EVE produces diagrams of the viral hits in each contig (only the best hit). To see different contigs containing viral hits, select a specimen and virus.",
     style = {'textAlign': 'center'}),
    html.Br(),

    html.Div ([dcc.Dropdown( id='countries-dropdown',placeholder = 'Select a region',options=[{'label': k, 'value': k} for k in all_regions.keys()])], style = {'width':'20%', 'display': 'inline-block'},
    ),

    html.Div([dcc.Dropdown(id='numbers-dropdown', placeholder = 'Select a num'),], style = {'width':'20%', 'display': 'inline-block'},),

    html.Div([dcc.Dropdown(id='family-dropdown', options = [{'label': k, 'value': k} for k in vFam.keys()], placeholder = 'Select a family')], style = {'width': '30%', 'display': 'inline-block'}),

    html.Div([dcc.Dropdown(id='virus-dropdown', placeholder = 'Select a virus'),], style = {'width':'30%', 'display': 'inline-block'}),




    html.Br(),

    dbc.Button('submit', id = 'submit', n_clicks = 0),

    html.P(id = 'result', style = {'display':'none'}),

    html.Br(),
    html.Br(),

    html.Div(id = 'contig-div', children = [html.P(id = 'choose-contig', style = {'font-size':14}),
            dcc.Input(id='contig-num', type ='number', placeholder = 'choose a contig', min=0, max = 5, value = 1),], style = {'display':'none'}),

    html.Div (id = 'match-div', children = [html.P(id = 'choose-matches', style = {'font-size':14}), dcc.Input(id='matches-num', type ='number', placeholder = 'Number of inv', min=0,  value = 0, )],
            style = {'display':'none'}),

    html.Div (id = 'inv-div', children =[html.P(id = 'choose-inv', style = {'font-size':14}), dcc.Input(id='inv-num', type ='number', placeholder = 'Number of inv', min=0, value = 0)],
             style = {'display':'none', "margin-left": "20px"}),

    html.Div (id = 'left-div', children =[html.P(id = 'choose-left', style = {'font-size':14}), dcc.Input(id='left-num', type ='number', placeholder = 'Number of inv', min=0, value = 0)],
             style = {'display':'none', "margin-left": "20px"}),

    html.Div (id = 'right-div', children =[html.P(id = 'choose-right', style = {'font-size':14}), dcc.Input(id='right-num', type ='number', placeholder = 'Number of inv', min=0, value = 0)],
             style = {'display':'none', "margin-left": "20px"}),

    html.Hr(),


    html.Div([
          dbc.Button('prev', id='left-scroll',n_clicks=0),
          dbc.Button('next', id = 'right-scroll', n_clicks = 0),
          ]),

    html.Img(id='contigDiagram')


])


@app.callback(
    dash.dependencies.Output('choose-matches', 'children'),
    dash.dependencies.Output('match-div', 'style'),
    dash.dependencies.Input('result', 'children'),
    state = [dash.dependencies.State('matches-num', 'max')],
    prevent_initial_call = True
    )

def set_matchNum_caption(submit, matches):
    return f"Choose number of matches (max: {matches})", {'display':'inline-block', "margin-left": "20px"}

@app.callback(
    dash.dependencies.Output('choose-inv', 'children'),
    dash.dependencies.Output('inv-div', 'style'),
    dash.dependencies.Input('result', 'children'),
    state = [dash.dependencies.State('inv-num', 'max')],
    prevent_initial_call = True
    )

def set_invNum_caption(state, inv):
    return f"Choose number of inversions (max: {inv})", {'display':'inline-block', "margin-left": "20px"}

@app.callback(
    dash.dependencies.Output('choose-left', 'children'),
    dash.dependencies.Output('left-div', 'style'),
    dash.dependencies.Input('result', 'children'),
    state = [dash.dependencies.State('left-num', 'max')],
    prevent_initial_call = True
    )

def set_leftNum_caption(state, left):
    return f"Choose number of left flanks (max: {left})", {'display':'inline-block', "margin-left": "20px"}

@app.callback(
    dash.dependencies.Output('choose-right', 'children'),
    dash.dependencies.Output('right-div', 'style'),
    dash.dependencies.Input('result', 'children'),
    state = [dash.dependencies.State('right-num', 'max')],
    prevent_initial_call = True
    )

def set_rightNum_caption(state, right):
    return f"Choose number of right flanks (max: {right})", {'display':'inline-block', "margin-left": "20px"}

@app.callback(
    dash.dependencies.Output('choose-contig', 'children'),
    dash.dependencies.Output('contig-div', 'style'),
    dash.dependencies.Input('result', 'children'),
    state = [dash.dependencies.State('contig-num', 'max')],
    prevent_initial_call = True
    )

def set_contig_caption(state, contig_max):
    return f"Choose a contig number (max: {contig_max})", {'display':'inline-block'}


@app.callback(
    dash.dependencies.Output('virus-dropdown', 'options'),
    dash.dependencies.Input('family-dropdown', 'value'), prevent_initial_call = True)
def set_virus_options(selected_family):
    return [{'label': i, 'value': i} for i in vFam[selected_family]]


@app.callback(
    dash.dependencies.Output('numbers-dropdown', 'options'),
    dash.dependencies.Input('countries-dropdown', 'value'), prevent_initial_call = True)
def set_cities_options(selected_country):
    return [{'label': i, 'value': i} for i in all_regions[selected_country]]


@app.callback(
    dash.dependencies.Output('virus-dropdown', 'value'),
    dash.dependencies.Input('virus-dropdown', 'options'), prevent_initial_call = True)
def set_virus_value(available_options):
    return available_options[0]['value']

@app.callback(
    Output('contig-num', 'max'),
    Output('result','children'),
    Output('matches-num', 'max'),
    Output('inv-num', 'max'),
    Output('left-num', 'max'),
    Output('right-num', 'max'),
    [dash.dependencies.Input('submit', 'n_clicks'),
    dash.dependencies.Input('contig-num', 'value')],
    state = [dash.dependencies.State('countries-dropdown', 'value'),
    dash.dependencies.State('numbers-dropdown', 'value'),
    dash.dependencies.State('virus-dropdown', 'value'),
    #dash.dependencies.State('contig-num', 'value'),

    ], prevent_initial_call = True
    )

def pickContig(state, cNum, spec, specNum, virus):

    cNum -= 1
    dir = label2Specimen[(spec,int(specNum))]
    pathN = '/Volumes/Data2/specimens/'+ dir + '/results/xml/' + dir +'_hits.xml'

    nodeL = []

    parser = etree.XMLParser(remove_blank_text=True)
    root = etree.parse(pathN, parser).getroot()
    try:
        root.xpath(f"""/root/contig[@besthit="True"]/virushit[@stitle = '{virus}']""")[0]
    except:
        return 0, "does not exist", 0,0,0,0
    contigNames = root.xpath(f"""/root/contig[@besthit="True"][virushit[@stitle = '{virus}']]/@name""")
    #for contigName in contigNames:
    contig = root.xpath(f"""/root/contig[@besthit="True" and @name = '{contigNames[cNum]}'][virushit[@stitle = '{virus}']]""")[0]
    caption = contigNames[cNum]#f"{virus} insertion sites in {spec} {specNum} ({contigName})"
    nodeL.append(caption)
    flanks = contig.xpath('./flanks')[0]
    matches = len(flanks.xpath('./match'))
    inversions = len(flanks.findall('inversion'))
    hitsLeft = len(contig.xpath('./vectorhitleft'))
    hitsRight = len(contig.xpath('./vectorhitright'))
    #print(matches, inversions, hitsLeft, hitsRight)

    return len(contigNames), pathN, matches, inversions, hitsLeft, hitsRight

@app.callback(
    Output('contigDiagram', 'src'),
    Output('left-scroll', 'disabled'),
    Output('right-scroll', 'disabled'),
    [dash.dependencies.Input('result', 'children'),
    dash.dependencies.Input('contig-num', 'value'),
    dash.dependencies.Input('matches-num', 'value'),
    dash.dependencies.Input('inv-num', 'value'),
    dash.dependencies.Input('left-num', 'value'),
    dash.dependencies.Input('right-num', 'value'),
    dash.dependencies.Input('left-scroll', 'n_clicks'),
    dash.dependencies.Input('right-scroll', 'n_clicks'),
    ],
    state =
    [#dash.dependencies.State('result', 'children'),
    #dash.dependencies.State('image-num', 'value'),
    dash.dependencies.State('virus-dropdown', 'value'),
    dash.dependencies.State('countries-dropdown', 'value'),
    dash.dependencies.State('numbers-dropdown', 'value'),
    State('matches-num', 'max'),
    State('inv-num', 'max'),
    State('left-num', 'max'),
    State('right-num', 'max')
    ],
     prevent_initial_call = True
    )
def drawContig(path, contigNum, matchNum, invNum, leftFlanks, rightFlanks, prev, next, virus, spec, specNum, mMax, invMax, lMax, rMax):

    if path == "does not exist":
        fig = pyplot.figure()
        fig.text(.5, .5, 'No data found', ha='center')
        buf = io.BytesIO()
        pyplot.savefig(buf, format='png')
        pyplot.close()
        data = base64.b64encode(buf.getbuffer()).decode('utf8')
        pic = f"""data:image/png;base64,{data}"""
        return pic, True, True
    chromNames = {'NC_035107.1': 'Chr1', 'NC_035108.1': 'Chr2', 'NC_035109.1': 'Chr3'}
    totalMatch = 0
    totalInv = 0
    totalLeft = 0
    totalRight = 0
    labelFontSize = 6
    axesFontSize = 8
    thickness = 12

    contigNum -= 1
    #imageNum = 0

    imageNum = next - prev



    parser = etree.XMLParser(remove_blank_text=True)
    root = etree.parse(path, parser).getroot()

    contigNames = root.xpath(f"""/root/contig[@besthit="True"][virushit[@stitle = '{virus}']]/@name""")
    node = contigNames[contigNum]

    contig = root.xpath(f"""/root/contig[@name="{node}"][@besthit="True"]""")[0]

    virusHits = contig.xpath("""./virushit""")
    features = []
    contigLength = int(contig.xpath("""./@name""")[0].split('_')[3])
    for virusHit in virusHits:
        qstart = int(virusHit.xpath("""./qstart/text()""")[0])
        qend = int(virusHit.xpath("""./qend/text()""")[0])
        sstart = int(virusHit.xpath("""./sstart/text()""")[0])
        send = int(virusHit.xpath("""./send/text()""")[0])
        evalue = virusHit.xpath("""./evalue/text()""")[0]
        bitscore = virusHit.xpath("""./bitscore/text()""")[0]
        pident = virusHit.xpath("""./pident/text()""")[0]
        stitle = virusHit.xpath("""./@stitle""")[0]
        seqid = virusHit.xpath("""./@seqid""")[0][3:-1]
        #print(seqid)
        length = abs(sstart - send) + 1
        if send > sstart:
            strandStr = '+'
        else:
            strandStr = '-'

        if qend > qstart:
            strand = 1
        else:
            strand = -1
        gf = GraphicFeature(start=qstart, end=qend, strand=strand, thickness=thickness, linewidth=0,
                                               color='#ff0000', fontdict = {'size': labelFontSize},
                                               label = stitle.lstrip('|') + ' ' + str(sstart) + '-' + str(send) + ' (' + str(length) + ' bp; ' + str(pident) + '%)')  #; ' + str(evalue) + ')')
        features.append(gf)

    caption = f"{virus} insertion sites in {spec} {specNum} ({node})"
    hitsLeft = contig.xpath('./vectorhitleft')
    hitsRight = contig.xpath('./vectorhitright')
    hitsOverlap = contig.xpath('./vectorhitoverlap')
    vqstart = int(contig.xpath('./virushit/qstart/text()')[0])
    vqend = int(contig.xpath('./virushit/qend/text()')[-1])
    flanks = contig.xpath('./flanks')[0]
    hitsDone = []

    matches = flanks.xpath('./match')
    #numMatches.append(len(matches))

    if matchNum > 0:
        totalMatch = math.ceil(len(matches)/matchNum)
        #print(totalMatch)
        if len(matches) > 0:
            matchCount = matchNum*imageNum +1
            for match in matches[imageNum*matchNum : imageNum*matchNum + matchNum]:
                for hits, attribName in [(hitsLeft, 'leftid'), (hitsRight, 'rightid')]:
                    for v in hits:
                        if v.attrib['id'] == match.attrib[attribName]:
                            hitsDone.append(v.attrib['id'])
                            qstart = int(v.find('qstart').text)
                            qend = int(v.find('qend').text)
                            if qend > qstart:
                                strand = 1
                            else:
                                strand = -1
                            sstart = int(v.find('sstart').text)
                            send = int(v.find('send').text)
                            if send > sstart:
                                strandStr = '+'
                            else:
                                strandStr = '-'
                            seqid = v.attrib['seqid']
                            if seqid in chromNames:
                                seqid = chromNames[seqid]
                            elif seqid[:8] == 'NW_01873':
                                seqid = 'Scaffold ' + seqid[8:12]
                            if attribName == 'leftid':
                                distance = vqstart - qend
                            else:
                                distance = qstart - vqend
                            features.append(GraphicFeature(start=qstart, end=qend, strand=strand, thickness=thickness, linewidth=0,
                                                       color='#009051', fontdict = {'size': labelFontSize},
                                                       label = 'M' + str(matchCount) + ': ' + seqid + ' {0:,}-{1:,} '.format(sstart, send) + '(' + strandStr + '), Dist: ' + str(distance) + ' bp'))
                            break
                matchCount += 1

    inversions = flanks.findall('inversion')
    #numInv.append(len(inversions))

    if invNum > 0:
        totalInv = math.ceil(len(inversions)/invNum)
        if len(inversions) > 0:
            matchCount = invNum*imageNum +1
            for match in inversions[invNum*imageNum : invNum*imageNum+ invNum]:
                for hits, attribName in [(hitsLeft, 'leftid'), (hitsRight, 'rightid')]:
                    for v in hits:
                        if v.attrib['id'] == match.attrib[attribName]:
                            hitsDone.append(v.attrib['id'])
                            qstart = int(v.find('qstart').text)
                            qend = int(v.find('qend').text)
                            if qend > qstart:
                                strand = 1
                            else:
                                strand = -1
                            sstart = int(v.find('sstart').text)
                            send = int(v.find('send').text)
                            if send > sstart:
                                strandStr = '+'
                            else:
                                strandStr = '-'
                            seqid = v.attrib['seqid']
                            if seqid in chromNames:
                                seqid = chromNames[seqid]
                            elif seqid[:8] == 'NW_01873':
                                seqid = 'Scaffold ' + seqid[8:12]
                            if attribName == 'leftid':
                                distance = vqstart - qend
                            else:
                                distance = qstart - vqend
                            features.append(GraphicFeature(start=qstart, end=qend, strand=strand, thickness=thickness, linewidth=0,
                                                       color='indigo', fontdict = {'size': labelFontSize},
                                                       label = 'I' + str(matchCount) + ': ' + seqid + ' {0:,}-{1:,} '.format(sstart, send) + '(' + strandStr + '), Dist: ' + str(distance) + ' bp'))
                            break
                matchCount += 1

    #numLeft.append(len(hitsLeft))
    #numRight.append(len(hitsRight))

    aaCoverage = [0] * contigLength
    for hits in (hitsLeft, hitsRight): #, hitsOverlap):
        if len(hits) > 0:
            hitsDict = {}  # consolidate query hit: [subject hits]
            for v in hits:
                pos = (v.find('qstart').text, v.find('qend').text)
                for x in range(int(pos[0]), int(pos[1])):
                    aaCoverage[x] += 1
                if v.attrib['id'] in hitsDone:
                    continue

                if pos not in hitsDict:
                    hitsDict[pos] = [(v.find('sstart').text, v.find('send').text, v.attrib['seqid'])] #, featureString)]
                else:
                    hitsDict[pos].append((v.find('sstart').text, v.find('send').text, v.attrib['seqid'])) #, featureString))
            posSort = list(hitsDict.keys())
            if hits != hitsRight:
                posSort.sort(key = lambda pos: int(pos[1]), reverse = True)  # closest first
            else:
                posSort.sort(key = lambda pos: int(pos[0]))  # closest first
            #if len(matches) == 0 and len(inversions) == 0:
            if leftFlanks > 0 or rightFlanks > 0:
                if hits == hitsLeft:
                    flankNum = leftFlanks
                    if leftFlanks == 0:
                        totalLeft = 0
                    else:
                        totalLeft = math.ceil(len(hitsLeft)/leftFlanks)
                elif hits == hitsRight:
                    flankNum = rightFlanks
                    if rightFlanks == 0:
                        totalRight = 0
                    else:
                        totalRight = math.ceil(len(hitsRight)/rightFlanks)

                for q in posSort[flankNum*imageNum : flankNum*imageNum+ flankNum]:
                    qstart = int(q[0])
                    qend = int(q[1])
                    if qend > qstart:
                        strand = 1
                    else:
                        strand = -1
                    s = hitsDict[q][0]
                    sstart = int(s[0])
                    send = int(s[1])
                    if send > sstart:
                        strandStr = '+'
                    else:
                        strandStr = '-'
                    seqid = s[2]
                    if seqid in chromNames:
                        seqid = chromNames[seqid]
                    elif seqid[:8] == 'NW_01873':
                        seqid = 'Scaffold ' + seqid[8:12]
                    if hits == hitsOverlap:
                        color = '#00ff00'
                    else:
                        color = '#0000ff'
                    if hits == hitsLeft:
                        distance = vqstart - qend
                    else:
                        distance = qstart - vqend
                    features.append(GraphicFeature(start=qstart, end=qend, strand=strand, thickness=thickness, linewidth=0,
                                               color=color, fontdict = {'size': labelFontSize},
                                               label = seqid + ' {0:,}-{1:,} '.format(sstart, send) + '(' + strandStr + '), Dist: ' + str(distance) + ' bp'))
    record = GraphicRecord(sequence_length = contigLength, features = features)

    fig, (ax1, ax2) = pyplot.subplots(2, 1, sharex = True, figsize = (10, 6), gridspec_kw = {'height_ratios': [5, 1]})

    record.plot(max_label_length = 80, ax = ax1, with_ruler = False)
    ax2.fill_between(range(contigLength), aaCoverage, step = 'mid', alpha = 0.3)
    ax2.tick_params(axis='both', which='major', labelsize=axesFontSize)
    ax2.set_ylim(bottom = 0, top = max(aaCoverage + [1]))
    ax2.set_yticks([0, max(aaCoverage + [1]) // 2, max(aaCoverage + [1])])
    ax2.set_ylabel('Aa hits', fontsize = axesFontSize)
    fig.text(.5, .01, caption, ha='center')

    buf = io.BytesIO()
    pyplot.savefig(buf, format='png')
    pyplot.close()
    data = base64.b64encode(buf.getbuffer()).decode('utf8')
    pic = f"""data:image/png;base64,{data}"""
    #totalMatches.append(max(totalMatch, totalInv, totalLeft, totalRight))

    imgmax = max(totalMatch, totalInv, totalLeft, totalRight) - 1

    if imgmax <= 0:
        out = [True, True]
    elif imageNum <= 0:
        out = [True, False]
    elif imageNum >= imgmax:
        out = [False, True]
    else:
        out = [False, False]

    return pic, out[0], out[1]

if __name__ == '__main__':
    app.run_server(debug=True)
